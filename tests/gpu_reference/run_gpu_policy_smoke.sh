#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LIGGGHTS_BINARY="${LIGGGHTS_BINARY:-${ROOT_DIR}/src/lmp_hdf5mpi}"
WORK_DIR="${TMPDIR:-/tmp}/liggghts_gpu_policy_smoke.$$"

cleanup() {
  rm -rf "${WORK_DIR}"
}
trap cleanup EXIT

mkdir -p "${WORK_DIR}"

require_binary() {
  if [[ ! -x "${LIGGGHTS_BINARY}" ]]; then
    echo "LIGGGHTS binary is not executable: ${LIGGGHTS_BINARY}" >&2
    exit 2
  fi
}

write_base_input() {
  local path="$1"
  local extra_commands="$2"

  cat > "${path}" <<EOF
units si
atom_style atomic
boundary f f f
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 single 0.5 0.5 0.5
mass 1 1.0
pair_style soft 1.0
pair_coeff * * 1.0
fix int all nve
${extra_commands}
run 0
EOF
}

assert_log_contains() {
  local log="$1"
  local expected="$2"

  if ! grep -Fq "${expected}" "${log}"; then
    echo "Missing expected log text:" >&2
    echo "  ${expected}" >&2
    echo "Log file: ${log}" >&2
    exit 1
  fi
}

run_liggghts() {
  local input="$1"

  rm -f "${WORK_DIR}/log.liggghts"
  (cd "${WORK_DIR}" && "${LIGGGHTS_BINARY}" -in "${input}")
}

print_log_on_failure() {
  local log="${WORK_DIR}/log.liggghts"

  if [[ -f "${log}" ]]; then
    echo "---- ${log} ----" >&2
    tail -n 80 "${log}" >&2
    echo "----------------" >&2
  fi
}

run_success_case() {
  local name="$1"
  local commands="$2"
  shift 2
  local input="${WORK_DIR}/${name}.in"
  local log="${WORK_DIR}/log.liggghts"

  write_base_input "${input}" "${commands}"
  if ! run_liggghts "${input}"; then
    echo "Expected success case failed: ${name}" >&2
    print_log_on_failure
    exit 1
  fi

  for expected in "$@"; do
    assert_log_contains "${log}" "${expected}"
  done
}

run_failure_case() {
  local name="$1"
  local commands="$2"
  local expected="$3"
  local input="${WORK_DIR}/${name}.in"
  local log="${WORK_DIR}/log.liggghts"

  write_base_input "${input}" "${commands}"
  if run_liggghts "${input}"; then
    echo "Expected failure case succeeded unexpectedly: ${name}" >&2
    exit 1
  fi

  assert_log_contains "${log}" "${expected}"
}

require_binary

run_success_case "auto_defaults" \
  "gpu_mode auto" \
  "GPU_DEM mode set to auto; strict rejects host fallbacks, auto reports them and runs on CPU" \
  "GPU_DEM feature probe: mode=auto precision=mixed neighbor=auto device=auto streams=on" \
  "GPU_DEM unsupported: only atom_style sphere/granular is supported" \
  "GPU_DEM auto mode: falling back to CPU for this run"

run_success_case "explicit_policy" \
  "gpu_device 0
gpu_precision mixed
gpu_neighbor verlet
gpu_streams off
gpu_mode auto" \
  "GPU_DEM device policy set to 0" \
  "GPU_DEM precision policy set to mixed" \
  "GPU_DEM neighbor policy set to verlet" \
  "GPU_DEM streams policy set to off" \
  "GPU_DEM feature probe: mode=auto precision=mixed neighbor=verlet device=0 streams=off" \
  "GPU_DEM auto mode: falling back to CPU for this run"

run_failure_case "strict_rejects_hard_blocker" \
  "gpu_mode strict" \
  "ERROR: gpu_mode strict rejected this run; see GPU_DEM unsupported messages above"

run_failure_case "invalid_streams" \
  "gpu_streams maybe
gpu_mode auto" \
  "ERROR: gpu_streams expects 'on' or 'off'"

echo "GPU_DEM policy smoke tests passed using ${LIGGGHTS_BINARY}"
