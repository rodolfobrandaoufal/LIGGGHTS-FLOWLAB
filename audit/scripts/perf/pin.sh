#!/bin/bash
# Pin each MPI rank to one hardware thread in 16..31 (rank r -> cpu 16+r mod 16).
# Usage: mpirun -np N --bind-to none pin.sh <binary> <args...>
r=${OMPI_COMM_WORLD_RANK:-0}
b=${PIN_BASE:-16}
exec taskset -c $((b + r % (32 - b))) "$@"
