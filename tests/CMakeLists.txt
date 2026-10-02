# ----------------------------------------------------------------------
#    CTest registration of the LIGGGHTS regression suites (roadmap A1).
#    Included from src/CMakeLists.txt (LIGGGHTS_ENABLE_TESTING=ON, default).
#
#    Contributing author: LIGGGHTS modernization branch
# ----------------------------------------------------------------------
#
# Every test runs one of the existing driver scripts in this directory
# against the liggghts_bin target of the current build.
#
# Labels:
#   unit        fast, single-feature checks (seconds)
#   physics     checks against analytic or cross-run expectations
#   mpi         runs on more than one MPI rank
#   regression  guards a fixed bug or the tutorial decks
#   bitwise     byte-identity against LIGGGHTS_TEST_REFERENCE_BIN
#   io          dump hdf5 / dump mesh/hdf5 (skipped without HDF5 or h5py)
#
# Optional binaries (cache variables); a test that needs one is SKIPPED
# (return code 77) when it is not available:
#   LIGGGHTS_TEST_REFERENCE_BIN  reference build for the bitwise tests
#   LIGGGHTS_TEST_SQ_BIN         superquadric build (default: this build if ENABLE_SQ=ON)
#   LIGGGHTS_TEST_STRICT_BIN     separate build with LIGGGHTS_NO_CONTACT_MODEL_FALLBACK=ON
#   LIGGGHTS_TEST_ASAN_BIN       ASan+UBSan build (default: this build if LIGGGHTS_SANITIZE is set)
#   LIGGGHTS_TEST_CPUS           taskset CPU list used inside the scripts (default: 0-15)
#
# Run:  ctest --output-on-failure                (all)
#       ctest -L 'regression|physics|mpi' -j4    (CI selection)
#       ctest -LE bitwise                        (without the reference tests)

SET(LIGGGHTS_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
SET(LIGGGHTS_TEST_WORK "${CMAKE_CURRENT_BINARY_DIR}/work")
FILE(MAKE_DIRECTORY "${LIGGGHTS_TEST_WORK}")

SET(LIGGGHTS_TEST_REFERENCE_BIN "" CACHE FILEPATH "Reference liggghts binary for the bitwise-identity tests (empty: tests skipped)")
SET(LIGGGHTS_TEST_SQ_BIN "" CACHE FILEPATH "Superquadric (ENABLE_SQ) liggghts binary for the superquadric tests (empty: this build if ENABLE_SQ=ON)")
SET(LIGGGHTS_TEST_STRICT_BIN "" CACHE FILEPATH "Separate liggghts binary built with LIGGGHTS_NO_CONTACT_MODEL_FALLBACK=ON (empty: test skipped)")
SET(LIGGGHTS_TEST_ASAN_BIN "" CACHE FILEPATH "ASan+UBSan liggghts binary (empty: this build if LIGGGHTS_SANITIZE is set)")
SET(LIGGGHTS_TEST_CPUS "" CACHE STRING "CPU list for taskset inside the test scripts (empty: scripts' default 0-15)")

FIND_PROGRAM(LIGGGHTS_BASH bash)
FIND_PROGRAM(LIGGGHTS_PYTHON3 python3)
IF(NOT LIGGGHTS_BASH OR NOT LIGGGHTS_PYTHON3)
  MESSAGE(WARNING "bash or python3 not found: LIGGGHTS tests are not registered")
  RETURN()
ENDIF()

SET(BIN "$<TARGET_FILE:liggghts_bin>")

SET(SQ_BIN "${LIGGGHTS_TEST_SQ_BIN}")
IF(NOT SQ_BIN AND ENABLE_SQ)
  SET(SQ_BIN "${BIN}")
ENDIF()
# strict mode is compared against this (fallback-enabled) build, so it is
# always a separate binary; a strict main build fails dispatch_suite by design
SET(STRICT_BIN "${LIGGGHTS_TEST_STRICT_BIN}")
SET(ASAN_BIN "${LIGGGHTS_TEST_ASAN_BIN}")
IF(NOT ASAN_BIN AND LIGGGHTS_SANITIZE)
  SET(ASAN_BIN "${BIN}")
ENDIF()
SET(REF_BIN "${LIGGGHTS_TEST_REFERENCE_BIN}")

SET(test_env "ASAN_OPTIONS=detect_leaks=0")
IF(LIGGGHTS_TEST_CPUS)
  LIST(APPEND test_env "LIGGGHTS_TEST_CPUS=${LIGGGHTS_TEST_CPUS}")
ENDIF()

# liggghts_test(<name> LABELS <l...> TIMEOUT <s> [PROCESSORS <n>] [ENV <k=v...>] COMMAND <cmd...>)
FUNCTION(LIGGGHTS_TEST name)
  CMAKE_PARSE_ARGUMENTS(T "" "TIMEOUT;PROCESSORS" "LABELS;ENV;COMMAND" ${ARGN})
  ADD_TEST(NAME ${name} COMMAND ${T_COMMAND})
  SET(env ${test_env} ${T_ENV})
  IF(NOT T_PROCESSORS)
    SET(T_PROCESSORS 1)
  ENDIF()
  SET_TESTS_PROPERTIES(${name} PROPERTIES
    LABELS "${T_LABELS}" TIMEOUT ${T_TIMEOUT} PROCESSORS ${T_PROCESSORS}
    ENVIRONMENT "${env}" SKIP_RETURN_CODE 77
    WORKING_DIRECTORY "${LIGGGHTS_TEST_WORK}")
ENDFUNCTION()

# a test that cannot run in this configuration: reported as "Skipped" with the reason
FUNCTION(LIGGGHTS_SKIPPED_TEST name reason)
  ADD_TEST(NAME ${name} COMMAND ${LIGGGHTS_BASH} -c "echo 'SKIPPED: ${reason}'; exit 77")
  SET_TESTS_PROPERTIES(${name} PROPERTIES LABELS "${ARGN}" SKIP_RETURN_CODE 77)
ENDFUNCTION()

SET(no_ref  "no reference binary; configure with -DLIGGGHTS_TEST_REFERENCE_BIN=<liggghts>")
SET(no_sq   "no superquadric binary; build with -DENABLE_SQ=ON or set -DLIGGGHTS_TEST_SQ_BIN=<liggghts>")

#---------------------------------------------------------------------------
# A2: v_ properties reach the force law (props)
LIGGGHTS_TEST(props_suite LABELS physics mpi TIMEOUT 900 PROCESSORS 2
  ENV "WORKDIR=${LIGGGHTS_TEST_WORK}/props"
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/props/run_all.sh ${BIN})

#---------------------------------------------------------------------------
# A3: fix adapt/liggghts (adapt); superquadric cases need an SQ binary
IF(SQ_BIN)
  SET(adapt_sq "${SQ_BIN}")
ELSE()
  SET(adapt_sq "-")
ENDIF()
LIGGGHTS_TEST(adapt_suite LABELS physics mpi TIMEOUT 1500 PROCESSORS 8
  ENV "WD=${LIGGGHTS_TEST_WORK}/adapt"
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/adapt/run_all.sh ${BIN} ${adapt_sq} -)

#---------------------------------------------------------------------------
# A4/A5: dispatch, fallback, cohesion compatibility, whitelist coverage
# (the optional reference and strict binaries are registered as separate tests below)
LIGGGHTS_TEST(dispatch_suite LABELS regression mpi TIMEOUT 900 PROCESSORS 2
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/dispatch/run_all.sh ${BIN})
IF(STRICT_BIN)
  LIGGGHTS_TEST(dispatch_strict_errors LABELS regression TIMEOUT 600
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/dispatch/check_strict_errors.sh ${STRICT_BIN} ${BIN} ${LIGGGHTS_TEST_WORK}/dispatch_strict)
ELSE()
  LIGGGHTS_SKIPPED_TEST(dispatch_strict_errors
    "no strict binary; set -DLIGGGHTS_TEST_STRICT_BIN=<liggghts built with LIGGGHTS_NO_CONTACT_MODEL_FALLBACK=ON>" regression)
ENDIF()

#---------------------------------------------------------------------------
# cross-fix integration (A2 refresh x A5 property names)
LIGGGHTS_TEST(integration_suite LABELS physics TIMEOUT 600
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/integration/run_all.sh ${BIN})

#---------------------------------------------------------------------------
# A6/A7 (cleanup): split into single checks so that optional parts skip cleanly
LIGGGHTS_TEST(cleanup_a6_removed LABELS unit regression TIMEOUT 60
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/cleanup/check_a6_removed.sh ${LIGGGHTS_SRC_DIR})
LIGGGHTS_TEST(cleanup_rolling_luding LABELS unit physics TIMEOUT 600
  COMMAND ${LIGGGHTS_PYTHON3} ${LIGGGHTS_TESTS_DIR}/cleanup/rolling_luding/check_rolling_luding.py
          ${BIN} ${LIGGGHTS_TEST_WORK}/rolling_luding)
LIGGGHTS_TEST(cleanup_p0_13_error_message LABELS unit regression TIMEOUT 300
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/cleanup/p0_13/check_error_message.sh ${BIN} 120)
LIGGGHTS_TEST(cleanup_p0_12_startup LABELS unit regression mpi TIMEOUT 300 PROCESSORS 2
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/cleanup/p0_12/check_startup_sanitizer.sh ${BIN})
IF(SQ_BIN)
  LIGGGHTS_TEST(cleanup_f18_sq LABELS regression TIMEOUT 600
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/cleanup/f18/check_f18.sh ${SQ_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(cleanup_f18_sq "${no_sq}" regression)
ENDIF()
IF(ASAN_BIN)
  LIGGGHTS_TEST(cleanup_p0_12_startup_asan LABELS unit regression mpi TIMEOUT 600 PROCESSORS 2
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/cleanup/p0_12/check_startup_sanitizer.sh ${ASAN_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(cleanup_p0_12_startup_asan
    "no ASan binary; build with -DLIGGGHTS_SANITIZE=address,undefined or set -DLIGGGHTS_TEST_ASAN_BIN=<liggghts>" unit regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# A8: dump hdf5 / dump mesh/hdf5 (hdf5 agent's runner; it exits 77 = skipped
# when h5py is missing or the binary has no HDF5 support)
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/hdf5/run_all.sh")
  LIGGGHTS_TEST(hdf5_suite LABELS io mpi TIMEOUT 900 PROCESSORS 4
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/hdf5/run_all.sh ${BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(hdf5_suite "tests/hdf5/run_all.sh not present" io mpi)
ENDIF()

#---------------------------------------------------------------------------
# A9 hygiene (gitignore, tutorials, tan_luding guard); bitwise parts only with a reference
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/hygiene/run_all.sh")
  LIGGGHTS_TEST(hygiene_suite LABELS regression TIMEOUT 1200 PROCESSORS 2
    ENV "WORK=${LIGGGHTS_TEST_WORK}/hygiene"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/hygiene/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(hygiene_suite "tests/hygiene/run_all.sh not present" regression)
ENDIF()

#---------------------------------------------------------------------------
# Phase B physics options (roadmap B1-B4, B6, B7). Each suite verifies its opt-in
# model against a stated reference; with a reference binary it also checks that
# the default input is bitwise unchanged.
FOREACH(suite tangential adhesion easo_dt normal wear quick)
  IF(EXISTS "${LIGGGHTS_TESTS_DIR}/${suite}/run_all.sh")
    LIGGGHTS_TEST(${suite}_suite LABELS physics mpi TIMEOUT 1500 PROCESSORS 2
      ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
      COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/${suite}/run_all.sh ${BIN} ${REF_BIN})
  ELSE()
    LIGGGHTS_SKIPPED_TEST(${suite}_suite "tests/${suite}/run_all.sh not present" physics mpi)
  ENDIF()
ENDFOREACH()

#---------------------------------------------------------------------------
# C1 load balancing (balance / fix balance): correctness on 1-16 ranks incl. mesh
# walls and contact history. BALANCE_QUICK=1 restricts it to 4 ranks.
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/balance/run_all.sh")
  LIGGGHTS_TEST(balance_suite LABELS physics mpi TIMEOUT 2400 PROCESSORS 16
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/balance/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(balance_suite "tests/balance/run_all.sh not present" physics mpi)
ENDIF()

#---------------------------------------------------------------------------
# C4 size-class neighbor multi for granular lists: nsq/bin/multi give identical
# contacts and history on 1-4 ranks, incl. mesh walls, growth and insertion
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/neigh/run_all.sh")
  LIGGGHTS_TEST(neigh_suite LABELS regression mpi TIMEOUT 900 PROCESSORS 4
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/neigh/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(neigh_suite "tests/neigh/run_all.sh not present" regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# C2/B8 OpenMP threading: deterministic mode bitwise identical across thread
# counts and to the serial kernel; exits 77 on a build without OpenMP
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/omp/run_all.sh")
  LIGGGHTS_TEST(omp_suite LABELS physics regression TIMEOUT 1800 PROCESSORS 8
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/omp/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(omp_suite "tests/omp/run_all.sh not present" physics regression)
ENDIF()

#---------------------------------------------------------------------------
# K-01/K-02: edinburgh (FMA NaN in wall contacts, unbounded loop) and
# thornton_ning (overshoot past pull-off) stop cleanly or run correctly
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/legacy/run_all.sh")
  LIGGGHTS_TEST(legacy_suite LABELS physics regression mpi TIMEOUT 1200 PROCESSORS 2
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/legacy/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(legacy_suite "tests/legacy/run_all.sh not present" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# C3 newton pair on with contact history: single contacts bitwise equal to newton
# off; history survives reneighbouring, migration, restart; 1-8 ranks
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/newton/run_all.sh")
  LIGGGHTS_TEST(newton_suite LABELS physics mpi TIMEOUT 1200 PROCESSORS 8
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/newton/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(newton_suite "tests/newton/run_all.sh not present" physics mpi)
ENDIF()

#---------------------------------------------------------------------------
# X-02/X-03 restart continuity (time, insertion RNG) and delete_atoms with
# contact history (compress renumbering, stale history copies)
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/restart/run_all.sh")
  LIGGGHTS_TEST(restart_suite LABELS physics regression mpi TIMEOUT 1200 PROCESSORS 8
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/restart/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(restart_suite "tests/restart/run_all.sh not present" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# X-04 pair-flip history signs and degenerate-input (FMA) guards
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/signfma/run_all.sh")
  LIGGGHTS_TEST(signfma_suite LABELS physics regression mpi TIMEOUT 900 PROCESSORS 2
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/signfma/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(signfma_suite "tests/signfma/run_all.sh not present" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# X-01 mesh walls in periodic boxes: rank-independent results, single coplanar
# edge contact (coplanar_legacy restores the old double count)
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/meshpbc/run_all.sh")
  LIGGGHTS_TEST(meshpbc_suite LABELS physics regression mpi TIMEOUT 1500 PROCESSORS 9
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/meshpbc/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(meshpbc_suite "tests/meshpbc/run_all.sh not present" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# phase E leftovers: superquadric acos and wall heat-conduction guards,
# insert/stream restart continuity, fp-contract option
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/misc/run_all.sh")
  LIGGGHTS_TEST(misc_suite LABELS physics regression mpi TIMEOUT 1200 PROCESSORS 4
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/misc/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(misc_suite "tests/misc/run_all.sh not present" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# X-05 pair history cleared at separation (default == rebuild every step) and
# multicontact radius of i; legacy keywords reproduce the earlier behaviour
IF(EXISTS "${LIGGGHTS_TESTS_DIR}/x05/run_all.sh")
  LIGGGHTS_TEST(x05_suite LABELS physics regression mpi TIMEOUT 1500 PROCESSORS 2
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/x05/run_all.sh ${BIN} ${REF_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(x05_suite "tests/x05/run_all.sh not present" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# superquadric a1/a2 pair-flip swap, compute pair/gran/local side effect, newton on
# (needs a superquadric build: this build with ENABLE_SQ or LIGGGHTS_TEST_SQ_BIN)
IF(SQ_BIN AND EXISTS "${LIGGGHTS_TESTS_DIR}/sq/run_all.sh")
  LIGGGHTS_TEST(sq_suite LABELS physics regression mpi TIMEOUT 1800 PROCESSORS 2
    ENV "TMPDIR=${LIGGGHTS_TEST_WORK}"
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/sq/run_all.sh ${SQ_BIN})
ELSE()
  LIGGGHTS_SKIPPED_TEST(sq_suite "no superquadric binary (ENABLE_SQ or LIGGGHTS_TEST_SQ_BIN)" physics regression mpi)
ENDIF()

#---------------------------------------------------------------------------
# B5 model matrix: every static whitelist combination, byte comparison against
# the reference binary (needs LIGGGHTS_TEST_REFERENCE_BIN)
IF(REF_BIN AND EXISTS "${LIGGGHTS_TESTS_DIR}/kernel/run_all.sh")
  LIGGGHTS_TEST(kernel_model_matrix LABELS regression bitwise TIMEOUT 1800 PROCESSORS 4
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/kernel/run_all.sh ${BIN} ${REF_BIN} ${LIGGGHTS_TEST_WORK}/kernel)
ELSE()
  LIGGGHTS_SKIPPED_TEST(kernel_model_matrix "no reference binary (LIGGGHTS_TEST_REFERENCE_BIN)" regression bitwise)
ENDIF()

#---------------------------------------------------------------------------
# whitelist entry LUDING/TANGENTIAL_HISTORY/ROLLING_LUDING runs on the static path
LIGGGHTS_TEST(whitelist_luding_rolling_luding LABELS unit regression TIMEOUT 120
  COMMAND ${BIN} -in ${LIGGGHTS_TESTS_DIR}/build/in.luding_rolling_luding -log none)
SET_TESTS_PROPERTIES(whitelist_luding_rolling_luding PROPERTIES
  PASS_REGULAR_EXPRESSION "Loop time of [^\n]* for 2000 steps"
  FAIL_REGULAR_EXPRESSION "is not in the static contact-model whitelist;ERROR")

#---------------------------------------------------------------------------
# tests/regression: negative deck, must stop with the controlled error (P0-04, F-19)
IF(SQ_BIN)
  LIGGGHTS_TEST(regression_asphere_scheme4_requires_implicit LABELS regression TIMEOUT 120
    COMMAND ${SQ_BIN} -in ${LIGGGHTS_TESTS_DIR}/regression/in.asphere_scheme4_requires_implicit -log none)
  SET_TESTS_PROPERTIES(regression_asphere_scheme4_requires_implicit PROPERTIES
    PASS_REGULAR_EXPRESSION "ERROR: integration_scheme 4 requires fix couple/cfd/force/implicit \\(")
ELSE()
  LIGGGHTS_SKIPPED_TEST(regression_asphere_scheme4_requires_implicit "${no_sq}" regression)
ENDIF()

#---------------------------------------------------------------------------
# tutorial smoke: every Tutorials_public deck for 10 steps (skips SQ/HDF5 decks
# the binary cannot run; known pre-existing failures are XFAIL)
LIGGGHTS_TEST(tutorials_smoke LABELS regression TIMEOUT 900 PROCESSORS 4
  ENV "LIGGGHTS_TEST_JOBS=4"
  COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/tutorials/run_tutorials.sh ${BIN} 10 120 ${LIGGGHTS_TEST_WORK}/tutorials)

#---------------------------------------------------------------------------
# bitwise identity against a reference binary (chute_wear np1/np2 dumps, packing thermo)
# cleanup/bitwise carries its own decks; the dispatch and adapt variants read
# the audit cases (audit/cases/npdep, audit/cases/smoke)
GET_FILENAME_COMPONENT(LIGGGHTS_AUDIT_CASES "${LIGGGHTS_TESTS_DIR}/../audit/cases" ABSOLUTE)
IF(REF_BIN)
  LIGGGHTS_TEST(bitwise_cleanup LABELS regression bitwise mpi TIMEOUT 1200 PROCESSORS 2
    COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/cleanup/bitwise/check_bitwise.sh ${BIN} ${REF_BIN} ${LIGGGHTS_TEST_WORK}/bitwise_cleanup)
  IF(EXISTS "${LIGGGHTS_AUDIT_CASES}/npdep/lmp_release_np1/in.chute_wear" AND EXISTS "${LIGGGHTS_AUDIT_CASES}/smoke/lmp_release__packing/in.packing")
    LIGGGHTS_TEST(bitwise_dispatch LABELS regression bitwise mpi TIMEOUT 1200 PROCESSORS 2
      COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/dispatch/check_bitwise.sh ${REF_BIN} ${BIN} ${LIGGGHTS_TEST_WORK}/bitwise_dispatch)
    LIGGGHTS_TEST(bitwise_adapt_identity LABELS regression bitwise mpi TIMEOUT 1200 PROCESSORS 2
      COMMAND ${LIGGGHTS_BASH} ${LIGGGHTS_TESTS_DIR}/adapt/check_identity.sh ${REF_BIN} ${BIN} ${LIGGGHTS_TEST_WORK}/bitwise_adapt)
  ELSE()
    LIGGGHTS_SKIPPED_TEST(bitwise_dispatch "${LIGGGHTS_AUDIT_CASES} (npdep, smoke) not present" regression bitwise mpi)
    LIGGGHTS_SKIPPED_TEST(bitwise_adapt_identity "${LIGGGHTS_AUDIT_CASES} (npdep, smoke) not present" regression bitwise mpi)
  ENDIF()
ELSE()
  LIGGGHTS_SKIPPED_TEST(bitwise_dispatch "${no_ref}" regression bitwise mpi)
  LIGGGHTS_SKIPPED_TEST(bitwise_cleanup "${no_ref}" regression bitwise mpi)
  LIGGGHTS_SKIPPED_TEST(bitwise_adapt_identity "${no_ref}" regression bitwise mpi)
ENDIF()
