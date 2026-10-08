#!/bin/bash
# line/branch coverage of src/ by the CTest suite (gcc --coverage -O0, gcovr)
W=/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work
D=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/coverage; mkdir -p $D
# (started after the perf campaign)
cd $W/clone/src && cmake -B $W/builds/cov -DLIGGGHTS_TEST_REFERENCE_BIN= -DLIGGGHTS_TEST_SQ_BIN= > /dev/null
cd $W/builds/cov && find . -name "*.gcda" -delete && ctest -j4 --timeout 7200 > $D/ctest_cov.log 2>&1
PYTHONPATH=$W/pylib python3 -m gcovr -r $W/clone/src $W/builds/cov --gcov-ignore-errors=all --gcov-ignore-parse-errors=all --txt-summary --json-summary $D/summary.json --txt $D/coverage_by_file.txt > $D/gcovr_summary.txt 2>&1
echo DONE > $D/coverage.done
