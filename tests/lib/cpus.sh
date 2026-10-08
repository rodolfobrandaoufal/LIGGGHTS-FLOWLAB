# Shared helper for the test runners (LIGGGHTS modernization branch, R2-T04).
# liggghts_cpus <default list>: the taskset CPU list a suite uses when neither
# its own variable nor LIGGGHTS_TEST_CPUS is set. The default is kept when
# taskset accepts it on this machine; otherwise (e.g. "14-27" on a 4-CPU CI
# runner, where taskset fails with EINVAL) the CPUs this process may run on.
liggghts_cpus() {
  if taskset -c "$1" true 2>/dev/null; then
    echo "$1"
  else
    taskset -cp $$ | sed 's/.*: *//'
  fi
}
