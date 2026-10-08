/* Minimal flat sampling profiler (perf is unavailable: perf_event_paranoid=4, no root).
 * LD_PRELOAD=sprof.so SPROF_OUT=prefix SPROF_HZ=2000 ./lmp ...
 * Samples the interrupted RIP on SIGPROF (ITIMER_PROF = process CPU time) and, at exit,
 * writes "<offset_in_main_binary_hex> <count>" lines (PIE offset, usable with addr2line / nm).
 * Samples outside the main binary (libc, libmpi, libhdf5) are aggregated by mapping name. */
#define _GNU_SOURCE
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <ucontext.h>
#include <unistd.h>
#define TSZ (1<<22)
static unsigned long *keys; static unsigned *cnt; static unsigned long total, dropped;
static void handler(int sig, siginfo_t *si, void *uc_) {
  ucontext_t *uc = (ucontext_t*)uc_;
  unsigned long pc = (unsigned long)uc->uc_mcontext.gregs[REG_RIP];
  unsigned long h = (pc * 0x9E3779B97F4A7C15UL) >> 42;
  for (int k = 0; k < 64; k++) {
    unsigned long s = (h + k) & (TSZ-1);
    if (keys[s] == pc) { cnt[s]++; total++; return; }
    if (keys[s] == 0) { keys[s] = pc; cnt[s] = 1; total++; return; }
  }
  dropped++;
}
__attribute__((constructor)) static void init(void) {
  if (!getenv("SPROF_OUT")) return;
  /* ITIMER_PROF survives execve(): only arm the timer in the target binary (taskset/stdbuf exec chain) */
  { char e[4096]; ssize_t k = readlink("/proc/self/exe", e, sizeof e - 1); e[k > 0 ? k : 0] = 0;
    const char *m = getenv("SPROF_MATCH") ? getenv("SPROF_MATCH") : "lmp_"; if (!strstr(e, m)) return; }
  keys = calloc(TSZ, sizeof(*keys)); cnt = calloc(TSZ, sizeof(*cnt));
  struct sigaction sa; memset(&sa, 0, sizeof sa);
  sa.sa_sigaction = handler; sa.sa_flags = SA_SIGINFO | SA_RESTART;
  sigaction(SIGPROF, &sa, NULL);
  int hz = getenv("SPROF_HZ") ? atoi(getenv("SPROF_HZ")) : 1000;
  struct itimerval it; it.it_interval.tv_sec = 0; it.it_interval.tv_usec = 1000000 / hz;
  it.it_value = it.it_interval; setitimer(ITIMER_PROF, &it, NULL);
}
__attribute__((destructor)) static void fini(void) {
  if (!keys) return;
  struct itimerval it; memset(&it, 0, sizeof it); setitimer(ITIMER_PROF, &it, NULL);
  char exe[4096]; ssize_t n = readlink("/proc/self/exe", exe, sizeof exe - 1); exe[n > 0 ? n : 0] = 0;
  /* read maps */
  FILE *m = fopen("/proc/self/maps", "r"); char line[8192];
  unsigned long lo[4096], hi[4096], off[4096]; char *name[4096]; int nm = 0;
  while (m && fgets(line, sizeof line, m) && nm < 4096) {
    unsigned long a, b, o; char perm[8], path[4096] = "";
    if (sscanf(line, "%lx-%lx %7s %lx %*s %*s %4095s", &a, &b, perm, &o, path) >= 4) {
      lo[nm] = a; hi[nm] = b; off[nm] = o; name[nm] = strdup(path); nm++; }
  }
  if (m) fclose(m);
  char fn[4200]; const char *r = getenv("OMPI_COMM_WORLD_RANK");
  snprintf(fn, sizeof fn, "%s.%s.%d.txt", getenv("SPROF_OUT"), r ? r : "0", getpid());
  FILE *f = fopen(fn, "w"); if (!f) return;
  fprintf(f, "# total %lu dropped %lu exe %s\n", total, dropped, exe);
  for (unsigned long s = 0; s < TSZ; s++) if (keys[s]) {
    unsigned long pc = keys[s]; int j;
    for (j = 0; j < nm; j++) if (pc >= lo[j] && pc < hi[j]) break;
    if (j < nm && strcmp(name[j], exe) == 0) fprintf(f, "M %lx %u\n", pc - lo[j] + off[j], cnt[s]);
    else fprintf(f, "L %s %u\n", j < nm && name[j][0] ? name[j] : "?", cnt[s]);
  }
  fclose(f);
}
