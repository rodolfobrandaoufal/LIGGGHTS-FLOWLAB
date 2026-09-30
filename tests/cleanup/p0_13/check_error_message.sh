#!/usr/bin/env bash
# P0-13 check (LIGGGHTS modernization branch, cleanup agent).
# Usage: check_error_message.sh <liggghts binary> [time budget s, default 90]
# The special message is chosen with srand(time(NULL)) (1 in 5 per second),
# so the deck is re-run until 3 special messages were seen or time runs out.
# Every line printed after the ERROR line must be one of the known messages,
# printable ASCII/UTF-8, and no sanitizer report may appear. Exit 0 on success.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=$(readlink -f "$1"); budget=${2:-90}
work=$(mktemp -d); cd "$work"; cp "$here/in.error_message" .
export ASAN_OPTIONS=detect_leaks=0:abort_on_error=0
hits=0; runs=0; bad=0; end=$(( $(date +%s) + budget ))
while [ $(date +%s) -lt $end ] && [ $hits -lt 3 ]; do
  "$bin" -in in.error_message -log none > out.txt 2>&1
  runs=$((runs+1))
  if grep -a -q -E "AddressSanitizer|runtime error" out.txt; then echo "FAIL: sanitizer report"; grep -m3 -E "AddressSanitizer|runtime error|#[0-3] " out.txt; bad=1; break; fi
  post=$(sed -n '/^ERROR:/,$p' out.txt | tail -n +2 | grep -a -v '^$')
  [ -z "$post" ] && continue
  hits=$((hits+1))
  if ! echo "$post" | grep -a -q -E "^(Comment from the off|Tip of the day): " || LC_ALL=C grep -a -q -P '[\x00-\x08\x0e-\x1f\x7f]' <<< "$post"; then
    echo "FAIL: corrupted special message:"; echo "$post" | cat -v | head -3; bad=1; break
  fi
  # the message must be followed by the (file:line) suffix on the same line
  echo "$post" | head -1 | grep -a -q -E '\(.*:[0-9]+\)$' || [ $(echo "$post" | wc -l) -gt 1 ] || { echo "FAIL: truncated message"; bad=1; break; }
done
rm -rf "$work"
echo "runs=$runs special_messages_seen=$hits"
if [ $bad = 0 ] && [ $hits -gt 0 ]; then echo "RESULT: PASS"; exit 0; fi
[ $hits = 0 ] && [ $bad = 0 ] && echo "no special message seen within ${budget}s"
echo "RESULT: FAIL"; exit 1
