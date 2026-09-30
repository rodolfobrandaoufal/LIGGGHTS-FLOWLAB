#!/bin/bash
# Compare generated code of DumpCustom pack_* / count() between HEAD and the working tree
# (accessor refactor x_component()/v_component()/...). Usage: accessor_codegen_check.sh <scratch_dir>
set -e
R=$(git -C "$(dirname "$0")" rev-parse --show-toplevel); S=${1:-/tmp/accessor_check}; mkdir -p $S/head
git -C $R show HEAD:src/dump_custom.cpp > $S/head/dump_custom.cpp; git -C $R show HEAD:src/dump_custom.h > $S/head/dump_custom.h
mpicxx -O2 -std=c++11 -c $S/head/dump_custom.cpp -I$S/head -I$R/src -o $S/dc_head.o
mpicxx -O2 -std=c++11 -c $R/src/dump_custom.cpp -I$R/src -o $S/dc_new.o
echo "calls to *_component in new object (0 = fully inlined):"; objdump -d -C $S/dc_new.o | grep -E "call.*_component\(" | wc -l
for fn in pack_x pack_vx pack_xs pack_xu_triclinic pack_tqz; do
  for o in dc_head dc_new; do n=$(objdump -d --no-show-raw-insn -C $S/$o.o | awk -v f="DumpCustom::$fn(int)>:" 'index($0,f){p=1;next} p&&/^$/{exit} p{c++} END{print c+0}'); printf "%s %s %s\n" $fn $o $n; done
done
