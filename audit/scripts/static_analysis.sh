#!/usr/bin/env bash
# clang-tidy (pip wheel, LLVM 22.1.8) + cppcheck (pip wheel, 2.17.1) over the
# modified/new CPU files, modified tree and HEAD.  Tools live in ~/.local/bin
# (installed with `pip install --user clang-tidy cppcheck`).
set -u
ROOT=/media/storage/LIGGGHTS-PUBLIC-v6
LOGS=$ROOT/audit/logs/static
mkdir -p $LOGS
CT=~/.local/bin/clang-tidy
CPPC=~/.local/bin/cppcheck
MPI_INC="-I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi"
CHECKS='-*,bugprone-*,performance-*,modernize-*,cppcoreguidelines-*,concurrency-*,mpi-*,-modernize-use-trailing-return-type,-cppcoreguidelines-avoid-magic-numbers,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-cppcoreguidelines-pro-type-vararg,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-cppcoreguidelines-avoid-c-arrays,-modernize-avoid-c-arrays,-bugprone-easily-swappable-parameters'
SCOPE_H='(contact_models|granular_styles|utils|pair_gran_base|fix_wall_gran_base|normal_model_luding|tangential_model_history|tangential_model_no_history|rolling_model_cdt|rolling_model_epsd|rolling_model_luding|cohesion_model_sjkr|cohesion_model_easo_capillary_viscous|global_properties|fix_property_global|fix_nve|velocity|neighbor|fix_neighlist_mesh|dump_custom|compute_property_atom|aligned_particle_soa|contact_model_crtp_api|cohesion_model_generalized_adhesion|fix_adapt_liggghts|dump_hdf5|dump_mesh_hdf5)\.h$'
CPPS_MOD="pair_gran_proxy fix_wall_gran global_properties fix_property_global fix_nve fix_nve_asphere_base velocity neighbor fix_neighlist_mesh dump_custom compute_property_atom fix_adapt_liggghts dump_hdf5 dump_mesh_hdf5"

mkdb() {  # $1 build dir -> writes $LOGS/<tag>_db/compile_commands.json with g++ + MPI includes
  local b=$1 tag=$2
  mkdir -p $LOGS/${tag}_db
  python3 - "$b/compile_commands.json" "$LOGS/${tag}_db/compile_commands.json" "$MPI_INC" <<'EOF'
import json, sys
src, dst, mpi = sys.argv[1:4]
cc = json.load(open(src))
for e in cc:
    e['command'] = e['command'].replace('/usr/bin/mpicxx', '/usr/bin/g++ ' + mpi, 1) \
                               .replace('-Wno-literal-suffix', '')
json.dump(cc, open(dst, 'w'), indent=1)
EOF
}

run_tidy() {  # $1 tag, $2 srcdir
  local tag=$1 src=$2
  for f in $CPPS_MOD; do [ -f $src/$f.cpp ] && echo $src/$f.cpp; done | \
    xargs -P 16 -I{} sh -c "$CT -p $LOGS/${tag}_db --quiet --extra-arg=--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/11 --checks='$CHECKS' --header-filter='$SCOPE_H' {} > $LOGS/tidy_${tag}_\$(basename {} .cpp).txt 2>&1"
  cat $LOGS/tidy_${tag}_*.txt | grep -E '^/.*: (warning|error): ' | sort -u > $LOGS/tidy_${tag}_all.txt
}

mkdb $ROOT/build_audit/release modified
mkdb $ROOT/build_audit/baseline baseline
run_tidy modified $ROOT/build_audit/src_modified
run_tidy baseline $ROOT/build_audit/baseline_src/src

# cppcheck: only the scope files, with project includes
for tag in modified baseline; do
  if [ $tag = modified ]; then src=$ROOT/build_audit/src_modified; defs="-DLIGGGHTS_HDF5 -DLAMMPS_SMALLBIG"; else src=$ROOT/build_audit/baseline_src/src; defs="-DLAMMPS_SMALLBIG"; fi
  files=$(for f in $CPPS_MOD; do [ -f $src/$f.cpp ] && echo $src/$f.cpp; done)
  $CPPC --enable=all --inconclusive --std=c++17 --language=c++ -j 16 $defs -I$src \
        --suppress=missingIncludeSystem --suppress=missingInclude --suppress=unusedFunction \
        --template='{file}:{line}: {severity}: {id}: {message}' $files \
        > $LOGS/cppcheck_${tag}_stdout.txt 2> $LOGS/cppcheck_${tag}_raw.txt
  names="$(echo $CPPS_MOD | sed 's/ /|/g')|$(echo "$SCOPE_H" | sed 's/^(//; s/)\\.h\$$//')"
  grep -E "/($names)\.(h|cpp):" $LOGS/cppcheck_${tag}_raw.txt | sort -u > $LOGS/cppcheck_${tag}_scope.txt
done
echo done
