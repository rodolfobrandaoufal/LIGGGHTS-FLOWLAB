#!/bin/bash
# round-2 paired A/B campaign (idle machine); results in audit2/data/perf/ab/
cd /media/storage/LIGGGHTS-PUBLIC-v6/audit2/scripts/perf
OUT=/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/perf/ab; mkdir -p $OUT
for e in flat25k flat200k sbf25k sbf200k sbf200k_np4 pred_normal pred_full sync twist lastforce incremental newton_np4 omp_build_1t omp_det hdf5 dumpcustom lammps25k lammps200k_np4 multi_poly_np1 multi_poly_np4; do
  python3 ab.py $e 6 $OUT > $OUT/$e.stdout 2>&1
done
for e in b01_idle b01_nobal; do python3 ab.py $e 8 $OUT > $OUT/$e.stdout 2>&1; done
echo DONE > $OUT/campaign.done
