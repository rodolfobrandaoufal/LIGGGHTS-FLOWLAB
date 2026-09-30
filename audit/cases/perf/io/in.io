# Task 4c: output cost. vars: restart, fmt (hdf5|custom|custombin|none), nev (dump interval), nsteps
atom_style granular
atom_modify map array
boundary p p f
newton off
communicate single vel yes
units si
read_restart ${restart}
neighbor 0.0005 bin
neigh_modify delay 0
timestep 1e-5
include ../bed/in.common_props
fix integr all nve/sphere
include in.dump_${fmt}
thermo_style custom step atoms cpu
thermo ${nthermo}
run ${nsteps}
