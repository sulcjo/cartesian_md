#!/bin/bash

output_path=/run/timeshift/backup/IOCB/MSM/fitted_trajectories

<<COMMENT
# Terminal resis have different contexts (connected, free with -/+), they have to be cut

for i in run_{20..29}; do
gmx make_ndx -f ${i}/md.tpr -o ${i}/index.ndx<<EOF
r2-81
name 19 pdz
19&2
q
EOF
done

# ^ for FD3, 608 atoms

for i in run_{20..29}; do
gmx make_ndx -f ${i}/md.tpr -o ${i}/index.ndx<<EOF
r2-81
name 19 pdz
19&2
q
EOF
done

# ^ for PDZ3, 608 atoms

for i in run_{20..29}; do
gmx make_ndx -f ${i}/md.tpr -o ${i}/index.ndx<<EOF
r39-118
name 19 pdz
19&2
q
EOF
done

# ^ for FD4, 608 atoms

COMMENT


#
##
### Throw away first 10000 ps of each
### Prepare index.ndx with pdz group in advance for fd3 and fd4

mkdir $output_path/pdz
mkdir $output_path/fd3
mkdir $output_path/fd4

PDZ_path=/run/timeshift/backup/IOCB/MSM/pdz/charmm_chwater
    # run_1...run_19
FD3_path=/run/timeshift/backup/IOCB/MSM/pdz_l_trp
    # run_1...run_19
FD4_path=/run/timeshift/backup/IOCB/MSM/trp_l_pdz_closed
    # run_1...run_19


mkdir $output_path/pdz/run_1
# Prepare global reference trajectory

gmx trjconv -f $PDZ_path/run_1/md_noPBC.xtc -s $PDZ_path/run_1/md.tpr -o $output_path/pdz/run_1/fit.pdb -fit rot+trans -b 10000 -n $PDZ_path/run_1/index.ndx <<EOF
4
pdz_&_Protein-H
EOF
echo 0 | gmx traj -f $output_path/pdz/run_1/fit.pdb -s $output_path/pdz/run_1/fit.pdb -ox $output_path/pdz/run_1.xvg

# Fit all PDZ3 trajectories to global reference
for i in run_{2..19}; do
    mkdir $output_path/pdz/$i
    echo "pdz_&_Protein-H" | gmx trjconv -f $PDZ_path/${i}/md_noPBC.xtc -s $PDZ_path/${i}/md.tpr -o $output_path/pdz/${i}/pdz_only.pdb -n $PDZ_path/${i}/index.ndx -b 10000

    gmx trjconv -f $output_path/pdz/${i}/pdz_only.pdb -s $output_path/pdz/run_1/fit.pdb -o $output_path/pdz/${i}/fit.pdb -fit rot+trans -b 10000 <<EOF
    4
    0
EOF

echo 0 | gmx traj -f $output_path/pdz/${i}/fit.pdb -s $output_path/pdz/${i}/fit.pdb -ox $output_path/pdz/${i}.xvg

done

# Fit all FD3 trajectories to global reference
for i in run_{1..19}; do
    mkdir $output_path/fd3/$i

    echo "pdz_&_Protein-H" | gmx trjconv -f $FD3_path/${i}/md_noPBC.xtc -s $FD3_path/${i}/md.tpr -o $output_path/fd3/${i}/pdz_only.pdb -n $FD3_path/${i}/index.ndx -b 10000
    gmx trjconv -f $output_path/fd3/${i}/pdz_only.pdb -s $output_path/pdz/run_1/fit.pdb -o $output_path/fd3/${i}/fit.pdb -fit rot+trans -b 10000 <<EOF
    4
    0
EOF
echo 0 | gmx traj -f $output_path/fd3/${i}/fit.pdb -s $output_path/fd3/${i}/fit.pdb -ox $output_path/fd3/${i}.xvg

done

# Fit all FD4 trajectories to global reference
for i in run_{1..19}; do
    mkdir $output_path/fd4/$i

    echo "pdz_&_Protein-H" | gmx trjconv -f $FD4_path/${i}/md_noPBC.xtc -s $FD4_path/${i}/md.tpr -o $output_path/fd4/${i}/pdz_only.pdb -n $FD4_path/${i}/index.ndx -b 10000
    gmx trjconv -f $output_path/fd4/${i}/pdz_only.pdb -s $output_path/pdz/run_1/fit.pdb -o $output_path/fd4/${i}/fit.pdb -fit rot+trans -b 10000 <<EOF
    4
    0
EOF
echo 0 | gmx traj -f $output_path/fd4/${i}/fit.pdb -s $output_path/fd4/${i}/fit.pdb -ox $output_path/fd4/${i}.xvg
done



# Assign vectors for cartesian
#######
for variant in pdz fd3 fd4; do
    cd $variant
    for i in run_{1..19}; do
        cartesian.py --f ${i}.xvg --o ${i}_${variant}_vectors.json --resi ${i}/fit.pdb
    done
cd ..
done
