#!/bin/bash


if [ -z $SUPERPOSE_PATH ]; then
    export SUPERPOSE_PATH="/home/sulcjo/Downloads/cath-superpose.ubuntu-20.04"
fi

echo "
########################
Welcome to Cartesian Allostery Tool. =^.^=/
Will look for trajectory .pdbs in this directory.
Those will be aligned using SSAP algorithm and prepared for analysis.

Make sure you're SUPERPOSE_PATH environment variable is set
/// export SUPERPOSE_PATH=/PATH/TO/SUPERPOSE/BINARY
########################
"

read -p 'Trajectory 1 name: name1.pdb ... ' traj1
read -p 'Trajectory 2 name: name2.pdb ... ' traj2
read -p 'Max concurent alignment processes to run: <integer>, more means more resources used, default 30 ' max_concurrent_procs

rm -r $traj1 ${traj1}_ssap $traj2 ${traj2}_ssap &> /dev/null



# Prepare, separate trajectories into single .pdbs
mkdir $traj1 $traj2 ${traj1}_ssap ${traj2}_ssap temp

# Speed up trajectory splitting with asynchronicity, continue when the loop is done

echo 0 | gmx trjconv -f $traj1.pdb -s $traj1.pdb -o $traj1/$traj1.pdb -sep
echo 0 | gmx trjconv -f $traj2.pdb -s $traj2.pdb -o $traj2/$traj2.pdb -sep

run_loop() {
        # These loops run asynchronously and sometimes throw out errors for some iterations
        # use && for each command so they don't create empty frames
        header=$(head -5 $1/${1}${2}.pdb)
        echo "$header" > ${1}_ssap/#header_temp_${2}


        mkdir temp/${2} # TEMP files have to be separated so it runs correctly
        $SUPERPOSE_PATH --pdb-infile $CATH_TOOLS_PDB_PATH/REFERENCE.pdb --pdb-infile $CATH_TOOLS_PDB_PATH/${1}${2}.pdb --do-the-ssaps temp/${2} --sup-to-pdb-files-dir ${1}_ssap &>> ${1}_ssap.log --align-refining heavy &&

        # Create temporary .pdb file (aligned), modify
        cat ${1}_ssap/#header_temp_${2} ${1}_ssap/${1}${2}.pdb > ${1}_ssap/#${1}${2}.pdb &&
        echo "ENDMDL" >> ${1}_ssap/#${1}${2}.pdb &&
        # Replace the original with the temporary modified one
        mv ${1}_ssap/#${1}${2}.pdb ${1}_ssap/${1}${2}.pdb
}


ssap_trajectory() {
    # $1 is the name for the trajectory file

    # Set environment for SSAP
    export CATH_TOOLS_PDB_PATH=$PWD/$1
    num_pdb=$(ls -1 $CATH_TOOLS_PDB_PATH/*.pdb | wc -l)

    rm ${1}_ssap.log &> /dev/null
    rm ${1}_ssap_traj.pdb &> /dev/null
    # Start alignment procedure for the trajectory (align all to 0)
    #sp="/-\|" # Set-up progress spiner
    #echo -n ' '


    #chmod 444 $CATH_TOOLS_PDB_PATH/REFERENCE.pdb #ssap tries to overwrite this, which is bad for asynchronous
    rm -r temp
    mkdir temp
    for i in $(seq 0 $(($num_pdb-1))); do
        # Run alignment asynchronously in background
        run_loop $1 $i &
        # But only allow max_concurrent_procs processes to run at a single time
        joblist=($(jobs -p))
        while (( ${#joblist[*]} >= $max_concurrent_procs )); do
            sleep 0.1
            joblist=($(jobs -p))
        done
    done
    #chmod 777 $CATH_TOOLS_PDB_PATH/REFERENCE.pdb

    # Check if all loops finished, then continue
    while [ $(ls -1 ${1}/*.pdb | wc -l) -ne $num_pdb ]; do
        sleep 0.1
    done

    # Combine into a contiguous trajectory .pdb file
    echo "Constructing contiguous trajectory"
    sp="/-\|" # Set-up progress spiner
    echo -n ' '
    for i in $(seq 0 $num_pdb); do
        sed -i '/END/d' ${1}_ssap/${1}${i}.pdb
        echo "ENDMDL" >> ${1}_ssap/${1}${i}.pdb
        cat ${1}_ssap/${1}${i}.pdb >> ${1}_ssap_traj.pdb
        printf "\b${sp:i++%${#sp}:1}" # Spin the progress bar
    done

    echo "END" >> ${1}_ssap_traj.pdb
    rm ${1}_ssap/REFERENCE.pdb
    rm -r temp
}

echo "Will now align ${traj1} to its reference structure 0 using SSAP algorithm"
cp ${traj1}/${traj1}0.pdb ${traj1}/REFERENCE.pdb
ssap_trajectory ${traj1}

echo "Will now align ${traj2} to ${traj1} reference structure 0 using SSAP algorithm"
cp ${traj1}/${traj1}0.pdb ${traj2}/REFERENCE.pdb
ssap_trajectory ${traj2}

echo "
########################
Trajectory alignment complete, check both
in your favorite visualization software.
########################
"
