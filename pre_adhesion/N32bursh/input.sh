for file in *.data
do
    name=${file%.data}
    mpirun -np 16 lmp_mpi -in ./in.PDMS -var name $name
done
