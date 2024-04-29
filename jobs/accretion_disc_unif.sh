#BSUB -J AccDiscUnif
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o AccDiscUnif.%J
make Accretiondiscuniform
