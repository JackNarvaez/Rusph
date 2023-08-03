#BSUB -J Sedov
#BSUB -W 16:00
#BSUB -n 40
#BSUB -q gpu
#BSUB -o Sedov.%J
bash './sedov_blast_wave/sedov.sh'
