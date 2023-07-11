#BSUB -J Sedov
#BSUB -W 1:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o Sedov.%J
bash './sedov_blast_wave/sedov.sh'
