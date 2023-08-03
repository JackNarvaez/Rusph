#BSUB -J Turbulence
#BSUB -W 16:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o Turbulence.%J
bash './turbulent_gas/turbulence.sh'