#BSUB -J Kelvin_Helmholtz
#BSUB -W 16:00
#BSUB -n 40
#BSUB -q gpu
#BSUB -o KelvinHelmholtz.%J
bash './kelvin_helmholtz/kh.sh'
