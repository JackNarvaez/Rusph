#BSUB -J Kelvin_Helmholtz
#BSUB -W 2:00
#BSUB -n 40
#BSUB -o KelvinHelmholtz.%J
bash '../kelvin_helmholtz/kh.sh'