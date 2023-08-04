Np=40
for ii in $(seq 1 ${Np}); do
     d=$(RAYON_NUM_THREADS=${ii} cargo run -p speed_up);
     echo ${ii} $d
done > "./Data/speed_up/scaling.txt"
