Np=10
cargo run -p sedov_blast_wave --bin init_dist_sedov_blast_wave
for ii in $(seq 1 ${Np}); do
     d=$(RAYON_NUM_THREADS=${ii} cargo run -p speed_up);
     echo ${ii} $d
done > "./Data/speed_up/scaling.txt"
