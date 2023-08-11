all: Toystar

Sedov: ./tests/sedov_blast_wave/sedov.sh ./tests/sedov_blast_wave/input ./tests/sedov_blast_wave/Cargo.toml 
	@bash $<

Toystar: ./tests/toy_star/toy_star.sh
	@bash $<

Kelvinhelmholtz: ./tests/kelvin_helmholtz/kh.sh
	@bash $<

Turbulence: ./tests/turbulent_gas/turbulence.sh
	@bash $<

clean:
	@rm -rf target;\
	rm -f *.lock
