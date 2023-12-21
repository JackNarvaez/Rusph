all: Toystar

Sedov: ./tests/sedov_blast_wave/sedov.sh
	@mkdir $@
	@bash $<

Toystar: ./tests/toy_star/toy_star.sh
	@mkdir $@
	@bash $<

Kelvinhelmholtz: ./tests/kelvin_helmholtz/kh.sh
	@mkdir $@
	@bash $<

Turbulence: ./tests/turbulent_gas/turbulence.sh
	@mkdir $@
	@bash $<

clean:
	@rm -rf target;\
	rm -f *.lock
