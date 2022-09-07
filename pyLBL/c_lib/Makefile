libabsorption.so: absorption.o voigt.o spectra.o spectral_database.o
	$(CC) -shared -o $@ $^ -lsqlite3 -lm

%.o: %.c
	$(CC) -I${CONDA_PREFIX}/include -g -O2 -fPIC -o $@ -c $<

clean:
	rm -f libabsorption.so absorption.o voigt.o spectra.o spectral_database.o
