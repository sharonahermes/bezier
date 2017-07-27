all: src/bezier/_speedup/libspeedup.so

src/bezier/_speedup/libspeedup.so: src/bezier/_speedup/speedup.f90
	gfortran -shared -fPIC \
	  src/bezier/_speedup/speedup.f90 \
	  -o src/bezier/_speedup/libspeedup.so
	rm speedup.mod

clean:
	rm -f src/bezier/_speedup/libspeedup.so

.PHONY: all clean
