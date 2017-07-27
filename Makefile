all: src/bezier/libspeedup.so

src/bezier/libspeedup.so: src/bezier/speedup.f90
	gfortran -shared -fPIC \
	  src/bezier/speedup.f90 \
	  -o src/bezier/libspeedup.so

clean:
	rm -f src/bezier/libspeedup.so speedup.mod

.PHONY: all clean
