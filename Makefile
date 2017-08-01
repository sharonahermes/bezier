all: src/bezier/_speedup/speedup.so

src/bezier/_speedup/speedup.so: src/bezier/_speedup/speedup.o src/bezier/_speedup/speedup.pyx
	python make_cython_so.py build_ext --inplace
	mv speedup.so src/bezier/_speedup/
	rm -fr build/ src/bezier/_speedup/speedup.c \
	  src/bezier/_speedup/speedup.mod src/bezier/_speedup/speedup.o

src/bezier/_speedup/speedup.o: src/bezier/_speedup/speedup.f90
	gfortran -shared -fPIC -c src/bezier/_speedup/speedup.f90
	mv speedup.o speedup.mod src/bezier/_speedup/

clean:
	rm -f src/bezier/_speedup/speedup.so

.PHONY: all clean
