all: src/bezier/_speedup.so

src/bezier/_speedup.so: src/bezier/_speedup.pyf src/bezier/speedup.f90 .f2py_f2cmap
	f2py --verbose -c --opt='-O3' \
	    src/bezier/_speedup.pyf src/bezier/speedup.f90
	mv _speedup*.so src/bezier

brute-force-pyf: src/bezier/speedup.f90 .f2py_f2cmap
	f2py src/bezier/speedup.f90 -m _speedup \
	  -h src/bezier/_speedup.pyf --overwrite-signature

clean:
	rm -f src/bezier/_speedup*.so

.PHONY: all clean brute-force-pyf
