# Compile C library and Cpython wrapper
# This will work for Linux but it's not tested in other platforms.

# cLoadNeuralynx.pyx loadNeuralynx.o 
all: cLoadNeuralynx.pyx loadNeuralynx.o 
	python setup.py build_ext --inplace

loadNeuralynx.o: loadNeuralynx.c
	gcc -c -O2 -Wall loadNeuralynx.c


#cython -a cloaddata.pyx
#gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.6 -o cloaddata.so cloaddata.c
