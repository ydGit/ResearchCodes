aperture : aperture.o
	gcc -L/usr/local/lib -lgsl -lgslcblas -lm -lglib-2.0 -o aperture aperture.o

aperture.o : aperture.c 
	gcc -Wall -I/usr/local/include -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -c aperture.c