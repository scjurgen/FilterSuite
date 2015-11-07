all:		    mkfilter mkshape mkaverage gencode genplot

mkfilter:   mkfilter.o complex.o
		    g++ mkfilter.o complex.o -lm
		    mv a.out mkfilter

mkshape:    mkshape.o complex.o
		    g++ mkshape.o complex.o -lm
		    mv a.out mkshape

mkaverage:  mkaverage.o complex.o
		    g++ mkaverage.o complex.o -lm
		    mv a.out mkaverage

gencode:    gencode.o complex.o readdata.o
		    g++ gencode.o complex.o readdata.o -lm
		    mv a.out gencode

genplot:    genplot.o complex.o readdata.o
		    g++ genplot.o complex.o readdata.o -L/usr/local/lib -lgd -lm
		    mv a.out genplot

%.o:		    %.c mkfilter.h complex.h
		    g++ -I/usr/local/lib -O -c $<

clean:
		    rm -f mkfilter.o mkshape.o mkaverage.o gencode.o genplot.o complex.o readdata.o gencode mkfilter mkaverage mkshape genplot

