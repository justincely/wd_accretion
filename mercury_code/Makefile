alll: compile run convert close

compile:
	g77 -o mercury6 mercury6_3.for
	g77 -o element6 element6.for
	g77 -o close6 close6.for

osx:
	/sw/bin/gfortran -o mercury6 mercury6_3.for
	/sw/bin/gfortran -o element6 element6.for
	/sw/bin/gfortran -o close6 close6.for


run:
	./mercury6

convert:
	./element6

close:
	./close6
