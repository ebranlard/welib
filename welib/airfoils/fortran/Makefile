all: prog test

prog:
	gfortran KarmanTrefftzGeometry.f90 -o KarmanTrefftzGeometry

test:
	# xc, yc, tau_deg
	echo "-0.1 0. 5" | ./KarmanTrefftzGeometry
	mv geom*.dat data/
