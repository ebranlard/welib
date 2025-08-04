all:
	gfortran -c PrecisionMod.f90
	gfortran -c FortLibTest.f90
	gfortran -c SurfaceVorticity2DTools.f90
	gfortran test_cylinder2d_nolift.f90 PrecisionMod.o FortLibTest.o SurfaceVorticity2DTools.o -o cylinder2d_nolift.exe
	./cylinder2d_nolift.exe


all2:
	gfortran -c PrecisionMod.f90
	gfortran -c FortLibTest.f90
	gfortran -c SurfaceVorticity2DTools.f90
	gfortran test_ellipse_2d_lift.f90 PrecisionMod.o FortLibTest.o SurfaceVorticity2DTools.o -o ellipse_2d_lift.exe
	./ellipse_2d_lift.exe



# all3:
# 	gfortran -c PrecisionMod.f90
# 	gfortran -c FortLibTest.f90
# 	gfortran -c SurfaceVorticity2DTools.f90
# 	gfortran test_ellipse_2d_no_lift.f90 PrecisionMod.o FortLibTest.o SurfaceVorticity2DTools.o -o ellipse_2d_no_lift.exe
# 	./ellipse_2d_no_lift.exe
# 
