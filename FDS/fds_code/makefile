CFLAGS=-Ofast
CFLAGS=-O0
# $(CFLAGS)
all:
	gfortran -fopenmp -shared $(CFLAGS) -g -o fds.so prec.f90 mpis.f90 cons.f90 devc.f90 type.f90 mesh.f90 func.f90 data.f90 smvv.f90 irad.f90 turb.f90 soot.f90 \
     ieva.f90 pois.f90 scrc.f90 radi.f90 evac.f90 gsmv.f90 geom.f90 part.f90 vege.f90 ctrl.f90 samr.f90 dump.f90 hvac.f90 mass.f90 \
     read.f90 wall.f90 fire.f90 divg.f90 velo.f90 pres.f90 init.f90 main.f90 -fPIC
