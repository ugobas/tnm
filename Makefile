#       Makefile
#
#       Copyright 2009 Unknown <fons@arnold.cbm.uam.es>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#        gaussj0.o \
#        Weighted_flux.o \
#        EC.o \


####################
# Progam files
COBJ= main_tnm.o \
	read.o \
	atoms.o \
	seqres.o \
	vector.o \
	allocate_tnm.o \
	align_tnm.o \
	dof_tnm.o \
	buildup.o \
	interactions_tnm.o \
	anharmonic_tnm.o \
	contacts.o \
	screened_interactions.o \
	kinetic_tnm.o \
	nma.o \
	Fit_B_factors.o \
	ridge_regression.o \
	confchange.o \
	simulation.o \
	force2confchange.o \
	mutation.o \
	Fit.o \
	optimization.o \
	output_tnm.o \
	allostery.o \
	site_dynamics.o \
	allocate.o \
	choldc.o ludcmp.o lubksb.o \
	d_diagonalize.o \
	nrutil.o \
	pythag.o \
	tred2.o \
	jacobi.o \
	tqli.o \
	Profit_aux.o \
	NeedlemanWunsch.o \
	McLachlan_float.o \
	McLachlan.o \
	EC.o \
	gasdev.o \
	ran2.o \
	random3.o 
# shadow_interactions.o 
#	unfolding.o \

PROG=tnm

####################
# Compiler flags
CC=gcc
FC=f77
FFLAG=-g -ffixed-line-length-132
#CFLAGS= -Wall -std=c99 -pedantic -g -pg -fbounds-check -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
CFLAGS= -g -O2 -Wall -std=c99 -pedantic -pg -fbounds-check -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
# segmentation fault with -O3 !!
LDFLAGS=-lm -lg2c
LDFLAGS=-lm -L/usr/lib64/libg2c.so.0.0.0


all: $(FOBJ) $(COBJ)
	$(CC) $(FOBJ) $(COBJ) $(F77LIB) -o $(PROG) $(LDFLAGS)

intel: CC=icc
intel: FC=ifort
intel: LDFLAGS=-lm
intel: CFLAGS=-Wall -static -O3 -xHOST -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -vec-report
intel: all

intel-static: CC=icc
intel-static: FC=ifort
intel-static: LDFLAGS=-lm -static
intel-static: CFLAGS=-Wall -static -O3 -xHOST -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -vec-report
intel-static: all

gnu: CC=g77
gnu: FC=f77
gnu: LDFLAGS=-lm -lg2c
gnu: CFLAGS=-Wall -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
gnu: FFLAGS= -O3 -ffixed-line-length-132
gnu: all

gnu-static: CC=g77
gnu-static: FC=f77
gnu-static: LDFLAGS=-lm -lg2c -static
gnu-static: CFLAGS=-Wall -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
gnu-static: FFLAGS= -O3 -ffixed-line-length-132
gnu-static: all


%.o: %.cpp
	$(CC) $(LDFLAGS) $(CFLAGS) -c $< -o $@

%.o: %.f
	$(FC)  $(FFLAG) -c $< -o $@

clean:
	rm -fr $(COBJ) $(FOBJ) $(PROG)
