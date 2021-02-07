!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2012, Felipe Zapata                                    *
!***********************************************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! The Nose-Hoover chain of thermostat is based on the following paper
! Journal of Physical Chemistry B, 2001, 105, 7598
! The implemetation required first initialized the thermostat, then
! before calling the first part of velocity verlet the NHC is call,
! then the positions and velocities are updated using the verlet_first
! and after that the verlet_second is called and after this step
! the NHC subroutine is call. The order must not be changed
!
! Felipe Zapata, November 1, 2012
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine NhcThermo(vel)

#ifdef _HDF5_
use mh5, only: mh5_put_dset
#endif

implicit real*8(a-h,o-z)
#include "prgm.fh"
#include "Molcas.fh"
parameter(ROUTINE='NhcThermo')
#include "warnings.fh"
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
parameter(nh=6)
integer natom, i, j
real*8 Ekin, kb
parameter(kb=CONST_BOLTZMANN_/(CONV_AU_TO_KJ_*1.0d3))
real*8 NHC(nh), vel(*)
real*8, allocatable :: Mass(:)

!nh stands for the number of variables in the thermostat
! NHC = Q1,Q2,X1,X2,VX1,VX2,Scale
!
! The parameter kb is the Boltzmann constant in E_h/K

! READ PARAMETERS FROM RUNFILE

call Get_nAtoms_Full(natom)

call mma_allocate(Mass,natom)

! Read Thermostat Variables
call Get_NHC(NHC,nh)

Q1 = NHC(1)
Q2 = NHC(2)
X1 = NHC(3)
X2 = NHC(4)
Vx1 = NHC(5)
Vx2 = NHC(6)

! Initialize the Mass variable
call GetMassDx(Mass,natom)

Ekin = 0.0d0

do i=1,natom
  do j=1,3
    Ekin = Ekin+5.0D-01*Mass(i)*(vel(3*(i-1)+j)**2.d0)
  end do
end do

DT2 = DT*0.5d0
DT4 = DT*2.5D-1
DT8 = DT*1.25D-1

G2 = (Q1*Vx1*Vx1-TEMP*kb)/Q2
Vx2 = Vx2+G2*DT4
Vx1 = Vx1*exp(-Vx2*DT8)
G1 = (2.d0*Ekin-3.d0*dble(natom)*TEMP*kb)/Q1
Vx1 = Vx1+G1*DT4
Vx1 = Vx1*exp(-Vx2*DT8)
sc = exp(-Vx1*DT2)

do i=1,natom
  do j=1,3
    vel(3*(i-1)+j) = vel(3*(i-1)+j)*sc
  end do
end do

Ekin = Ekin*sc*sc

X1 = X1+Vx1*DT2
X2 = X2+Vx2*DT2
Vx1 = Vx1*exp(-Vx2*DT8)
G1 = (2.d0*Ekin-3.d0*dble(natom)*TEMP*kb)/Q1
Vx1 = Vx1+G1*DT4
Vx1 = Vx1*exp(-Vx2*DT8)
G2 = (Q1*(Vx1**2.d0)-TEMP*kb)/Q2
Vx2 = Vx2+G2*DT4

NHC(3) = X1
NHC(4) = X2
NHC(5) = Vx1
NHC(6) = Vx2

call Put_NHC(NHC,nh)
#ifdef _HDF5_
call mh5_put_dset(dyn_nh,NHC)
#endif

call mma_deallocate(Mass)

return

end subroutine NhcThermo
