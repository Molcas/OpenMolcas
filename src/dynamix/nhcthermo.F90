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
use Dynamix_Globals, only: dyn_nh
#endif
use Dynamix_Globals, only: DT, TEMP, nh, iQ1, iQ2, iX1, iX2, iVx1, iVx2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: kBoltzmann, auTokJ, Zero, Two, Three, Half, Quart
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vel(*)
integer(kind=iwp) :: i, j, natom
real(kind=wp) :: Ekin, DT2, DT4, DT8, G1, G2, NHC(nh), Q1, Q2, sc, Vx1, Vx2, X1, X2
real(kind=wp), allocatable :: Mass(:)
real(kind=wp), parameter :: kb = kBoltzmann/(auTokJ*1.0e3_wp)

!nh stands for the number of variables in the thermostat
! NHC = Q1,Q2,X1,X2,VX1,VX2,Scale
!
! The parameter kb is the Boltzmann constant in E_h/K

! READ PARAMETERS FROM RUNFILE

call Get_nAtoms_Full(natom)

call mma_allocate(Mass,natom)

! Read Thermostat Variables
call Get_NHC(NHC,nh)

Q1 = NHC(iQ1)
Q2 = NHC(iQ2)
X1 = NHC(iX1)
X2 = NHC(iX2)
Vx1 = NHC(iVx1)
Vx2 = NHC(iVx2)

! Initialize the Mass variable
call GetMassDx(Mass,natom)

Ekin = Zero

do i=1,natom
  do j=1,3
    Ekin = Ekin+Half*Mass(i)*(vel(3*(i-1)+j)**2)
  end do
end do

DT2 = DT*Half
DT4 = DT*Quart
DT8 = DT*0.125_wp

G2 = (Q1*Vx1*Vx1-TEMP*kb)/Q2
Vx2 = Vx2+G2*DT4
Vx1 = Vx1*exp(-Vx2*DT8)
G1 = (Two*Ekin-Three*natom*TEMP*kb)/Q1
Vx1 = Vx1+G1*DT4
Vx1 = Vx1*exp(-Vx2*DT8)
sc = exp(-Vx1*DT2)

vel(1:3*natom) = vel(1:3*natom)*sc
Ekin = Ekin*sc*sc

X1 = X1+Vx1*DT2
X2 = X2+Vx2*DT2
Vx1 = Vx1*exp(-Vx2*DT8)
G1 = (Two*Ekin-Three*natom*TEMP*kb)/Q1
Vx1 = Vx1+G1*DT4
Vx1 = Vx1*exp(-Vx2*DT8)
G2 = (Q1*(Vx1**2)-TEMP*kb)/Q2
Vx2 = Vx2+G2*DT4

NHC(iX1) = X1
NHC(iX2) = X2
NHC(iVx1) = Vx1
NHC(iVx2) = Vx2

call Put_NHC(NHC,nh)
#ifdef _HDF5_
call mh5_put_dset(dyn_nh,NHC)
#endif

call mma_deallocate(Mass)

return

end subroutine NhcThermo
