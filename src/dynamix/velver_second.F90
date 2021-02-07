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
! Copyright (C) 2006, Igor Schapiro                                    *
!***********************************************************************
!
! *********************************************************************
! *                                                                   *
! * Second part of the velocity Verlet algorithm, which calculates    *
! * velocities of the next time step from the velocities from the     *
! * half step and the new forces.                                     *
! *                                                                   *
! * The algorithm is based on the example F4 of molecular dynamics    *
! * bible from Allen and Tildesley.                                   *
! *                                                                   *
! * 15/10/2006                                                        *
! * Igor Schapiro                                                     *
! *                                                                   *
! *********************************************************************

! REFERENCE:
! SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.

subroutine VelVer_Second(irc)

#ifdef _HDF5_
use mh5, only: mh5_put_dset
#endif

implicit real*8(a-h,o-z)
#include "prgm.fh"
#include "warnings.fh"
#include "Molcas.fh"
parameter(ROUTINE='VV_Second')
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
external IsFreeUnit
integer i, j, natom, file, IsFreeUnit, nsAtom
character filname*80, line*80
real*8 conv, tolerance, DT2, time, kb
real*8 Ekin, Epot, Etot, Etot0, Ekin_target
real*8 EtotLastPoint
logical hybrid, found
character caption*15, lastline*80
character, allocatable :: atom(:)*2
real*8, allocatable :: Mass(:), vel(:), force(:), xyz(:)

! The parameter conv converts the gradients (Hartree/Bohr) to
! (i)  forces (Hartree/Bohr)    => -1.0d0
! (ii) forces (kJ/mole/Agstrom) => -4961.475514610d0

parameter(conv=-1.0d0)
!parameter(conv=-CONV_AU_TO_KJ_PER_MOLE_/Angstrom)
parameter(kb=CONST_BOLTZMANN_/(CONV_AU_TO_KJ_*1.0d3))

if (IPRINT == INSANE) write(6,*) ' Entering ',ROUTINE

write(6,*) '*** Second step of the Velocity Verlet algorithm ***'

filname = 'comqum.dat'
call F_INQUIRE(filname,hybrid)

if (hybrid) then
  write(6,'(/,5X,A)') 'Perform QM/MM Molecular Dynamics'
  call DxRdNAtomHbrd(natom)
else
  call DxRdNAtomStnd(natom)
end if

call mma_allocate(atom,natom)
call mma_allocate(Mass,natom)
call mma_allocate(vel,natom*3)
call mma_allocate(force,natom*3)
call mma_allocate(xyz,natom*3)

if (hybrid) then
  call DxRdHbrd(natom,atom,xyz,force)
else
  call DxRdStnd(natom,atom,xyz,force)
end if

call Get_Velocity(vel,3*natom)
call GetMassDx(Mass,natom)

! Check if reduced dimensionality
if (POUT /= 0) then
  call project_out_for(force,natom)
elseif (PIN /= natom*3) then
  call project_in_for(force,natom)
end if

! Definition of the time step

DT2 = DT/2.0d0
Ekin = 0.0d0

call Get_dScalar('MD_Time',time)

do i=1,natom
  do j=1,3
    vel(3*(i-1)+j) = vel(3*(i-1)+j)+DT2*force(3*(i-1)+j)/Mass(i)
  end do
end do
!  Calling the thermostats for canonical ensemble
if (THERMO == 2) then
  call NhcThermo(vel)
end if

! Check if reduced dimensionality (should not be needed)
if (POUT /= 0) then
  call project_out_vel(vel,natom)
end if

! Final kinetic energy
do i=1,natom
  do j=1,3
    Ekin = Ekin+5.0d-01*Mass(i)*(vel(3*(i-1)+j)**2)
  end do
end do

call Add_Info('EKin',[EKin],1,6)

! Write out the velocities

caption = 'Velocities'
lastline = ''
call DxPtTableCo(caption,time,natom,atom,vel,lastline,Mass,force)

! Read in the potential energy

if (hybrid) then
  file = IsFreeUnit(81)
  filname = 'fixenergy.out'
  call Molcas_Open(file,filname)
  read(file,*) line
  do while (i <= 80)
    if (line(i:i) == '$') then
      read(file,'(E20.13)') Epot
      i = 80
    elseif (i == 80 .and. line(i:i) /= '$') then
      write(6,*) 'No energy found'
      call Abend()
    end if
    i = i+1
  end do
  close(file)
else
  call Get_dScalar('Last Energy',Epot)
end if

Etot = Epot+Ekin

! Output

if (hybrid) then
  write(6,'(//,5X,A,F8.1)') 'Final QM/MM Energy at time ',time
else
  write(6,'(//,5X,A,F8.1)') 'Final Energy at time ',time
end if
write(6,'(5X,A,/)') '============================'
write(6,'(5X,A,6X,D19.12,1X,A)') 'Kinetic energy',Ekin,'a.u.'
write(6,'(5X,A,4X,D19.12,1X,A)') 'Potential Energy',Epot,'a.u.'
write(6,'(5X,A,8X,D19.12,1X,A)') 'Total Energy',Etot,'a.u.'

!--------------------------------------------------------------------C
! CANONICAL ENSEMBLE
!--------------------------------------------------------------------C
tempNow = 2.d0*Ekin/(3.d0*dble(natom)*kb)
if (THERMO == 2) then
  tempNow = 2.d0*Ekin/(3.d0*dble(natom)*kb)
  write(6,'(//,5X,A)') 'Canonical Ensemble'
  write(6,'(5X,A)') 'The temperature is control with a '
  write(6,'(5X,A)') 'Nose-Hoover chain of thermostats'
  write(6,'(5X,A,/)') '========================'
  write(6,'(5X,A,5X,E11.4,1X,A,//)') 'instantaneous temperature',tempNow,'kelvin'

end if
!--------------------------------------------------------------------C

!---------------------------------------------------------------------C
!     MICRO-CANONICAL ENSEMBLE
!---------------------------------------------------------------------C
if (THERMO == 1) then
  call Get_dScalar('MD_Etot0',Etot0)
  write(6,'(//,5X,A)') 'Micro-Canonical Ensemble'
  write(6,'(5X,A,/)') '========================'
  write(6,'(5X,A,7X,D11.4,1X,A)') 'Target Total Energy',Etot0,'a.u.'
  write(6,'(5X,A,D11.4,1X,A)') 'Deviation from this Energy',abs(Etot0-Etot),'a.u.'

  ! Check if the total energy is conserved and scale the velocities
  ! if necessary.

  ! 1.0K * k_B
  tolerance = 1.0d0*kb
  tolerance = 1.5d0*natom*tolerance
  if (abs(Etot0-Etot) > tolerance) then
    Ekin_target = Etot0-Epot
    scalfac = sqrt(Ekin_target/Ekin)
    do i=1,natom
      do j=1,3
        vel(3*(i-1)+j) = scalfac*vel(3*(i-1)+j)
      end do
    end do
    write(6,'(5X,A)') 'is larger then Scaling-Threshold XX.'
    write(6,'(5X,A)') 'Velocity scaling is necessary.'
    write(6,401) 'Velocity scaling factor',scalfac
    Ekin = Ekin_target
  else
    write(6,'(5X,A)') 'is smaller then Scaling-Threshold XX.'
    write(6,'(5X,A)') 'Velocity scaling is not necessary.'
  end if
end if
!---------------------------------------------------------------------C
call Put_dScalar('MD_Time',time)

Etot = Epot+Ekin

!------ Tully reascaling velocities in case of HOP ---------------------C
call qpg_iscalar('hopped',Found)
if (found) call get_lscalar('hopped',Found)
if (found) then
  call Get_iScalar('Unique atoms',nsAtom)
  call get_dScalar('MD_Etot',EtotLastPoint)
  Ekin_target = EtotLastPoint-Epot
  if (Ekin_target < 0.0) then
    Ekin_target = 0.0
    write(6,*) 'warning negative kin energy rescaled to 0.0'
  end if
  !write(6,*) 'DEBUGS'
  !write(6,*) Epot,EtotLastPoint,Etot
  !write(6,*) 'EtotLastPoint'
  !write(6,*) EtotLastPoint
  !write(6,*) 'scalfac'
  !write(6,*) scalfac
  !write(6,*) 'nsAtom'
  !write(6,*) nsAtom
  scalfac = sqrt(Ekin_target/Ekin)
  write(6,*) 'Velocities before Hop:'
  do i=1,nsAtom
    write(6,*) vel(i*3-2),vel(i*3-1),vel(i*3)
  end do
  do i=1,natom
    do j=1,3
      vel(3*(i-1)+j) = scalfac*vel(3*(i-1)+j)
    end do
  end do

  write(6,*) 'Velocities after Hop:'
  do i=1,nsAtom
    write(6,*) vel(i*3-2),vel(i*3-1),vel(i*3)
  end do
  Etot = Epot+Ekin_target
end if

call Put_dScalar('MD_Etot',Etot)
#ifdef _HDF5_
call mh5_put_dset(dyn_etot,Etot)
#endif
call DxEnergies(time,Epot,Ekin,Etot)

call DxWtVel(vel,3*natom)
call Put_Velocity(vel,3*natom)
#ifdef _HDF5_
call mh5_put_dset(dyn_vel,vel)
#endif

call mma_deallocate(atom)
call mma_deallocate(Mass)
call mma_deallocate(vel)
call mma_deallocate(force)
call mma_deallocate(xyz)

! The return code is set in order to continue the loop

irc = _RC_ALL_IS_WELL_

return

401 format(5x,a,3x,e11.4)

end subroutine VelVer_Second
