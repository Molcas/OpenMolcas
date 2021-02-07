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
! * First part of the velocity Verlet algorithm, which calculates the *
! * new positions for the next time step and the new velocities for   *
! * a half time step. The algorithm is based on the example F4 from   *
! * molecular dynamics bible of Allen and Tildesley.                  *
! *                                                                   *
! * 15/10/2006                                                        *
! * Igor Schapiro                                                     *
! *                                                                   *
! *********************************************************************

subroutine VelVer_First(irc)

#ifdef _HDF5_
use mh5, only: mh5_put_dset
#endif

implicit real*8(a-h,o-z)
#include "prgm.fh"
#include "warnings.fh"
#include "Molcas.fh"
parameter(ROUTINE='VV_First')
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
external IsFreeUnit
integer natom, i, j, irc, file, IsFreeUnit, ipCoord
real*8 DT_2, DTSQ2, Ekin, time, totimpl, RMS
character caption*15, lastline*80, filname*80
logical hybrid, qmmm
integer natom2
real*8, allocatable :: vel(:), xyz(:), force(:)
real*8, allocatable :: Mass(:), tstxyz(:)
character, allocatable :: atom(:)*2, atom2(:)*2
real*8, allocatable :: force2(:), xyz2(:)

if (IPRINT == INSANE) write(6,*) ' Entering ',ROUTINE

write(6,*) '*** First step of the Velocity Verlet algorithm ***'

! Check for QM/MM calculation

filname = 'comqum.dat'
call F_INQUIRE(filname,hybrid)
hybrid = .false.

! Read atom, their coordinates and forces

if (hybrid) then
  write(6,'(/,5X,A)') 'Perform QM/MM Molecular Dynamics'
  call DxRdNAtomHbrd(natom)
else
  call DxRdNAtomStnd(natom)
end if

call mma_allocate(vel,natom*3)
call mma_allocate(xyz,natom*3)
call mma_allocate(force,natom*3)
call mma_allocate(tstxyz,natom*3)
call mma_allocate(Mass,natom)
call mma_allocate(atom,natom)

if (hybrid) then
  call DxRdHbrd(natom,atom,xyz,force)
else
  call DxRdStnd(natom,atom,xyz,force)
end if
call Get_dScalar('MD_time',time)

! Read the velocities
call Get_Velocity(vel,3*natom)

! Initialize the Mass variable
call GetMassDx(Mass,natom)

! Check if reduced dimensionality
if (POUT /= 0) then
  call project_out_for(force,natom)
elseif (PIN /= natom*3) then
  call project_in_for(force,natom)
end if
!--------------------------------------------------------------------C
! CANONICAL ENSEMBLE
!--------------------------------------------------------------------C
if (THERMO == 2) then
  call NhcThermo(vel)
end if
!--------------------------------------------------------------------C
! Write out the old coordinates

call DxRdNAtomStnd(natom2)
call mma_allocate(xyz2,natom*3)
call mma_allocate(force2,natom*3)
call mma_allocate(atom2,natom)
call DxRdStnd(natom2,atom2,xyz2,force2)
!if (iPrint > VERBOSE) then
caption = 'Old Coordinates'
lastline = ''
call DxPtTableCo(caption,time,natom2,atom2,xyz2,lastline,Mass,force)
!end if

!Definition of the time step

DT_2 = DT/2.0d0
DTSQ2 = DT*DT_2

Ekin = 0.0d0
RMS = 0.0d0
totimpl = 0.0d0

do i=1,natom
  do j=1,3
    !** Root mean square deviation ****************************
    tstxyz(3*(i-1)+j) = xyz(3*(i-1)+j)
    !**********************************************************
    xyz(3*(i-1)+j) = xyz(3*(i-1)+j)+DT*vel(3*(i-1)+j)+DTSQ2*force(3*(i-1)+j)/Mass(i)
    !**********************************************************
    RMS = RMS+(tstxyz(3*(i-1)+j)-xyz(3*(i-1)+j))**2
    !**********************************************************
    Ekin = Ekin+5.0D-01*Mass(i)*(vel(3*(i-1)+j)**2)
    vel(3*(i-1)+j) = vel(3*(i-1)+j)+DT_2*force(3*(i-1)+j)/Mass(i)
    totimpl = totimpl+vel(3*(i-1)+j)*Mass(i)
  end do
end do

! Check if reduced dimensionality (should not be needed)
if (POUT /= 0) then
  call project_out_vel(vel,natom)
elseif (PIN /= natom*3) then
  call project_in_vel(vel,natom)
end if

call Add_Info('EKin',[EKin],1,6)

RMS = sqrt(RMS/natom)

! Update the Link-Atom position
! (Temporary solution uses ipCoord instead of xyz)

qmmm = .false.
call DecideOnESPF(qmmm)
if (qmmm) then
  call GetMem('Coordinates','ALLO','REAL',ipCoord,3*natom)
  call dcopy_(3*natom,xyz,1,Work(ipCoord),1)
  call LA_Morok(natom,Work(ipCoord),2)
  call dcopy_(3*natom,Work(ipCoord),1,xyz,1)
  call Free_Work(ipCoord)
end if

! Output

if (iPrint >= USUAL) then
  write(6,400) 'Molecular Dynamics specifications (time = ',time,' a.u.)'
  write(6,'(5X,A,/)') '=========================================='
  write(6,402) 'Kinetic energy',Ekin,'a.u.'
  write(6,405) 'Total linear momentum ',totimpl,'a.u.'
  write(6,402) 'RMS deviation ',RMS,'Bohr'
end if

time = time+DT

! Set the initialization flag and save the current time-value

call Put_dScalar('MD_Time',time)
#ifdef _HDF5_
call mh5_put_dset(dyn_time,time)
#endif

! Write out the new coordinates

caption = 'New Coordinates'
lastline = ''
call DxPtTableCo(caption,time,natom,atom,xyz,lastline,Mass,force)

! Write coordinates to output file

!#ifdef _DEBUGPRINT_
!write(6,*) ' Dynamix calls 2 DxCoord.'
!#endif
call DxCoord(natom,atom,xyz,hybrid)
!#ifdef _DEBUGPRINT_
!write(6,*) ' Dynamix back from 2 DxCoord.'
!#endif

! Save the new coordinates

if (hybrid) then
  call dscal_(natom*3,Angstrom,xyz,1)
  file = IsFreeUnit(81)
  filname = 'prmcrd2'
  call Molcas_Open(file,filname)
  write(file,403) natom
  write(file,404) (xyz(i),i=1,natom*3)
  close(file)
else
  call Put_Coord_Full(xyz,natom)
# ifdef _HDF5_
  call mh5_put_dset(dyn_geom,xyz)
# endif
end if

call Put_Velocity(vel,3*natom)
#ifdef _HDF5_
call mh5_put_dset(dyn_vel,vel)
#endif

call mma_deallocate(vel)
call mma_deallocate(xyz)
call mma_deallocate(force)
call mma_deallocate(tstxyz)
call mma_deallocate(Mass)
call mma_deallocate(atom)
call mma_deallocate(xyz2)
call mma_deallocate(force2)
call mma_deallocate(atom2)

! The return code is set in order to continue the loop

irc = _RC_ALL_IS_WELL_

return

400 format(5x,a,f8.1,a)
402 format(5x,a14,8x,d11.4,1x,a)
403 format(/,i5)
404 format(6f12.7)
405 format(5x,a22,d11.4,1x,a)

end subroutine VelVer_First
