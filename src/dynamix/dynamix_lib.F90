!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! Library of routines for Dynamix module
!
!***********************************************************************
!      SUBROUTINE DxRdNAtomStnd(natom)
!      SUBROUTINE DxRdNAtomHbrd(natom)
!      SUBROUTINE DxRdStnd(natom,atom,xyz,force)
!      SUBROUTINE DxRdHbrd(natom,atom,xyz,force)
!      SUBROUTINE DxPtTableCo(caption,time,natom,atom,xyz,lastline,M,fo)
!      SUBROUTINE DxPtTableWithoutMassForce(caption,time,natom,atom,xyz)
!      SUBROUTINE DxRdOut(pcoo,POUT,natom)
!      SUBROUTINE DxRdVel(vel,natom)
!      SUBROUTINE DxWtVel(vel,natom3)
!      SUBROUTINE DxCoord(natom,atom,xyz,hybrid)
!      SUBROUTINE DxEnergies(time,Epot,Ekin,Etot)
!      SUBROUTINE Put_Velocity(vel,natom3)
!      SUBROUTINE Get_Velocity(vel,natom3)
!      SUBROUTINE Get_NHC(NHC,nh)
!      SUBROUTINE Put_NHC(NHC,nh)
!***********************************************************************
!
! Read the number of atoms. This needs to be done separately, because
! dynamic-size matrices depend on the result.

subroutine DxRdNAtomStnd(natom)

implicit none
integer natom

call Get_nAtoms_Full(natom)

end subroutine DxRdNAtomStnd

!***********************************************************************
!
! Read the number of atoms. This needs to be done separately, because
! dynamic-size matrices depend on the result.

subroutine DxRdNAtomHbrd(natom)

implicit none
external IsFreeUnit
integer natom, file, IsFreeUnit
character filname*80

file = 81
file = IsFreeUnit(file)
filname = 'fixforce.dmx'
call Molcas_Open(file,filname)
read(file,100) natom
close(file)

return

100 format(i6)

end subroutine DxRdNAtomHbrd

!***********************************************************************
!
! Read in the atom names, coordinates and forces from RUNFILE for
! the standard QM molecular dynamics.
!
! IS 14/06-2007

subroutine DxRdStnd(natom,atom,xyz,force)

implicit none
#include "Molcas.fh"
#include "WrkSpc.fh"
integer natom
real*8 xyz(natom*3), force(natom*3), conv
character atom(natom)*2
parameter(conv=-1.0d0)

! The parameter conv converts the gradients (Hartree/Bohr)
!                              to forces (Hartree/Bohr)

call Get_Name_Full(atom)
call Get_Coord_Full(xyz,natom)

! Read the gradients from RUNFILE and convert them to forces.

call Get_Grad_Full(force,natom)
call dscal_(3*natom,conv,force,1)

end subroutine DxRdStnd

!***********************************************************************
!
! Read in the atom names, coordinates and forces from files
! fixforce.dmx and prmcrd2.

subroutine DxRdHbrd(natom,atom,xyz,force)

implicit none
#include "Molcas.fh"
#include "constants2.fh"
external IsFreeUnit
integer i, j, natom, natom2, file, IsFreeUnit
real*8 xyz(natom*3), force(natom*3), conv, a2bohr
character filname*80, atom(natom)*2
parameter(a2bohr=1.0d0/Angstrom)
parameter(conv=Angstrom/CONV_AU_TO_KJ_PER_MOLE_)

! The parameter conv converts the forces from Hartree/Bohr
!                    to kJ/mole/Agstrom   => 1/4961.475514610d0
!             a2bohr converts the coordinates from Angstrom
!                    to Bohr              => 1/0.52917720859

file = 81
file = IsFreeUnit(file)
filname = 'fixforce.dmx'
call Molcas_Open(file,filname)
read(file,*)
do i=1,natom
  read(file,101) (force((i-1)*3+j),j=1,3),atom(i)
end do
close(file)

! Convert forces from kJ/mole/Angstrom to Hartree/Bohr

call dscal_(3*natom,conv,force,1)

! Read coordinates from fiel 'prmcrd2' and check for consistency
! with forces.

file = IsFreeUnit(file)
filname = 'prmcrd2'
call Molcas_Open(file,filname)
read(file,'(/,I6)') natom2
if (natom2 /= natom) stop 'Inconsistency between coordinates'
read(file,102) (xyz(i),i=1,natom*3)
close(file)

! Convert coordinates from Angstrom to Bohr

call dscal_(3*natom,a2bohr,xyz,1)

return

101 format(3e21.14,a)
102 format(6f12.7)

end subroutine DxRdHbrd

!***********************************************************************
!
! Prints a nice table in the output file

subroutine DxPtTableCo(caption,time,natom,atom,xyz,lastline,M,fo)

implicit none
#include "Molcas.fh"
integer i, j, natom, k
character caption*15, lastline*80, atom(natom)*2
real*8 time, xyz(natom*3), M(natom), fo(natom*3)

do i=1,3
  write(6,*)
end do

write(6,100) caption,' (time = ',time,' a.u.):'
write(6,102) '----------------------------------------------------------------------------------------------'
write(6,102) '      No. Atom    X          Y          Z        Mass       F(x)         F(y)         F(z)'
write(6,102) '----------------------------------------------------------------------------------------------'

do i=1,natom
  write(6,101) i,atom(i),(xyz(3*(i-1)+j),j=1,3),M(i),(fo(3*(i-1)+k),k=1,3)
end do

write(6,102) '----------------------------------------------------------------------------------------------'
write(6,102) trim(lastline)

do i=1,3
  write(6,*)
end do

return

100 format(a22,a7,f8.1,a)
101 format(6x,i4,a3,3(1x,f10.6),1x,es9.2,3(1x,es12.5))
102 format(1x,a)

end subroutine DxPtTableCo

!***********************************************************************
!
! Prints a nice table in the output file for the initial velcoities
!
subroutine DxPtTableWithoutMassForce(caption,time,natom,atom,xyz)

implicit none
#include "Molcas.fh"
integer i, j, natom
character caption*15, atom(natom)*2
real*8 time, xyz(natom*3)

do i=1,3
  write(6,*) ' '
end do

write(6,102) caption,' (time = ',time,' a.u.):'
write(6,*) '----------------------------------------------'
write(6,*) '     No.  Atom    X          Y          Z     '
write(6,*) '----------------------------------------------'

do i=1,natom
  write(6,103) '      ',i,atom(i),(xyz(3*(i-1)+j),j=1,3)
end do

write(6,*) '----------------------------------------------'

do i=1,3
  write(6,*) ' '
end do

return

102 format(a22,a7,f8.1,a)
103 format(a6,i4,a3,3f11.6)

end subroutine DxPtTableWithoutMassForce

!***********************************************************************
!
! This Subroutine reads in the coordinates to project out from the file 'out.00N.xyz'
!
subroutine DxRdOut(pcoo,POUT,natom)

implicit none
#include "stdalloc.fh"
#include "real.fh"
external IsFreeUnit
integer i, j, p, file, natom, POUT, IsFreeUnit, mn
real*8 pcoo(POUT,natom*3)
character*80 filname
character*180 OutLine, Get_Ln
real*8, allocatable :: UMat(:,:), VMat(:,:), S(:), SqrtM(:)

file = 81
file = IsFreeUnit(file)
do P=1,POUT
  write(filname,'(A,I3.3,A)') 'out.',p,'.xyz'
  call Molcas_Open(file,filname)
  do i=1,natom
    OutLine = Get_Ln(file)
    call UpCase(OutLine)
    do j=1,3
      pcoo(p,3*(i-1)+j) = 0.0d0
      call Get_F(j,pcoo(p,3*(i-1)+j),1)
    end do
  end do
end do
close(file)

! Orthonormalize the vectors in mass-weighted coordinates

call mma_Allocate(SqrtM,natom)
call GetMassDx(SqrtM,natom)
do i=1,natom
  SqrtM(i) = sqrt(SqrtM(i))
  j = 3*(i-1)+1
  pcoo(:,j:j+2) = pcoo(:,j:j+2)*SqrtM(i)
end do
mn = min(POUT,3*natom)
call mma_Allocate(UMat,POUT,mn)
call mma_Allocate(VMat,mn,3*natom)
call mma_Allocate(S,mn)
call large_svd(POUT,3*natom,pcoo,UMat,VMat,S)
call dgemm_('N','N',POUT,3*natom,mn,One,UMat,POUT,VMat,mn,Zero,pcoo,POUT)
do i=1,natom
  j = 3*(i-1)+1
  pcoo(:,j:j+2) = pcoo(:,j:j+2)/SqrtM(i)
end do
call mma_deAllocate(UMat)
call mma_deAllocate(VMat)
call mma_deAllocate(S)
call mma_deAllocate(SqrtM)

return

end subroutine DxRdOut

!***********************************************************************
!
!     This Subroutine reads in the coordinates to keep in from the file 'in.00N.xyz'
!
subroutine DxRdIn(pcoo,PIN,natom)

implicit none
#include "stdalloc.fh"
#include "real.fh"
external IsFreeUnit
integer i, j, p, file, natom, PIN, IsFreeUnit, mn
real*8 pcoo(PIN,natom*3)
character*80 filname
character*180 OutLine, Get_Ln
real*8, allocatable :: UMat(:,:), VMat(:,:), S(:), SqrtM(:)

file = 81
file = IsFreeUnit(file)
do P=1,PIN
  write(filname,'(A,I3.3,A)') 'in.',p,'.xyz'
  call Molcas_Open(file,filname)
  do i=1,natom
    OutLine = Get_Ln(file)
    call UpCase(OutLine)
    do j=1,3
      pcoo(p,3*(i-1)+j) = 0.0d0
      call Get_F(j,pcoo(p,3*(i-1)+j),1)
    end do
  end do
end do
close(file)

! Orthonormalize the vectors in mass-weighted coordinates

call mma_Allocate(SqrtM,natom)
call GetMassDx(SqrtM,natom)
do i=1,natom
  SqrtM(i) = sqrt(SqrtM(i))
  j = 3*(i-1)+1
  pcoo(:,j:j+2) = pcoo(:,j:j+2)*SqrtM(i)
end do
mn = min(PIN,3*natom)
call mma_Allocate(UMat,PIN,mn)
call mma_Allocate(VMat,mn,3*natom)
call mma_Allocate(S,mn)
call large_svd(PIN,3*natom,pcoo,UMat,VMat,S)
call dgemm_('N','N',PIN,3*natom,mn,One,UMat,PIN,VMat,mn,Zero,pcoo,PIN)
do i=1,natom
  j = 3*(i-1)+1
  pcoo(:,j:j+2) = pcoo(:,j:j+2)/SqrtM(i)
end do
call mma_deAllocate(UMat)
call mma_deAllocate(VMat)
call mma_deAllocate(S)
call mma_deAllocate(SqrtM)

return

end subroutine DxRdIn

!***********************************************************************
!
!     This Subroutine reads in the velocities from the file 'velocity.xyz'
!
subroutine DxRdVel(vel,natom)

implicit none
#include "Molcas.fh"
external IsFreeUnit
integer i, j, file, natom, IsFreeUnit
real*8 vel(natom*3)
character*80 filname
character*180 VelLine, Get_Ln

file = 81
file = IsFreeUnit(file)
filname = 'velocity.xyz'
call Molcas_Open(file,filname)
do i=1,natom
  VelLine = Get_Ln(file)
  call UpCase(VelLine)
  do j=1,3
    vel(3*(i-1)+j) = 0.0d0
    call Get_F(j,vel(3*(i-1)+j),1)
  end do
end do
!read(file,100) (vel(i),i=1,3*natom3)
close(file)

return

!100 format(3e18.10)

end subroutine DxRdVel

!***********************************************************************
!
! This Subroutine writes the velocities to the file 'velocity.xyz'
!
subroutine DxWtVel(vel,natom3)

implicit none
#include "Molcas.fh"
external IsFreeUnit
integer i, file, natom3, IsFreeUnit
real*8 vel(natom3)
character*80 filname

file = 81
file = IsFreeUnit(file)
filname = 'velocity.xyz'
call Molcas_Open(file,filname)
write(file,100) (vel(i),i=1,natom3)
close(file)

return

100 format(3e18.10)

end subroutine DxWtVel

!***********************************************************************
!
! Write the coordinates to a file in xyz format
!
subroutine DxCoord(natom,atom,xyz,hybrid)

implicit none
#include "Molcas.fh"
#include "MD.fh"
#include "constants2.fh"
external IsFreeUnit
integer i, j, file, natom, IsFreeUnit
real*8 xyz(natom*3)
character atom(natom)*2, filename*9
logical hybrid, Exist

if (.not. hybrid) then
  file = 82
  file = IsFreeUnit(file)
  filename = 'md.xyz'
  call OpnFl(filename,file,Exist)
  call Append_file(file)
  write(file,'(I5,/)') natom
  do i=1,natom
    write(file,100) atom(i),(Angstrom*xyz(3*(i-1)+j),j=1,3)
  end do
  close(file)
else
  file = IsFreeUnit(file)
  filename = 'md.prmcrd'
  call OpnFl(filename,file,Exist)
  call Append_file(file)
  write(file,'(/,I6)') natom
  write(file,'(6F12.7)') (Angstrom*xyz(i),i=1,3*natom)
  close(file)
  file = IsFreeUnit(file)
  filename = 'vmd.mdcrd'
  call OpnFl(filename,file,Exist)
  call Append_file(file)
  write(file,'(10F8.3)') (Angstrom*xyz(i),i=1,3*natom)
  close(file)
end if

return

100 format(1x,a2,3f15.8)

end subroutine DxCoord

!***********************************************************************
!
! Write out a summary of energies in CSV format
!
subroutine DxEnergies(time,Epot,Ekin,Etot)
implicit none
external IsFreeUnit
#include "WrkSpc.fh"
integer file, nEnergies, i, n, ipEnergies
integer IsFreeUnit
real*8 time, Epot, Ekin, Etot
character filename*12, frmt*24
logical Exist, RootCheck

filename = 'md.energies'

RootCheck = .false.

call Qpg_iScalar('Relax CASSCF root',RootCheck)

file = 82
file = IsFreeUnit(file)
call OpnFl(filename,file,Exist)
call Append_file(file)
if (.not. Exist) then
  write(file,'(3X,A4,10X,A4,2(17X,A4))') 'time','Epot','Ekin','Etot'
end if
frmt = '(F8.2, (2X,D19.12))'

if (RootCheck) then
  call Get_iScalar('Number of roots',nEnergies)
  call GetMem('MS energies','ALLO','REAL',ipEnergies,nEnergies)
  call Get_dArray('Last energies',Work(ipEnergies),nEnergies)
  n = nEnergies+3
  write(frmt(7:7),'(I1)') n
  write(file,frmt) time,Epot,Ekin,Etot,(Work(ipEnergies-1+i),i=1,nEnergies)
  call GetMem('MS energies','FREE','REAL',ipEnergies,nEnergies)
else
  n = 3
  write(frmt(7:7),'(I1)') n
  write(file,frmt) time,Epot,Ekin,Etot
end if
close(file)

return

end subroutine DxEnergies

!***********************************************************************
!
!     Writes the velocities on RUNFILE
!
subroutine Put_Velocity(vel,natom3)

implicit none
#include "WrkSpc.fh"
#include "Molcas.fh"
integer natom3
real*8 vel(natom3)

call Put_dArray('Velocities',vel,natom3)
!call GetMem('Velocities','FREE','REAL',vel,natom3)

return

end subroutine Put_Velocity

!***********************************************************************
!
! Reads the velocities from RUNFILE
!
subroutine Get_Velocity(vel,natom3)

implicit none
#include "WrkSpc.fh"
#include "Molcas.fh"
integer natom3
real*8 vel(natom3)

!call GetMem('Velocities','ALLO','REAL',vel,natom3)
call Get_dArray('Velocities',vel,natom3)

return

end subroutine Get_Velocity

!***********************************************************************
!
!     Reads the extra degrees of freedom from RUNFILE
!
subroutine Get_NHC(NHC,nh)

implicit none
#include "WrkSpc.fh"
#include "Molcas.fh"
integer nh
real*8 NHC(nh)

!call GetMem('NOSEHOOVER','ALLO','REAL',NHC,nh)
call Get_dArray('NOSEHOOVER',NHC,nh)

return

end subroutine Get_NHC

!***********************************************************************
!
! Writes the extra degrees of Freedom on RUNFILE
!
subroutine Put_NHC(NHC,nh)

implicit none
#include "WrkSpc.fh"
#include "Molcas.fh"
integer nh
real*8 NHC(nh)

call Put_dArray('NOSEHOOVER',NHC,nh)
!call GetMem('NOSEHOOVER','FREE','REAL',NHC,nh)

return

end subroutine Put_NHC
