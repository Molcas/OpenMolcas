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

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: natom

call Get_nAtoms_Full(natom)

end subroutine DxRdNAtomStnd

!***********************************************************************
!
! Read the number of atoms. This needs to be done separately, because
! dynamic-size matrices depend on the result.

subroutine DxRdNAtomHbrd(natom)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: natom
integer(kind=iwp) :: filenum
integer(kind=iwp), external :: IsFreeUnit
character(len=80) :: filename

filenum = IsFreeUnit(81)
filename = 'fixforce.dmx'
call Molcas_Open(filenum,filename)
read(filenum,100) natom
close(filenum)

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

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom
character(len=2), intent(out) :: atom(natom)
real(kind=wp), intent(out) :: xyz(natom*3), force(natom*3)
real(kind=wp), parameter :: conv = -One

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

use Constants, only: Angstrom, auTokJ, rNAVO, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom
character(len=2), intent(out) :: atom(natom)
real(kind=wp), intent(out) :: xyz(natom*3), force(natom*3)
integer(kind=iwp) :: filenum, i, j, natom2
character(len=80) :: filename
real(kind=wp), parameter :: a2bohr = One/Angstrom, conv = Angstrom/(rNAVO*auTokJ)
integer(kind=iwp), external :: IsFreeUnit

! The parameter conv converts the forces from Hartree/Bohr
!                    to kJ/mole/Angstrom  => 1/4961.47525891
!             a2bohr converts the coordinates from Angstrom
!                    to Bohr              => 1/0.52917721090

filenum = IsFreeUnit(81)
filename = 'fixforce.dmx'
call Molcas_Open(filenum,filename)
read(filenum,*)
do i=1,natom
  read(filenum,101) (force((i-1)*3+j),j=1,3),atom(i)
end do
close(filenum)

! Convert forces from kJ/mole/Angstrom to Hartree/Bohr

call dscal_(3*natom,conv,force,1)

! Read coordinates from fiel 'prmcrd2' and check for consistency
! with forces.

filenum = IsFreeUnit(filenum)
filename = 'prmcrd2'
call Molcas_Open(filenum,filename)
read(filenum,'(/,I6)') natom2
if (natom2 /= natom) call WarningMessage(2,'Inconsistency between coordinates')
read(filenum,102) (xyz(i),i=1,natom*3)
close(filenum)

! Convert coordinates from Angstrom to Bohr

call dscal_(3*natom,a2bohr,xyz,1)

return

101 format(3es21.14,a)
102 format(6f12.7)

end subroutine DxRdHbrd

!***********************************************************************
!
! Prints a nice table in the output file

subroutine DxPtTableCo(caption,time,natom,atom,xyz,lastline,M,fo)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: natom
real(kind=wp), intent(in) :: time, xyz(natom*3), M(natom), fo(natom*3)
character(len=15), intent(in) :: caption
character(len=2), intent(in) :: atom(natom)
character(len=80), intent(in) :: lastline
integer(kind=iwp) :: i, j, k

write(u6,*)
write(u6,*)
write(u6,*)

write(u6,100) caption,' (time = ',time,' a.u.):'
write(u6,102) '----------------------------------------------------------------------------------------------'
write(u6,102) '      No. Atom    X          Y          Z        Mass       F(x)         F(y)         F(z)'
write(u6,102) '----------------------------------------------------------------------------------------------'

do i=1,natom
  write(u6,101) i,atom(i),(xyz(3*(i-1)+j),j=1,3),M(i),(fo(3*(i-1)+k),k=1,3)
end do

write(u6,102) '----------------------------------------------------------------------------------------------'
write(u6,102) trim(lastline)

write(u6,*)
write(u6,*)
write(u6,*)

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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: natom
real(kind=wp), intent(in) :: time, xyz(natom*3)
character(len=15), intent(in) :: caption
character(len=2), intent(in) :: atom(natom)
integer(kind=iwp) :: i, j

write(u6,*)
write(u6,*)
write(u6,*)

write(u6,102) caption,' (time = ',time,' a.u.):'
write(u6,*) '----------------------------------------------'
write(u6,*) '     No.  Atom    X          Y          Z     '
write(u6,*) '----------------------------------------------'

do i=1,natom
  write(u6,103) '      ',i,atom(i),(xyz(3*(i-1)+j),j=1,3)
end do

write(u6,*) '----------------------------------------------'

write(u6,*)
write(u6,*)
write(u6,*)

return

102 format(a22,a7,f8.1,a)
103 format(a6,i4,a3,3f11.6)

end subroutine DxPtTableWithoutMassForce

!***********************************************************************
!
! This Subroutine reads in the coordinates to project out from the file 'out.00N.xyz'
!
subroutine DxRdOut(pcoo,POUT,natom)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: POUT, natom
real(kind=wp), intent(out) :: pcoo(POUT,natom*3)
integer(kind=iwp) :: i, j, p, filenum, mn
character(len=80) :: filename
character(len=180) :: OutLine
real(kind=wp), allocatable :: UMat(:,:), VMat(:,:), S(:), SqrtM(:)
integer(kind=iwp), external :: IsFreeUnit
character(len=180), external :: Get_Ln

filenum = 81
do p=1,POUT
  pcoo(p,:) = Zero
  write(filename,'(A,I3.3,A)') 'out.',p,'.xyz'
  filenum = IsFreeUnit(filenum)
  call Molcas_Open(filenum,filename)
  do i=1,natom
    OutLine = Get_Ln(filenum)
    call UpCase(OutLine)
    do j=1,3
      call Get_F(j,pcoo(p,3*(i-1)+j),1)
    end do
  end do
  close(filenum)
end do

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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: PIN, natom
real(kind=wp), intent(out) :: pcoo(PIN,natom*3)
integer(kind=iwp) :: i, j, p, filenum, mn
character(len=80) :: filename
character(len=180) :: InLine
real(kind=wp), allocatable :: UMat(:,:), VMat(:,:), S(:), SqrtM(:)
integer(kind=iwp), external :: IsFreeUnit
character(len=180), external :: Get_Ln

filenum = 81
do p=1,PIN
  pcoo(p,:) = Zero
  write(filename,'(A,I3.3,A)') 'in.',p,'.xyz'
  filenum = IsFreeUnit(filenum)
  call Molcas_Open(filenum,filename)
  do i=1,natom
    InLine = Get_Ln(filenum)
    call UpCase(InLine)
    do j=1,3
      call Get_F(j,pcoo(p,3*(i-1)+j),1)
    end do
  end do
  close(filenum)
end do

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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom
real(kind=wp), intent(out) :: vel(natom*3)
integer(kind=iwp) :: i, j, filenum
character(len=80) :: filename
character(len=180) :: VelLine
integer(kind=iwp), external :: IsFreeUnit
character(len=180), external :: Get_Ln

filenum = IsFreeUnit(81)
filename = 'velocity.xyz'
call Molcas_Open(filenum,filename)
vel(:) = Zero
do i=1,natom
  VelLine = Get_Ln(filenum)
  call UpCase(VelLine)
  do j=1,3
    call Get_F(j,vel(3*(i-1)+j),1)
  end do
end do
!read(filenum,100) (vel(i),i=1,3*natom3)
close(filenum)

return

!100 format(3es18.10)

end subroutine DxRdVel

!***********************************************************************
!
! This Subroutine writes the velocities to the file 'velocity.xyz'
!
subroutine DxWtVel(vel,natom3)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom3
real(kind=wp), intent(in) :: vel(natom3)
integer(kind=iwp) :: i, filenum
character(len=80) :: filename
integer(kind=iwp), external :: IsFreeUnit

filenum = IsFreeUnit(81)
filename = 'velocity.xyz'
call Molcas_Open(filenum,filename)
write(filenum,100) (vel(i),i=1,natom3)
close(filenum)

return

100 format(3es18.10)

end subroutine DxWtVel

!***********************************************************************
!
! Write the coordinates to a file in xyz format
!
subroutine DxCoord(natom,atom,xyz,hybrid)

use Constants, only: Angstrom
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom
character(len=2), intent(in) :: atom(natom)
real(kind=wp), intent(in) :: xyz(natom*3)
logical(kind=iwp), intent(in) :: hybrid
integer(kind=iwp) :: i, j, filenum
character(len=9) :: filename
logical(kind=iwp) :: exists
integer(kind=iwp), external :: IsFreeUnit

if (.not. hybrid) then
  filenum = IsFreeUnit(81)
  filename = 'md.xyz'
  call OpnFl(filename,filenum,exists)
  call Append_file(filenum)
  write(filenum,'(I5,/)') natom
  do i=1,natom
    write(filenum,100) atom(i),(Angstrom*xyz(3*(i-1)+j),j=1,3)
  end do
  close(filenum)
else
  filenum = IsFreeUnit(81)
  filename = 'md.prmcrd'
  call OpnFl(filename,filenum,exists)
  call Append_file(filenum)
  write(filenum,'(/,I6)') natom
  write(filenum,'(6F12.7)') (Angstrom*xyz(i),i=1,3*natom)
  close(filenum)
  filenum = IsFreeUnit(filenum)
  filename = 'vmd.mdcrd'
  call OpnFl(filename,filenum,exists)
  call Append_file(filenum)
  write(filenum,'(10F8.3)') (Angstrom*xyz(i),i=1,3*natom)
  close(filenum)
end if

return

100 format(1x,a2,3f15.8)

end subroutine DxCoord

!***********************************************************************
!
! Write out a summary of energies in CSV format
!
subroutine DxEnergies(time,Epot,Ekin,Etot)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: time, Epot, Ekin, Etot
integer(kind=iwp) :: filenum, i, n, nEnergies
logical(kind=iwp) :: exists, RootCheck
character(len=12) :: filename
character(len=24) :: frmt
real(kind=wp), allocatable :: Energies(:)
integer(kind=iwp), external :: IsFreeUnit

filename = 'md.energies'

RootCheck = .false.

call Qpg_iScalar('Relax CASSCF root',RootCheck)

filenum = IsFreeUnit(81)
call OpnFl(filename,filenum,exists)
call Append_file(filenum)
if (.not. exists) then
  write(filenum,'(3X,A4,10X,A4,2(17X,A4))') 'time','Epot','Ekin','Etot'
end if
frmt = '(f8.2,   (2x,es19.12))'

if (RootCheck) then
  call Get_iScalar('Number of roots',nEnergies)
  call mma_allocate(Energies,nEnergies,label='MS energies')
  call Get_dArray('Last energies',Energies,nEnergies)
  n = nEnergies+3
  write(frmt(7:9),'(i3)') n
  write(filenum,frmt) time,Epot,Ekin,Etot,(Energies(i),i=1,nEnergies)
  call mma_deallocate(Energies)
else
  n = 3
  write(frmt(7:9),'(i3)') n
  write(filenum,frmt) time,Epot,Ekin,Etot
end if
close(filenum)

return

end subroutine DxEnergies

!***********************************************************************
!
! Writes the velocities on RUNFILE
!
subroutine Put_Velocity(vel,natom3)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom3
real(kind=wp), intent(in) :: vel(natom3)

call Put_dArray('Velocities',vel,natom3)

return

end subroutine Put_Velocity

!***********************************************************************
!
! Reads the velocities from RUNFILE
!
subroutine Get_Velocity(vel,natom3)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom3
real(kind=wp), intent(out) :: vel(natom3)

call Get_dArray('Velocities',vel,natom3)

return

end subroutine Get_Velocity

!***********************************************************************
!
!     Reads the extra degrees of freedom from RUNFILE
!
subroutine Get_NHC(NHC,nh)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh
real(kind=wp), intent(out) :: NHC(nh)

call Get_dArray('NOSEHOOVER',NHC,nh)

return

end subroutine Get_NHC

!***********************************************************************
!
! Writes the extra degrees of Freedom on RUNFILE
!
subroutine Put_NHC(NHC,nh)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh
real(kind=wp), intent(in) :: NHC(nh)

call Put_dArray('NOSEHOOVER',NHC,nh)

return

end subroutine Put_NHC
