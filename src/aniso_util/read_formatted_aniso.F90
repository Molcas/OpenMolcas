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

subroutine read_formatted_aniso(input_file_name,nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(inout) :: nss, nstate
integer, intent(out) :: multiplicity(nstate)
real(kind=8), intent(out) :: eso(nss), esfs(nstate)
real(kind=8), intent(out) :: edmom(3,nstate,nstate)
real(kind=8), intent(out) :: angmom(3,nstate,nstate)
real(kind=8), intent(out) :: amfi(3,nstate,nstate)
complex(kind=8), intent(out) :: MM(3,nss,nss)
complex(kind=8), intent(out) :: MS(3,nss,nss)
complex(kind=8), intent(out) :: ML(3,nss,nss)
! electric dipole moment
complex(kind=8), intent(out) :: DM(3,nss,nss)
complex(kind=8), intent(out) :: U(nss,nss)
complex(kind=8), intent(out) :: HSO(nss,nss)
character(Len=180) :: input_file_name
! local variables:
integer :: l, j, j1, j2, LuAniso, IsFreeUnit
real(kind=8) :: g_e
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)
external :: IsFreeUnit
logical :: DBG

dbg = .false.

if (dbg) write(6,'(A)') 'Entering read_formatted_aniso'
call xFlush(6)
g_e = 2.00231930437180_wp
!   set to zero all arrays:
multiplicity = 0
eso = 0.0_wp
esfs = 0.0_wp
edmom = 0.0_wp
angmom = 0.0_wp
MM = (0.0_wp,0.0_wp)
MS = (0.0_wp,0.0_wp)
ML = (0.0_wp,0.0_wp)
DM = (0.0_wp,0.0_wp)
U = (0.0_wp,0.0_wp)
! read the file "aniso.input":
LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,trim(input_file_name))
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate,nss
if (dbg) write(6,'(A,2I6)') 'nstate, nss:',nstate,nss
call xFlush(6)
read(LuAniso,*) (eso(j),j=1,nss)
if (dbg) then
  write(6,'(A)') 'ESO:'
  write(6,'(5ES24.14)') (eso(j),j=1,nss)
end if
read(LuAniso,*) (multiplicity(j),j=1,nstate)
if (dbg) then
  write(6,'(A)') '(multiplicity(j),j=1,nstate)'
  write(6,'(50I3)') (multiplicity(j),j=1,nstate)
end if
call xFlush(6)

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')
! magnetic moment
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  do j1=1,nss
    do j2=1,nss
      MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
    end do
  end do
end do
call xFlush(6)

! spin moment
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  do j1=1,nss
    do j2=1,nss
      MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
    end do
  end do
end do

! spin-free energies
read(LuAniso,*) (esfs(j),j=1,nstate)
if (dbg) then
  write(6,'(A)') 'ESFS:'
  write(6,'(5ES24.14)') (esfs(j),j=1,nstate)
end if

! U matrix
tmpR = 0.0_wp
tmpI = 0.0_wp
do j1=1,nss
  read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
end do

do j1=1,nss
  do j2=1,nss
    U(j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
  end do
end do

! angmom
do l=1,3
  do j1=1,nstate
    read(LuAniso,*) (angmom(l,j1,j2),j2=1,nstate)
  end do
end do

! compute the orbital moment
do l=1,3
  do j1=1,nss
    do j2=1,nss
      ML(l,j1,j2) = -MM(l,j1,j2)-MS(l,j1,j2)*g_e
    end do
  end do
end do

! DMmom
do l=1,3
  tmpR = 0.0_wp
  tmpI = 0.0_wp
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  do j1=1,nss
    do j2=1,nss
      DM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
    end do
  end do
end do

! edmom
do l=1,3
  do j1=1,nstate
    read(LuAniso,*) (edmom(l,j1,j2),j2=1,nstate)
  end do
end do

! amfi
do l=1,3
  do j1=1,nstate
    read(LuAniso,*) (amfi(l,j1,j2),j2=1,nstate)
  end do
end do

! HSO matrix
tmpR = 0.0_wp
tmpI = 0.0_wp
do j1=1,nss
  read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
end do

do j1=1,nss
  do j2=1,nss
    HSO(j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),wp)
  end do
end do

call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

close(LuAniso)

return

end subroutine read_formatted_aniso
