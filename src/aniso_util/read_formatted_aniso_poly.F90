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

subroutine read_formatted_aniso_poly(input_file_name,nss,nstate,nLoc,eso,MM,MS,iReturn)
! nLoc is the maximal value of the array nss(1:nneq)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero
use Definitions, only: wp, iwp

implicit none
character(len=180), intent(in) :: input_file_name
integer(kind=iwp), intent(inout) :: nss, nstate
integer(kind=iwp), intent(in) :: nLoc
real(kind=wp), intent(out) :: eso(nLoc)
complex(kind=wp), intent(out) :: MM(3,nLoc,nLoc), MS(3,nLoc,nLoc)
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: j, j1, j2, l, LuAniso
integer(kind=iwp), allocatable :: multiplicity(:)
real(kind=wp), allocatable :: tmpR(:,:), tmpI(:,:)
integer(kind=iwp), external :: IsFreeUnit

#include "macros.fh"

! set to zero all arrays:
iReturn = 0
MM(:,:,:) = cZero
MS(:,:,:) = cZero
! read the file "aniso.input":
LuAniso = IsFreeUnit(40)
call molcas_open(LuAniso,input_file_name)
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate,nss
read(LuAniso,*) (eso(j),j=1,nss)
call mma_allocate(multiplicity,nstate,label='multiplicity')
read(LuAniso,*) (multiplicity(j),j=1,nstate)
unused_var(multiplicity)

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')
! magnetic moment
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  MM(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

! spin moment
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  MS(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do
call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

!if (iprint > 4) then
!  write(u6,'(10a12)') (('------------'),j=1,10)
!  write(u6,'(15x,a,i2,a)') 'aniso_.input'
!  write(u6,'(10a12)') (('------------'),j=1,10)
!  write(u6,'(5x,a,i6)') 'nstate = ',nstate
!  write(u6,'(5x,a,i6)') '   nss = ',nss
!  write(u6,'(5x,a)') ' eso(j): '
!  write(u6,'(10(f12.5,1x))') (eso(j),j=1,nss)
!  write(u6,'(5x,a,i2,a)') 'multiplicity(j):'
!  write(u6,'(40i3)') (multiplicity(j),j=1,nstate)
!  write(u6,'(5x,a,i2,a)') 'dipso(l,j1,j2):'
!  do l=1,3
!    write(u6,'(5x,a,i2)') 'axis= ',l
!    do j1=1,nss
!      write(u6,'(20(2f20.14,1x))') (MM(l,j1,j2),j2=1,nss)
!    end do
!  end do
!  write(u6,*)
!  write(u6,'(5x,a,i2,a)') 's_so(l,j1,j2):'
!  do l=1,3
!    write(u6,'(5x,a,i2)') 'axis= ',l
!    do j1=1,nss
!      write(u6,'(20(2f20.14,1x))') (MS(l,j1,j2),j2=1,nss)
!    end do
!  end do
!end if
close(LuAniso)

call mma_deallocate(multiplicity)

return

end subroutine read_formatted_aniso_poly
