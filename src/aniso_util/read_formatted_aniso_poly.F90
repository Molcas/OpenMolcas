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

use Constants, only: Zero, cZero
use Definitions, only: wp

implicit none
#include "stdalloc.fh"
integer, intent(inout) :: nss, nstate, nLoc, iReturn
! nLoc is the maximal value of the array nss(1:nneq)
real(kind=8), intent(out) :: eso(nLoc)
complex(kind=8), intent(out) :: MM(3,nLoc,nLoc)
complex(kind=8), intent(out) :: MS(3,nLoc,nLoc)
character(len=180) :: input_file_name
! local variables:
integer :: l, j, j1, j2, LuAniso, IsFreeUnit
integer :: multiplicity(nLoc)
real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:)
external :: IsFreeUnit

! set to zero all arrays:
iReturn = 0
multiplicity = 0
eso(:) = Zero
MM(:,:,:) = cZero
MS(:,:,:) = cZero
! read the file "aniso.input":
LuAniso = IsFreeUnit(40)
call molcas_open(LuAniso,input_file_name)
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate,nss
read(LuAniso,*) (eso(j),j=1,nss)
read(LuAniso,*) (multiplicity(j),j=1,nstate)

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')
! magnetic moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  do j1=1,nss
    do j2=1,nss
      MM(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),kind=wp)
    end do
  end do
end do

! spin moment
do l=1,3
  tmpR(:,:) = Zero
  tmpI(:,:) = Zero
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  do j1=1,nss
    do j2=1,nss
      MS(l,j1,j2) = cmplx(tmpR(j1,j2),tmpI(j1,j2),kind=wp)
    end do
  end do
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

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_integer_array(multiplicity)
#endif

end subroutine read_formatted_aniso_poly
