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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Sort_Localisation_1(CMO,U,nBas,nOcc)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: sort CMO columns according to U.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas, nOcc
real(kind=wp), intent(inout) :: CMO(nBas,nOcc)
real(kind=wp), intent(in) :: U(nOcc,nOcc)
integer(kind=iwp) :: i, j, jmax
real(kind=wp) :: Umax, Utst
integer(kind=iwp), allocatable :: I1(:), I2(:)
real(kind=wp), allocatable :: C(:,:)

! Allocations.
! ------------

call mma_allocate(I1,nOcc,label='Sr1I1')
call mma_allocate(I2,nOcc,label='Sr1I2')
call mma_allocate(C,nBas,nOcc,label='Sr1C')

! Find max U element in each row.
! -------------------------------

do i=1,nOcc
  I1(i) = i
end do

do i=1,nOcc
  jmax = 0
  Umax = -huge(Umax)
  do j=1,nOcc
    if (I1(j) == j) then
      Utst = abs(U(i,j))
      if (Utst > Umax) then
        jmax = j
        Umax = Utst
      end if
    end if
  end do
  if (jmax == 0) then
    call SysAbendMsg('Sort_Localisation_1','Error:','jmax=0')
  else
    I1(jmax) = 0
    I2(i) = jmax
  end if
end do

! Swap MOs according to I2.
! -------------------------

C(:,:) = CMO
do i=1,nOcc
  CMO(:,i) = C(:,I2(i))
end do

! De-allocate.
! ------------

call mma_deallocate(I1)
call mma_deallocate(I2)
call mma_deallocate(C)

end subroutine Sort_Localisation_1
