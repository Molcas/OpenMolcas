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

subroutine decomp_aniso(A,Jiso,Jsym,Jantisym,dbg)

use Constants, only: Zero, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: A(3,3)
real(kind=wp), intent(out) :: Jiso, Jsym(3,3), Jantisym(3,3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j
real(kind=wp) :: Dtmp(3,3)

Jiso = Zero
Jsym(:,:) = Zero
Jantisym(:,:) = Zero
!-------------------------------------
do i=1,3
  Jiso = Jiso+A(i,i)
end do
Jiso = Jiso/Three
!-------------------------------------
do i=1,3
  Jsym(i,i) = A(i,i)-Jiso
end do

! find the symmetric matrix:
do i=1,3
  do j=1,3
    if (i /= j) Jsym(i,j) = (A(i,j)+A(j,i))*Half
  end do
end do
! find the anti-symmetric matrix:
do i=1,3
  do j=1,3
    if (i /= j) Jantisym(i,j) = (A(i,j)-A(j,i))*Half
  end do
end do

if (dbg) then
  Dtmp(:,:) = Zero
  do i=1,3
    Dtmp(i,i) = Jiso+Jsym(i,i)+Jantisym(i,i)
    do j=1,3
      if (i /= j) Dtmp(i,j) = Jsym(i,j)+Jantisym(i,j)
    end do
  end do

  write(u6,*)
  write(u6,*) 'J recovered = '
  do i=1,3
    write(u6,'(3F24.14)') (Dtmp(i,j),j=1,3)
  end do
end if

end subroutine decomp_aniso
