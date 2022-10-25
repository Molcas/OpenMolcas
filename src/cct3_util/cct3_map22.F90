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

subroutine cct3_map22(a,b,dimp,dimq,dim_1,dim_2,p,nfact)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dim_1, dim_2, p, nfact
real(kind=wp), intent(in) :: a(dimp,dimq)
real(kind=wp), intent(out) :: b(dim_1,dim_2)
integer(kind=iwp) :: pp !, indx(2)

if (nfact == 1) then

  ! nfact = + 1

  !do qq=1,dimq
  !  indx(q) = qq
  !  do pp=1,dimp
  !    indx(p) = pp
  !    b(indx(1),indx(2)) = a(pp,qq)
  !  end do
  !end do

  if (p == 1) then
    !1* (2)
    b(1:dimp,1:dimq) = a
  else
    !2* (1)
    do pp=1,dimp
      b(1:dimq,pp) = a(pp,:)
    end do
  end if

else

  ! nfact = -1

  !do qq=1,dimq
  !  indx(q) = qq
  !  do pp=1,dimp
  !    indx(p) = pp
  !    b(indx(1),indx(2)) = -a(pp,qq)
  !  end do
  !end do

  if (p == 1) then
    !1* (2)
    b(1:dimp,1:dimq) = -a
  else
    !2* (1)
    do pp=1,dimp
      b(1:dimq,pp) = -a(pp,:)
    end do
  end if

end if

return

end subroutine cct3_map22
