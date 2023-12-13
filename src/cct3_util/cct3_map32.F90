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

subroutine cct3_map32(a,b,dimp,dimq,dimr,dim_1,dim_2,dim_3,p,q,nfact)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr, dim_1, dim_2, dim_3, p, q, nfact
real(kind=wp), intent(in) :: a(dimp,dimq,dimr)
real(kind=wp), intent(out) :: b(dim_1,dim_2,dim_3)
integer(kind=iwp) :: pp, qq, rr !, indx(3)

if (nfact == 1) then

  ! fact = + 1

  !do rr=1,dimr
  !  indx(r) = rr
  !  do qq=1,dimq
  !    indx(q) = qq
  !    do pp=1,dimp
  !      indx(p) = pp
  !      b(indx(1),indx(2),indx(3)) = a(pp,qq,rr)
  !    end do
  !  end do
  !end do

  if (p == 1) then
    !1**
    if (q == 2) then
      !12* (3)
      b(1:dimp,1:dimq,1:dimr) = a
    else
      !13* (2)
      do qq=1,dimq
        b(1:dimp,1:dimr,qq) = a(:,qq,:)
      end do
    end if
  else if (p == 2) then
    !2**
    if (q == 1) then
      !21* (3)
      do pp=1,dimp
        b(1:dimq,pp,1:dimr) = a(pp,:,:)
      end do
    else
      !23* (1)
      do rr=1,dimr
        b(rr,1:dimp,1:dimq) = a(:,:,rr)
      end do
    end if
  else if (p == 3) then
    !3**
    if (q == 1) then
      !31* (2)
      do pp=1,dimp
        b(1:dimq,1:dimr,pp) = a(pp,:,:)
      end do
    else
      !32* (1)
      do qq=1,dimq
        do pp=1,dimp
          b(1:dimr,qq,pp) = a(pp,qq,:)
        end do
      end do
    end if
  end if

else

  ! factor = -1

  !do rr=1,dimr
  !  indx(r) = rr
  !  do qq=1,dimq
  !    indx(q) = qq
  !    do pp=1,dimp
  !      indx(p) = pp
  !      b(indx(1),indx(2),indx(3)) = -a(pp,qq,rr)
  !    end do
  !  end do
  !end do

  if (p == 1) then
    !1**
    if (q == 2) then
      !12* (3)
      b(1:dimp,1:dimq,1:dimr) = -a
    else
      !13* (2)
      do qq=1,dimq
        b(1:dimp,1:dimr,qq) = -a(:,qq,:)
      end do
    end if
  else if (p == 2) then
    !2**
    if (q == 1) then
      !21* (3)
      do pp=1,dimp
        b(1:dimq,pp,1:dimr) = -a(pp,:,:)
      end do
    else
      !23* (1)
      do rr=1,dimr
        b(rr,1:dimp,1:dimq) = -a(:,:,rr)
      end do
    end if
  else if (p == 3) then
    !3**
    if (q == 1) then
      !31* (2)
      do pp=1,dimp
        b(1:dimq,1:dimr,pp) = -a(pp,:,:)
      end do
    else
      !32* (1)
      do qq=1,dimq
        do pp=1,dimp
          b(1:dimr,qq,pp) = -a(pp,qq,:)
        end do
      end do
    end if
  end if

end if

return

end subroutine cct3_map32
