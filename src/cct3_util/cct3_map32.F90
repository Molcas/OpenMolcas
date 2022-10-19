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
integer(kind=iwp) :: dimp, dimq, dimr, dim_1, dim_2, dim_3, p, q, nfact
real(kind=wp) :: a(dimp,dimq,dimr), b(dim_1,dim_2,dim_3)
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
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,qq,rr) = a(pp,qq,rr)
          end do
        end do
      end do
    else
      !13* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,rr,qq) = a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 2) then
    !2**
    if (q == 1) then
      !21* (3)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,pp,rr) = a(pp,qq,rr)
          end do
        end do
      end do
    else
      !23* (1)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(rr,pp,qq) = a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 3) then
    !3**
    if (q == 1) then
      !31* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,rr,pp) = a(pp,qq,rr)
          end do
        end do
      end do
    else
      !32* (1)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(rr,qq,pp) = a(pp,qq,rr)
          end do
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
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,qq,rr) = -a(pp,qq,rr)
          end do
        end do
      end do
    else
      !13* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,rr,qq) = -a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 2) then
    !2**
    if (q == 1) then
      !21* (3)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,pp,rr) = -a(pp,qq,rr)
          end do
        end do
      end do
    else
      !23* (1)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(rr,pp,qq) = -a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 3) then
    !3**
    if (q == 1) then
      !31* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,rr,pp) = -a(pp,qq,rr)
          end do
        end do
      end do
    else
      !32* (1)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(rr,qq,pp) = -a(pp,qq,rr)
          end do
        end do
      end do
    end if
  end if

end if

return

end subroutine cct3_map32
