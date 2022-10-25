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

subroutine cct3_map42(a,b,dimp,dimq,dimr,dims,dim_1,dim_2,dim_3,dim_4,p,q,r,nfact)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr, dims, dim_1, dim_2, dim_3, dim_4, p, q, r, nfact
real(kind=wp), intent(in) :: a(dimp,dimq,dimr,dims)
real(kind=wp), intent(out) :: b(dim_1,dim_2,dim_3,dim_4)
integer(kind=iwp) :: pp, qq, rr, ss !, indx(4)

if (nfact == 1) then

  ! factor + 1

  !do ss=1,dims
  !  indx(s) = ss
  !  do rr=1,dimr
  !    indx(r) = rr
  !    do qq=1,dimq
  !      indx(q) = qq
  !      do pp=1,dimp
  !        indx(p) = pp
  !        b(indx(1),indx(2),indx(3),indx(4)) = a(pp,qq,rr,ss)
  !      end do
  !    end do
  !  end do
  !end do

  if (p == 1) then
    !1***
    if (q == 2) then
      !12**
      if (r == 3) then
        !123* (4)
        b(1:dimp,1:dimq,1:dimr,1:dims) = a
      else
        !124* (3)
        do rr=1,dimr
          b(1:dimp,1:dimq,1:dims,rr) = a(:,:,rr,:)
        end do
      end if
    else if (q == 3) then
      !13**
      if (r == 2) then
        !132* (4)
        do qq=1,dimq
          b(1:dimp,1:dimr,qq,1:dims) = a(:,qq,:,:)
        end do
      else
        !134* (2)
        do ss=1,dims
          b(1:dimp,ss,1:dimq,1:dimr) = a(:,:,:,ss)
        end do
      end if
    else if (q == 4) then
      !14**
      if (r == 2) then
        !142* (3)
        do qq=1,dimq
          b(1:dimp,1:dimr,1:dims,qq) = a(:,qq,:,:)
        end do
      else
        !143* (2)
        do ss=1,dims
          do rr=1,dimr
            b(1:dimp,ss,rr,1:dimq) = a(:,:,rr,ss)
          end do
        end do
      end if
    end if

  else if (p == 2) then
    !2***
    if (q == 1) then
      !21**
      if (r == 3) then
        !213* (4)
        do pp=1,dimp
          b(1:dimq,pp,1:dimr,1:dims) = a(pp,:,:,:)
        end do
      else
        !214* (3)
        do rr=1,dimr
          do pp=1,dimp
            b(1:dimq,pp,1:dims,rr) = a(pp,:,rr,:)
          end do
        end do
      end if
    else if (q == 3) then
      !23**
      if (r == 1) then
        !231* (4)
        do rr=1,dimr
          b(rr,1:dimp,1:dimq,1:dims) = a(:,:,rr,:)
        end do
      else
        !234* (1)
        do ss=1,dims
          b(ss,1:dimp,1:dimq,1:dimr) = a(:,:,:,ss)
        end do
      end if
    else if (q == 4) then
      !24**
      if (r == 3) then
        !243* (1)
        do ss=1,dims
          do qq=1,dimq
            b(ss,1:dimp,1:dimr,qq) = a(:,qq,:,ss)
          end do
        end do
      else
        !241* (3)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,pp,1:dims,qq) = a(pp,qq,:,:)
          end do
        end do
      end if
    end if

  else if (p == 3) then
    !3***
    if (q == 1) then
      !31**
      if (r == 2) then
        !312* (4)
        do pp=1,dimp
          b(1:dimq,1:dimr,pp,1:dims) = a(pp,:,:,:)
        end do
      else
        !314* (2)
        do rr=1,dimr
          do pp=1,dimp
            b(1:dimq,1:dims,pp,rr) = a(pp,:,rr,:)
          end do
        end do
      end if
    else if (q == 2) then
      !32**
      if (r == 1) then
        !321* (4)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,qq,pp,1:dims) = a(pp,qq,:,:)
          end do
        end do
      else
        !324* (1)
        do ss=1,dims
          do pp=1,dimp
            b(ss,1:dimq,pp,1:dimr) = a(pp,:,:,ss)
          end do
        end do
      end if
    else if (q == 4) then
      !34**
      if (r == 1) then
        !341* (2)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,1:dims,pp,qq) = a(pp,qq,:,:)
          end do
        end do
      else
        !342* (1)
        do ss=1,dims
          do rr=1,dimr
            b(ss,rr,1:dimp,1:dimq) = a(:,:,rr,ss)
          end do
        end do
      end if
    end if

  else if (p == 4) then
    !4***
    if (q == 1) then
      !41**
      if (r == 3) then
        !413* (2)
        do ss=1,dims
          do pp=1,dimp
            b(1:dimq,ss,1:dimr,pp) = a(pp,:,:,ss)
          end do
        end do
      else
        !412* (3)
        do pp=1,dimp
          b(1:dimq,1:dimr,1:dims,pp) = a(pp,:,:,:)
        end do
      end if
    else if (q == 2) then
      !42**
      if (r == 1) then
        !421* (3)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,qq,1:dims,pp) = a(pp,qq,:,:)
          end do
        end do
      else
        !423* (1)
        do ss=1,dims
          do pp=1,dimp
            b(ss,1:dimq,1:dimr,pp) = a(pp,:,:,ss)
          end do
        end do
      end if
    else if (q == 3) then
      !43**
      if (r == 1) then
        !431* (2)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,1:dims,qq,pp) = a(pp,qq,:,:)
          end do
        end do
      else
        !432* (1)
        do rr=1,dimr
          do qq=1,dimq
            do pp=1,dimp
              b(1:dims,rr,qq,pp) = a(pp,qq,rr,:)
            end do
          end do
        end do
      end if
    end if

  end if

else

  ! factor = -1

  !do ss=1,dims
  !  indx(s) = ss
  !  do rr=1,dimr
  !    indx(r) = rr
  !    do qq=1,dimq
  !      indx(q) = qq
  !      do pp=1,dimp
  !        indx(p) = pp
  !        b(indx(1),indx(2),indx(3),indx(4)) = -a(pp,qq,rr,ss)
  !      end do
  !    end do
  !  end do
  !end do

  if (p == 1) then
    !1***
    if (q == 2) then
      !12**
      if (r == 3) then
        !123* (4)
        b(1:dimp,1:dimq,1:dimr,1:dims) = -a
      else
        !124* (3)
        do rr=1,dimr
          b(1:dimp,1:dimq,1:dims,rr) = -a(:,:,rr,:)
        end do
      end if
    else if (q == 3) then
      !13**
      if (r == 2) then
        !132* (4)
        do qq=1,dimq
          b(1:dimp,1:dimr,qq,1:dims) = -a(:,qq,:,:)
        end do
      else
        !134* (2)
        do ss=1,dims
          b(1:dimp,ss,1:dimq,1:dimr) = -a(:,:,:,ss)
        end do
      end if
    else if (q == 4) then
      !14**
      if (r == 2) then
        !142* (3)
        do qq=1,dimq
          b(1:dimp,1:dimr,1:dims,qq) = -a(:,qq,:,:)
        end do
      else
        !143* (2)
        do ss=1,dims
          do rr=1,dimr
            b(1:dimp,ss,rr,1:dimq) = -a(:,:,rr,ss)
          end do
        end do
      end if
    end if

  else if (p == 2) then
    !2***
    if (q == 1) then
      !21**
      if (r == 3) then
        !213* (4)
        do pp=1,dimp
          b(1:dimq,pp,1:dimr,1:dims) = -a(pp,:,:,:)
        end do
      else
        !214* (3)
        do rr=1,dimr
          do pp=1,dimp
            b(1:dimq,pp,1:dims,rr) = -a(pp,:,rr,:)
          end do
        end do
      end if
    else if (q == 3) then
      !23**
      if (r == 1) then
        !231* (4)
        do rr=1,dimr
          b(rr,1:dimp,1:dimq,1:dims) = -a(:,:,rr,:)
        end do
      else
        !234* (1)
        do ss=1,dims
          b(ss,1:dimp,1:dimq,1:dimr) = -a(:,:,:,ss)
        end do
      end if
    else if (q == 4) then
      !24**
      if (r == 3) then
        !243* (1)
        do ss=1,dims
          do qq=1,dimq
            b(ss,1:dimp,1:dimr,qq) = -a(:,qq,:,ss)
          end do
        end do
      else
        !241* (3)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,pp,1:dims,qq) = -a(pp,qq,:,:)
          end do
        end do
      end if
    end if

  else if (p == 3) then
    !3***
    if (q == 1) then
      !31**
      if (r == 2) then
        !312* (4)
        do pp=1,dimp
          b(1:dimq,1:dimr,pp,1:dims) = -a(pp,:,:,:)
        end do
      else
        !314* (2)
        do rr=1,dimr
          do pp=1,dimp
            b(1:dimq,1:dims,pp,rr) = -a(pp,:,rr,:)
          end do
        end do
      end if
    else if (q == 2) then
      !32**
      if (r == 1) then
        !321* (4)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,qq,pp,1:dims) = -a(pp,qq,:,:)
          end do
        end do
      else
        !324* (1)
        do ss=1,dims
          do pp=1,dimp
            b(ss,1:dimq,pp,1:dimr) = -a(pp,:,:,ss)
          end do
        end do
      end if
    else if (q == 4) then
      !34**
      if (r == 1) then
        !341* (2)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,1:dims,pp,qq) = -a(pp,qq,:,:)
          end do
        end do
      else
        !342* (1)
        do ss=1,dims
          do rr=1,dimr
            b(ss,rr,1:dimp,1:dimq) = -a(:,:,rr,ss)
          end do
        end do
      end if
    end if

  else if (p == 4) then
    !4***
    if (q == 1) then
      !41**
      if (r == 3) then
        !413* (2)
        do ss=1,dims
          do pp=1,dimp
            b(1:dimq,ss,1:dimr,pp) = -a(pp,:,:,ss)
          end do
        end do
      else
        !412* (3)
        do pp=1,dimp
          b(1:dimq,1:dimr,1:dims,pp) = -a(pp,:,:,:)
        end do
      end if
    else if (q == 2) then
      !42**
      if (r == 1) then
        !421* (3)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,qq,1:dims,pp) = -a(pp,qq,:,:)
          end do
        end do
      else
        !423* (1)
        do ss=1,dims
          do pp=1,dimp
            b(ss,1:dimq,1:dimr,pp) = -a(pp,:,:,ss)
          end do
        end do
      end if
    else if (q == 3) then
      !43**
      if (r == 1) then
        !431* (2)
        do qq=1,dimq
          do pp=1,dimp
            b(1:dimr,1:dims,qq,pp) = -a(pp,qq,:,:)
          end do
        end do
      else
        !432* (1)
        do rr=1,dimr
          do qq=1,dimq
            do pp=1,dimp
              b(1:dims,rr,qq,pp) = -a(pp,qq,rr,:)
            end do
          end do
        end do
      end if
    end if

  end if

end if

return

end subroutine cct3_map42
