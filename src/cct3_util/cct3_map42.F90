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
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,qq,rr,ss) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !124* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,qq,ss,rr) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 3) then
      !13**
      if (r == 2) then
        !132* (4)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,rr,qq,ss) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !134* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,ss,qq,rr) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 4) then
      !14**
      if (r == 2) then
        !142* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,rr,ss,qq) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !143* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,ss,rr,qq) = a(pp,qq,rr,ss)
              end do
            end do
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
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,pp,rr,ss) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !214* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,pp,ss,rr) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 3) then
      !23**
      if (r == 1) then
        !231* (4)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,pp,qq,ss) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !234* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,pp,qq,rr) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 4) then
      !24**
      if (r == 3) then
        !243* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,pp,rr,qq) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !241* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,pp,ss,qq) = a(pp,qq,rr,ss)
              end do
            end do
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
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,rr,pp,ss) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !314* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,ss,pp,rr) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 2) then
      !32**
      if (r == 1) then
        !321* (4)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,qq,pp,ss) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !324* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,qq,pp,rr) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 4) then
      !34**
      if (r == 1) then
        !341* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,ss,pp,qq) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !342* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,rr,pp,qq) = a(pp,qq,rr,ss)
              end do
            end do
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
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,ss,rr,pp) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !412* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,rr,ss,pp) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 2) then
      !42**
      if (r == 1) then
        !421* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,qq,ss,pp) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !423* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,qq,rr,pp) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 3) then
      !43**
      if (r == 1) then
        !431* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,ss,qq,pp) = a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !432* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,rr,qq,pp) = a(pp,qq,rr,ss)
              end do
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
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,qq,rr,ss) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !124* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,qq,ss,rr) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 3) then
      !13**
      if (r == 2) then
        !132* (4)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,rr,qq,ss) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !134* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,ss,qq,rr) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 4) then
      !14**
      if (r == 2) then
        !142* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,rr,ss,qq) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !143* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(pp,ss,rr,qq) = -a(pp,qq,rr,ss)
              end do
            end do
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
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,pp,rr,ss) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !214* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,pp,ss,rr) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 3) then
      !23**
      if (r == 1) then
        !231* (4)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,pp,qq,ss) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !234* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,pp,qq,rr) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 4) then
      !24**
      if (r == 3) then
        !243* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,pp,rr,qq) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !241* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,pp,ss,qq) = -a(pp,qq,rr,ss)
              end do
            end do
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
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,rr,pp,ss) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !314* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,ss,pp,rr) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 2) then
      !32**
      if (r == 1) then
        !321* (4)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,qq,pp,ss) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !324* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,qq,pp,rr) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 4) then
      !34**
      if (r == 1) then
        !341* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,ss,pp,qq) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !342* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,rr,pp,qq) = -a(pp,qq,rr,ss)
              end do
            end do
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
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,ss,rr,pp) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !412* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(qq,rr,ss,pp) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 2) then
      !42**
      if (r == 1) then
        !421* (3)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,qq,ss,pp) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !423* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,qq,rr,pp) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    else if (q == 3) then
      !43**
      if (r == 1) then
        !431* (2)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(rr,ss,qq,pp) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      else
        !432* (1)
        do ss=1,dims
          do rr=1,dimr
            do qq=1,dimq
              do pp=1,dimp
                b(ss,rr,qq,pp) = -a(pp,qq,rr,ss)
              end do
            end do
          end do
        end do
      end if
    end if

  end if

end if

return

end subroutine cct3_map42
