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

subroutine map32(a,b,dimp,dimq,dimr,dim_1,dim_2,dim_3,p,q,r,nfact)

integer dimp, dimq, dimr, dim_1, dim_2, dim_3, p, q, r, nfact
real*8 a(1:dimp,1:dimq,1:dimr)
real*8 b(1:dim_1,1:dim_2,1:dim_3)
!integer index(1:3)
integer pp, qq, rr

if (nfact == 1) then

  ! nfact = + 1

  !do rr=1,dimr
  !  index(r) = rr
  !  do qq=1,dimq
  !    index(q) = qq
  !    do pp=1,dimp
  !      index(p) = pp
  !      b(index(1),index(2),index(3)) = a(pp,qq,rr)
  !    end do
  !  end do
  !end do

  if (p == 1) then
    ! 1**
    if (q == 2) then
      ! 12* (3)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,qq,rr) = a(pp,qq,rr)
          end do
        end do
      end do
    else
      ! 13* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,rr,qq) = a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 2) then
    ! 2**
    if (q == 1) then
      ! 21* (3)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,pp,rr) = a(pp,qq,rr)
          end do
        end do
      end do
    else
      ! 23* (1)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(rr,pp,qq) = a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 3) then
    ! 3**
    if (q == 1) then
      ! 31* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,rr,pp) = a(pp,qq,rr)
          end do
        end do
      end do
    else
      ! 32* (1)
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
  !  index(r) = rr
  !  do qq=1,dimq
  !    index(q) = qq
  !    do pp=1,dimp
  !      index(p) = pp
  !      b(index(1),index(2),index(3)) = -a(pp,qq,rr)
  !    end do
  !  end do
  !end do

  if (p == 1) then
    ! 1**
    if (q == 2) then
      ! 12* (3)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,qq,rr) = -a(pp,qq,rr)
          end do
        end do
      end do
    else
      ! 13* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(pp,rr,qq) = -a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 2) then
    ! 2**
    if (q == 1) then
      ! 21* (3)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,pp,rr) = -a(pp,qq,rr)
          end do
        end do
      end do
    else
      ! 23* (1)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(rr,pp,qq) = -a(pp,qq,rr)
          end do
        end do
      end do
    end if
  else if (p == 3) then
    ! 3**
    if (q == 1) then
      ! 31* (2)
      do rr=1,dimr
        do qq=1,dimq
          do pp=1,dimp
            b(qq,rr,pp) = -a(pp,qq,rr)
          end do
        end do
      end do
    else
      ! 32* (1)
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
! Avoid unused argument warnings
if (.false.) call Unused_integer(r)

end subroutine map32
