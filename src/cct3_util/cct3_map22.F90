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

subroutine cct3_map22(a,b,dimp,dimq,dim1,dim2,p,q,nfact)

integer dimp, dimq, dim1, dim2, p, q, nfact
real*8 a(1:dimp,1:dimq)
real*8 b(1:dim1,1:dim2)
!integer index(1:2)
integer pp, qq

if (nfact == 1) then

  ! nfact = + 1

  !do qq=1,dimq
  !  index(q) = qq
  !  do pp=1,dimp
  !    index(p) = pp
  !    b(index(1),index(2)) = a(pp,qq)
  !  end do
  !end do

  if (p == 1) then
    !1* (2)
    do qq=1,dimq
      do pp=1,dimp
        b(pp,qq) = a(pp,qq)
      end do
    end do
  else
    !2* (1)
    do qq=1,dimq
      do pp=1,dimp
        b(qq,pp) = a(pp,qq)
      end do
    end do
  end if

else

  ! nfact = -1

  !do qq=1,dimq
  !  index(q) = qq
  !  do pp=1,dimp
  !    index(p) = pp
  !    b(index(1),index(2)) = -a(pp,qq)
  !  end do
  !end do

  if (p == 1) then
    !1* (2)
    do qq=1,dimq
      do pp=1,dimp
        b(pp,qq) = -a(pp,qq)
      end do
    end do
  else
    !2* (1)
    do qq=1,dimq
      do pp=1,dimp
        b(qq,pp) = -a(pp,qq)
      end do
    end do
  end if

end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(q)

end subroutine cct3_map22
