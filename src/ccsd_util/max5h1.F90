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

subroutine max5h1(imax,rmax,symp,symq,symr,syms,p,q,r,s,val)
! this routine adds new max amplitude and skips smallest
!
! imax - store of indices of 5 max (I/O)
! rmax - store of values of 5 max (I/O)
! symp - symmetry of p index (I)
! symq - symmetry of q index (I)
! symr - symmetry of r index (I)
! syms - symmetry of s index (I)
! p    - value of p index (I)
! q    - value of q index (I)
! r    - value of r index (I)
! s    - value of s index (I)
! val  - value of amplitude (I)

integer imax(1:8,1:5)
real*8 rmax(1:5)
integer symp, symq, symr, syms, p, q, r, s
real*8 val
! help variables
integer nhelp1, nhelp2, nhelp3

!1 find position of this value

do nhelp1=1,5
  if (abs(val) >= abs(rmax(nhelp1))) then
    goto 20
  end if
end do

!2 1-(nhelp1-1) stay untought

!3 push other records if necc.

20 continue
if (nhelp1 < 5) then
  do nhelp2=4,nhelp1,-1
    rmax(nhelp2+1) = rmax(nhelp2)
    do nhelp3=1,8
      imax(nhelp3,nhelp2+1) = imax(nhelp3,nhelp2)
    end do
  end do
end if

!4 add new one

rmax(nhelp1) = val
imax(1,nhelp1) = symp
imax(2,nhelp1) = symq
imax(3,nhelp1) = symr
imax(4,nhelp1) = syms
imax(5,nhelp1) = p
imax(6,nhelp1) = q
imax(7,nhelp1) = r
imax(8,nhelp1) = s

return

end subroutine max5h1
