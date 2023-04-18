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

subroutine DefW4cdab(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (cd|ab)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp) :: W(dima,dimb,dimc,dimd), Wx(cdLen,abLen)
integer(kind=iwp) :: a, ab, b, c, cd, d

if ((aSGrp == bSGrp) .and. (cSGrp == dSGrp)) then
  ! case (c=d|a=b)
  do a=1,dima
    ab = a*(a-1)/2
    do c=2,dimc
      cd = c*(c-1)/2
      do d=1,c-1
        cd = cd+1
        W(a,1:a-1,c,d) = W(a,1:a-1,c,d)+Wx(cd,ab+1:ab+a-1)
        W(a,1:a-1,d,c) = W(a,1:a-1,d,c)+Wx(cd,ab+1:ab+a-1)
        W(1:a,a,c,d) = W(1:a,a,c,d)+Wx(cd,ab+1:ab+a)
        W(1:a,a,d,c) = W(1:a,a,d,c)+Wx(cd,ab+1:ab+a)
      end do
    end do
    do c=1,dimc
      cd = c*(c+1)/2
      W(a,1:a-1,c,c) = W(a,1:a-1,c,c)+Wx(cd,ab+1:ab+a-1)
      W(1:a,a,c,c) = W(1:a,a,c,c)+Wx(cd,ab+1:ab+a)
    end do
  end do

else if ((aSGrp == bSGrp) .and. (cSGrp /= dSGrp)) then
  ! case (c,d|a=b)
  do a=1,dima
    ab = a*(a-1)/2
    cd = 0
    do d=1,dimd
      do c=1,dimc
        cd = cd+1
        W(a,1:a-1,c,d) = W(a,1:a-1,c,d)+Wx(cd,ab+1:ab+a-1)
        W(1:a,a,c,d) = W(1:a,a,c,d)+Wx(cd,ab+1:ab+a)
      end do
    end do
  end do

else if ((aSGrp /= bSGrp) .and. (cSGrp == dSGrp)) then
  ! case (c=d|a,b)
  ab = 0
  do b=1,dimb
    do c=2,dimc
      cd = c*(c-1)/2
      do d=1,c-1
        cd = cd+1
        W(:,b,c,d) = W(:,b,c,d)+Wx(cd,ab+1:ab+dima)
        W(:,b,d,c) = W(:,b,d,c)+Wx(cd,ab+1:ab+dima)
      end do
    end do
    do c=1,dimc
      cd = c*(c+1)/2
      W(:,b,c,c) = W(:,b,c,c)+Wx(cd,ab+1:ab+dima)
    end do
    ab = ab+dima
  end do

else
  ! case (c,d|a,b)
  ab = 0
  do b=1,dimb
    cd = 0
    do d=1,dimd
      do c=1,dimc
        cd = cd+1
        W(:,b,c,d) = W(:,b,c,d)+Wx(cd,ab+1:ab+dima)
      end do
    end do
    ab = ab+dima
  end do

end if

return

end subroutine DefW4cdab
