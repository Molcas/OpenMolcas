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

subroutine DefW4abcd(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (ab|cd)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(abLen,cdLen)
integer(kind=iwp) :: a, ab, b, c, cd, d

if ((aSGrp == bSGrp) .and. (cSGrp == dSGrp)) then
  ! case (a=b|c=d)
  do c=1,dimc
    cd = nTri_Elem(c-1)
    do a=1,dima
      ab = nTri_Elem(a-1)
      W(a,1:a-1,c,1:c-1) = W(a,1:a-1,c,1:c-1)+Wx(ab+1:ab+a-1,cd+1:cd+c-1)
      W(a,1:a-1,1:c,c) = W(a,1:a-1,1:c,c)+Wx(ab+1:ab+a-1,cd+1:cd+c)
      W(1:a,a,c,1:c-1) = W(1:a,a,c,1:c-1)+Wx(ab+1:ab+a,cd+1:cd+c-1)
      W(1:a,a,1:c,c) = W(1:a,a,1:c,c)+Wx(ab+1:ab+a,cd+1:cd+c)
      ab = ab+a
    end do
    cd = cd+c
  end do

else if ((aSGrp == bSGrp) .and. (cSGrp /= dSGrp)) then
  ! case (a=b|c,d)
  cd = 0
  do d=1,dimd
    do a=1,dima
      ab = nTri_Elem(a-1)
      W(a,1:a-1,:,d) = W(a,1:a-1,:,d)+Wx(ab+1:ab+a-1,cd+1:cd+dimc)
      W(1:a,a,:,d) = W(1:a,a,:,d)+Wx(ab+1:ab+a,cd+1:cd+dimc)
      ab = ab+a
    end do
    cd = cd+dimc
  end do

else if ((aSGrp /= bSGrp) .and. (cSGrp == dSGrp)) then
  ! case (a,b|c=d)
  do c=1,dimc
    cd = nTri_Elem(c-1)
    ab = 0
    do b=1,dimb
      W(:,b,c,1:c-1) = W(:,b,c,1:c-1)+Wx(ab+1:ab+dima,cd+1:cd+c-1)
      W(:,b,1:c,c) = W(:,b,1:c,c)+Wx(ab+1:ab+dima,cd+1:cd+c)
      ab = ab+dima
    end do
    cd = cd+c
  end do

else
  ! case (a,b|c,d)
  cd = 0
  do d=1,dimd
    ab = 0
    do b=1,dimb
      W(:,b,:,d) = W(:,b,:,d)+Wx(ab+1:ab+dima,cd+1:cd+dimc)
      ab = ab+dima
    end do
    cd = cd+dimc
  end do

end if

return

end subroutine DefW4abcd
