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

subroutine DefW4abdc(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (ab|dc)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(abLen,cdLen)
integer(kind=iwp) :: a, ab, b, c, dc

#include "macros.fh"
unused_var(cSGrp)
unused_var(dSGrp)

if (aSGrp == bSGrp) then
  ! case (a=b|d,c)
  dc = 0
  do c=1,dimc
    do a=1,dima
      ab = nTri_Elem(a-1)
      W(a,1:a-1,c,:) = W(a,1:a-1,c,:)+Wx(ab+1:ab+a-1,dc+1:dc+dimd)
      W(1:a,a,c,:) = W(1:a,a,c,:)+Wx(ab+1:ab+a,dc+1:dc+dimd)
      ab = ab+a
    end do
    dc = dc+dimd
  end do

else
  ! case (a,b|d,c)
  dc = 0
  do c=1,dimc
    ab = 0
    do b=1,dimb
      W(:,b,c,:) = W(:,b,c,:)+Wx(ab+1:ab+dima,dc+1:dc+dimd)
      ab = ab+dima
    end do
    dc = dc+dimd
  end do

end if

return

end subroutine DefW4abdc
