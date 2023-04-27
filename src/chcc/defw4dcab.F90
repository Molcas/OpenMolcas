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

subroutine DefW4dcab(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (dc|ab)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(cdLen,abLen)
integer(kind=iwp) :: a, b, c, d, ab, dc

#include "macros.fh"
unused_var(cSGrp)
unused_var(dSGrp)

if (aSGrp == bSGrp) then
  ! case (d,c|a=b)
  do a=1,dima
    ab = nTri_Elem(a-1)
    dc = 0
    do c=1,dimc
      do d=1,dimd
        dc = dc+1
        W(a,1:a-1,c,d) = W(a,1:a-1,c,d)+Wx(dc,ab+1:ab+a-1)
        W(1:a,a,c,d) = W(1:a,a,c,d)+Wx(dc,ab+1:ab+a)
      end do
    end do
  end do

else
  ! case (d,c|a,b)
  ab = 0
  do b=1,dimb
    dc = 0
    do c=1,dimc
      do d=1,dimd
        dc = dc+1
        W(:,b,c,d) = W(:,b,c,d)+Wx(dc,ab+1:ab+dima)
      end do
    end do
  end do

end if

return

end subroutine DefW4dcab
