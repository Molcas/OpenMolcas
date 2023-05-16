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

subroutine DefW4cdba(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (cd|ba)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(cdLen,abLen)
integer(kind=iwp) :: a, ba, c, cd, d

#include "macros.fh"
unused_var(aSGrp)
unused_var(bSGrp)

if (cSGrp == dSGrp) then
  ! case (c=d|b,a)
  ba = 0
  do a=1,dima
    do c=2,dimc
      cd = nTri_Elem(c-1)
      do d=1,c-1
        cd = cd+1
        W(a,:,c,d) = W(a,:,c,d)+Wx(cd,ba+1:ba+dimb)
        W(a,:,d,c) = W(a,:,d,c)+Wx(cd,ba+1:ba+dimb)
      end do
    end do
    do c=1,dimc
      cd = nTri_Elem(c)
      W(a,:,c,c) = W(a,:,c,c)+Wx(cd,ba+1:ba+dimb)
    end do
    ba = ba+dimb
  end do

else
  ! case (c,d|b,a)
  ba = 0
  do a=1,dima
    cd = 0
    do d=1,dimd
      do c=1,dimc
        cd = cd+1
        W(a,:,c,d) = W(a,:,c,d)+Wx(cd,ba+1:ba+dimb)
      end do
    end do
    ba = ba+dimb
  end do

end if

return

end subroutine DefW4cdba
