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

subroutine DefW4bacd(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (ba|cd)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(abLen,cdLen)
integer(kind=iwp) :: a, ba, c, cd, d

#include "macros.fh"
unused_var(aSGrp)
unused_var(bSGrp)

if (cSGrp == dSGrp) then
  ! case (b,a|c=d)
  do c=1,dimc
    cd = nTri_Elem(c-1)
    ba = 0
    do a=1,dima
      W(a,:,c,1:c-1) = W(a,:,c,1:c-1)+Wx(ba+1:ba+dimb,cd+1:cd+c-1)
      W(a,:,1:c,c) = W(a,:,1:c,c)+Wx(ba+1:ba+dimb,cd+1:cd+c)
      ba = ba+dimb
    end do
    cd = cd+c
  end do

else
  ! case (b,a|c,d)
  cd = 0
  do d=1,dimd
    ba = 0
    do a=1,dima
      W(a,:,:,d) = W(a,:,:,d)+Wx(ba+1:ba+dimb,cd+1:cd+dimc)
      ba = ba+dimb
    end do
    cd = cd+dimc
  end do

end if

return

end subroutine DefW4bacd
