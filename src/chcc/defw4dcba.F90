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

subroutine DefW4dcba(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (dc|ba)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(cdLen,abLen)
integer(kind=iwp) :: a, ba, c, d, dc

#include "macros.fh"
unused_var(aSGrp)
unused_var(bSGrp)
unused_var(cSGrp)
unused_var(dSGrp)

! case (c,d|b,a)
ba = 0
do a=1,dima
  dc = 0
  do c=1,dimc
    do d=1,dimd
      dc = dc+1
      W(a,:,c,d) = W(a,:,c,d)+Wx(dc,ba+1:ba+dimb)
    end do
  end do
  ba = ba+dimb
end do

return

end subroutine DefW4dcba
