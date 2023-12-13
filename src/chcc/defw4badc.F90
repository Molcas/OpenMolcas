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

subroutine DefW4badc(W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp)
! define W(a,b,c,d) from (ba|dc)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp), intent(inout) :: W(dima,dimb,dimc,dimd)
real(kind=wp), intent(in) :: Wx(abLen,cdLen)
integer(kind=iwp) :: a, ba, c, dc

#include "macros.fh"
unused_var(aSGrp)
unused_var(bSGrp)
unused_var(cSGrp)
unused_var(dSGrp)

! case (b,a|c,d)
dc = 0
do c=1,dimc
  ba = 0
  do a=1,dima
    W(a,:,c,:) = W(a,:,c,:)+Wx(ba+1:ba+dimb,dc+1:dc+dimd)
    ba = ba+dimb
  end do
  dc = dc+dimd
end do

return

end subroutine DefW4badc
