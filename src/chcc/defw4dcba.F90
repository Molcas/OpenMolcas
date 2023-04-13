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
integer(kind=iwp) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp) :: W(dima,dimb,dimc,dimd), Wx(cdLen,abLen)
integer(kind=iwp) :: a, b, ba, c, d, dc

! case (c,d|b,a)
ba = 0
do a=1,dima
  do b=1,dimb
    ba = ba+1
    dc = 0
    do c=1,dimc
      do d=1,dimd
        dc = dc+1
        W(a,b,c,d) = W(a,b,c,d)+Wx(dc,ba)
      end do
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(aSGrp)
  call Unused_integer(bSGrp)
  call Unused_integer(cSGrp)
  call Unused_integer(dSGrp)
end if

end subroutine DefW4dcba
