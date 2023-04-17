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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, dimc, dimd, abLen, cdLen, aSGrp, bSGrp, cSGrp, dSGrp
real(kind=wp) :: W(dima,dimb,dimc,dimd), Wx(cdLen,abLen)
integer(kind=iwp) :: a, b, ba, c, cd, d

#include "macros.fh"
unused_var(aSGrp)
unused_var(bSGrp)

if (cSGrp == dSGrp) then
  ! case (c=d|b,a)
  ba = 0
  do a=1,dima
    do b=1,dimb
      ba = ba+1
      do c=2,dimc
        cd = c*(c-1)/2
        do d=1,c-1
          cd = cd+1
          W(a,b,c,d) = W(a,b,c,d)+Wx(cd,ba)
          W(a,b,d,c) = W(a,b,d,c)+Wx(cd,ba)
        end do
      end do
      do c=1,dimc
        cd = c*(c+1)/2
        W(a,b,c,c) = W(a,b,c,c)+Wx(cd,ba)
      end do
    end do
  end do

else
  ! case (c,d|b,a)
  ba = 0
  do a=1,dima
    do b=1,dimb
      ba = ba+1
      cd = 0
      do d=1,dimd
        do c=1,dimc
          cd = cd+1
          W(a,b,c,d) = W(a,b,c,d)+Wx(cd,ba)
        end do
      end do
    end do
  end do

end if

return

end subroutine DefW4cdba
