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

real*8 function GAM(M)

implicit real*8(A-H,O-Z)
#include "crelop.fh"
save/CRELOP/

if (mod(M,2) == 0) goto 10
MA = (M+1)/2
G = 1.d0
if (MA == 1) goto 11
do I=2,MA
  G = G*dble(I-1)
end do
GAM = G
return
10 MA = M
G = SQPI
if (MA == 0) goto 11
do I=1,MA,2
  G = G*0.5d0*dble(I)
end do
11 GAM = G

return

end function GAM
