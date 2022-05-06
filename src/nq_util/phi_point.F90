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

subroutine Phi_point(iPhi,nPhi,Cos_Phi,Sin_Phi,w_Phi)

implicit real*8(a-h,o-z)
#include "real.fh"

q = Pi*(Two*dble(iPhi)-1.0d0)/dble(nPhi)
if (abs(cos(q)) > 1.0D-14) then
  Cos_Phi = cos(q)
else
  Cos_Phi = Zero
end if
if (abs(sin(q)) > 1.0D-14) then
  Sin_Phi = sin(q)
else
  Sin_Phi = Zero
end if
w_Phi = Two*Pi/dble(nPhi)

return

end subroutine Phi_point
