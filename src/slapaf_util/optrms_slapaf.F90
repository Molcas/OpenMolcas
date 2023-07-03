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

subroutine OptRMS_Slapaf(x,y,nAt,RMS,RMSMax)

use Symmetry_Info, only: VarR, VarT
implicit real*8(a-h,o-z)
#include "print.fh"
#include "real.fh"
real*8 x(3,nAt), y(3,nAt)

! Only align if energy is rotational and translational invariant.
! (no weighting)

if (.not. (VarR .or. VarT)) then
  call Superpose(x,y,nAt,RMS,RMSMax)
else

  ! Otherwise, just compute RMS

  RMS = Zero
  RMSMax = Zero
  do i=1,nAt
    disp = Zero
    do ixyz=1,3
      diff = x(ixyz,i)-y(ixyz,i)
      RMS = RMS+diff**2
      disp = disp+diff**2
    end do
    if (sqrt(disp) > RMSMax) RMSMax = sqrt(disp)
  end do
  RMS = sqrt(RMS/dble(nAt))
end if

return

end subroutine OptRMS_Slapaf
