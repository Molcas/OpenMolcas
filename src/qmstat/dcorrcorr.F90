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

! MP2 density correction.
subroutine DCorrCorr(Dens,DenCorr,Trace_Diff,iOrb,iOcc)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iOrb, iOcc
real(kind=wp), intent(inout) :: Dens(nTri_Elem(iOrb))
real(kind=wp), intent(in) :: DenCorr(nTri_Elem(iOrb)), Trace_Diff
real(kind=wp) :: T, Trace_HF

Trace_HF = real(iOcc*2,kind=wp)
T = Trace_HF/(Trace_HF-Trace_Diff)
Dens(:) = T*(Dens(:)-DenCorr(:))
!Trace = Zero
!kaunt = 0
!do i=1,iOrb
!  Trace = Trace+Dens(nTri_Elem(i))
!end do
!call triprt('KKK',' ',Dens,iorb)
!write(u6,*) 'QQQ:',Trace
return

end subroutine DCorrCorr
