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

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Dens(*), DenCorr(*), Trace_Diff
integer(kind=iwp) :: iOrb, iOcc
integer(kind=iwp) :: i, j, kaunt
real(kind=wp) :: T, Trace_HF

Trace_HF = real(iOcc*2,kind=wp)
kaunt = 0
T = Trace_HF/(Trace_HF-Trace_Diff)
do i=1,iOrb
  do j=1,i
    kaunt = kaunt+1
    Dens(kaunt) = T*(Dens(kaunt)-DenCorr(kaunt))
  end do
end do
!Trace = Zero
!kaunt = 0
!do i=1,iOrb
!  do j=1,i
!    kaunt = kaunt+1
!    if (i == j) Trace = Trace+Dens(kaunt)
!  end do
!end do
!call triprt('KKK',' ',Dens,iorb)
!write(u6,*) 'QQQ:',Trace
return

end subroutine DCorrCorr
