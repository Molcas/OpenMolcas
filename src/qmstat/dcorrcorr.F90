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

implicit real*8(a-h,o-z)
dimension Dens(*), DenCorr(*)

Trace_HF = dble(iOcc*2)
kaunt = 0
T = Trace_HF/(Trace_HF-Trace_Diff)
do i=1,iOrb
  do j=1,i
    kaunt = kaunt+1
    Dens(kaunt) = T*(Dens(kaunt)-DenCorr(kaunt))
  end do
end do
!Trace = 0.0d0
!kaunt = 0
!do i=1,iOrb
!  do j=1,i
!    kaunt = kaunt+1
!    if (i == j) Trace = Trace+Dens(kaunt)
!  end do
!end do
!call triprt('KKK',' ',Dens,iorb)
!write(6,*) 'QQQ:',Trace
return

end subroutine DCorrCorr
