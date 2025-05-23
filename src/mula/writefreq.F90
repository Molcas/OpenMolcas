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

subroutine WriteFreq(Freq,NormModes,l_NormModes,Title)

use Constants, only: auTocm, auToHz
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: l_NormModes, NormModes(l_NormModes)
real(kind=wp), intent(in) :: Freq(l_NormModes)
character(len=*), intent(in) :: Title
integer(kind=iwp) :: i, NumInt
real(kind=wp) :: frequency

NumInt = l_NormModes
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,'(2x,a)') Title
write(u6,'(2x,a)') '===================================================='
write(u6,'(2x,a)') ' mode        1/cm             GHz          hartrees'
write(u6,'(2x,a)') '----------------------------------------------------'
do i=1,NumInt
  frequency = Freq(i)
  write(u6,'(3X,I2,A1,3X,F14.8,F18.8,F12.8)') NormModes(i),'.',frequency*auTocm,frequency*auToHz*1.0e-6_wp,frequency
end do
write(u6,'(2x,a)') '===================================================='
write(u6,*)
write(u6,*)

end subroutine WriteFreq
