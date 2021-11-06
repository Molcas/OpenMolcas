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

#include "Constants_mula.fh"
real*8 frequency
integer NormModes(l_NormModes)
real*8 Freq(l_NormModes)
character*(*) Title
character*8 F1
#include "inout.fh"

F1 = '(a2,a)'
NumInt = l_NormModes
write(6,*)
write(6,*)
write(6,*)
write(6,F1) ' ',Title
write(6,F1) ' ','===================================================='
write(6,F1) ' ',' mode        1/cm             GHz          hartrees '
write(6,F1) ' ','----------------------------------------------------'
do i=1,NumInt
  frequency = Freq(i)
  write(6,'(A3,I2,A1,A3,F14.8,F18.8,F12.8)') ' ',NormModes(i),'.',' ',frequency*HarToRcm,frequency*HarToGHz,frequency
end do
write(6,F1) ' ','===================================================='
write(6,*)
write(6,*)

end subroutine WriteFreq
