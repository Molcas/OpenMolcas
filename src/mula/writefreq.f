************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine WriteFreq(Freq,NormModes,l_NormModes,Title)
C!
#include "Constants_mula.fh"
      Real*8      frequency
      Integer NormModes(l_NormModes)
      Real*8 Freq(l_NormModes)
      Character*(*) Title
      Character*8 F1
#include "inout.fh"
C!
      F1='(a2,a)'
      NumInt = l_NormModes
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,F1) ' ',Title
      Write(6,F1) ' ',
     &       '===================================================='
      Write(6,F1) ' ',
     &       ' mode        1/cm             GHz          hartrees '
      Write(6,F1) ' ',
     &       '----------------------------------------------------'
      Do i = 1,NumInt
      frequency = Freq(i)
      Write(6,'(A3,I2,A1,A3,F14.8,F18.8,F12.8)')
     &               ' ',NormModes(i),'.',' ',frequency*HarToRcm,
     &               frequency*HarToGHz,frequency
      End Do
      Write(6,F1) ' ',
     &       '===================================================='
      Write(6,*)
      Write(6,*)
C!
      End
