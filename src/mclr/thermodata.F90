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
       Subroutine ThermoData(in_Freq, in_nFreq)
!
       Use Constants, only: auTocm
       use Temperatures, only: DefTemp
       Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
!----- Compute thermodynamic data at different temperatures.
       Real*8 in_Freq(in_nFreq), Freq(MxAtom*3-6)

!
!----- Remove translational and rotational frequencies
!
       nFreq = 0
       Do i = 1, in_nFreq
          If (in_Freq(i).gt.20.0d0) Then
             nFreq = nFreq + 1
             Freq(nFreq)=in_Freq(i)
          End if
       End Do
!
!----- Is the system linear?
!
       nAtom=(nFreq+6)/3
       nTR=3*nAtom-nFreq           ! Number of trans and rot fg
!
       Do i = 1, nFreq
!         Convert frequecnies from cm-1 to hartree
!         Freq(i) = Freq(i) * 4.55633538D-06
          Freq(i) = Freq(i) / auTocm
       End Do
!
       Do i = 1, Size(DefTemp)
          Call Thermo_Vib(nFreq,Freq,DefTemp(i),nTR,i)
       End Do
       End
