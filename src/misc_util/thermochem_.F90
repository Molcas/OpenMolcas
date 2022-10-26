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
      Subroutine ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,      &
     &                 nUserPT,nsRot,iMult,nAtom,EVal,in_nFreq,         &
     &                 lSlapaf)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "constants.fh"
#include "WrkSpc.fh"
#include "real.fh"
      Integer in_nFreq, nFreq, iMult, nUserPT, nsRot, nAtom
      Real*8 UserT(64), UserP, EVal(*)
      Real*8 TotalM, TRotA, TRotB, TRotC
      Real*8 Freq(MxAtom*3-6), VibT(MxAtom*3-6), dFreqI
      Real*8 Energy
      Logical lSlapaf
!
! --- If nUserPT.EQ.0 (no User-defined Pressure and Temperatures)
! --- Compute thermodynamic data at  1.00 atm &  298.15 Kelvin.
!
      If (nUserPT.EQ.0) then
        nUserPT=1
        UserP=1.0d0
        UserT(1)=298.15d0
      EndIf
      Call Get_dScalar('Last energy',Energy)
!
! --- Is the system linear?
!
      nTR = 6
      If (TRotA.GT.1.0d99) then
        TRotA = 0.0d0
        nTR = nTR -1
      EndIf
      If (TRotB.GT.1.0d99) then
        TRotB = 0.0d0
        nTR = nTR -1
      EndIf
      If (TRotC.GT.1.0d99) then
        TRotC = 0.0d0
        nTR = nTR -1
      EndIf
!
! --- Remove translational and rotational frequencies
!
      nFreq = 0
      nTr2=nTr
      If(lSlapaf) nTr2=0
      Do i = 1, in_nFreq
        dFreqI = EVal(i)
        If (dFreqI.GT.20.0d0) Then
          nFreq = nFreq + 1
          Freq(nFreq) = dFreqI
        End if
      End Do
      If ((in_nFreq-nFreq-nTR2).GT.0 ) Write (6,*) ' *** Warning: ',    &
     &   (in_nFreq-nFreq-nTR2),' vibrational contributions removed.'
!
! --- Convert frequencies from cm-1 to hartree
!
      r_k = CONST_BOLTZMANN_
      r_J2au=1.0D-3 / CONV_AU_TO_KJ_
      rk = r_k * r_J2au ! Bolzmann constant in a.u./ K
      ZPVE  = 0.0d0
      Do i = 1, nFreq
         Freq(i) = Freq(i) / CONV_AU_TO_CM1_ ! = 219474.625
         VibT(i) = Freq(i) / rk
         ZPVE = ZPVE  + Freq(i)/2.0d0
      End Do
      Write (6,'(A)') ' Vibrational temperature (K): '
      Do i = 1, nFreq, 5
        If (i.LE.(nFreq-4)) Write (6,'(1X,5F9.2)') (VibT(j),j=i,i+4)
        If (i.EQ.(nFreq-3)) Write (6,'(1X,4F9.2)') (VibT(j),j=i,i+3)
        If (i.EQ.(nFreq-2)) Write (6,'(1X,3F9.2)') (VibT(j),j=i,i+2)
        If (i.EQ.(nFreq-1)) Write (6,'(1X,2F9.2)') (VibT(j),j=i,i+1)
        If (i.EQ.(nFreq))   Write (6,'(1X, F9.2)') VibT(i)
      EndDo
      Write (6,'(A,I2)')                                                &
     &            ' Number of trans. and rot. degrees of freedom: ',nTR
      Write (6,'(A,F9.3,A,F9.6,A)') ' ZPVE             ',               &
     &   ZPVE*6.27509541D+2,' kcal/mol     ',ZPVE,' au.'
      Write (6,'(A,13X,F15.6,A)') ' ZPVE corrected energy',             &
     &   ZPVE+Energy,' au.'
!
      Do i = 1, nUserPT
        Call Thermo_VibG(nFreq,Freq,UserT(i),UserP,TotalM,              &
     &                  nTR,nsRot,TRotA,TRotB,TRotC,iMult,Energy)
      End Do
!
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(nAtom)
      End
