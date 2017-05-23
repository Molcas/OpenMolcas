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
      Subroutine ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,
     &                 nUserPT,nsRot,iMult,nAtom,ipEVal,in_nFreq,
     &                 lSlapaf)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "constants.fh"
#include "WrkSpc.fh"
#include "real.fh"
      Integer in_nFreq, nFreq, iMult, nUserPT, nsRot, nAtom
      Real*8 UserT(64), UserP
      Real*8 TotalM, TRotA, TRotB, TRotC
      Real*8 Freq(MxAtom*3-6), VibT(MxAtom*3-6), dFreqI
      Real*8 Energy
      Logical lSlapaf
*
* --- If nUserPT.EQ.0 (no User-defined Pressure and Temperatures)
* --- Compute thermodynamic data at  1.00 atm &  298.15 Kelvin.
*
      If (nUserPT.EQ.0) then
        nUserPT=1
        UserP=1.0d0
        UserT(1)=298.15d0
      EndIf
      Call Get_dScalar('Last energy',Energy)
*
* --- Is the system linear?
*
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
*
* --- Remove translational and rotational frequencies
*
      nFreq = 0
      nTr2=nTr
      If(lSlapaf) nTr2=0
      Do i = 1, in_nFreq
        dFreqI = Work(ipEVal+i-1)
        If (dFreqI.GT.20.0d0) Then
          nFreq = nFreq + 1
          Freq(nFreq) = dFreqI
        End if
      End Do
      If ((in_nFreq-nFreq-nTR2).GT.0 ) Write (6,*) ' *** Warning: ',
     &   (in_nFreq-nFreq-nTR2),' vibrational contributions removed.'
*
* --- Convert frequencies from cm-1 to hartree
*
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
      Write (6,'(A,I2)')
     &            ' Number of trans. and rot. degrees of freedom: ',nTR
      Write (6,'(A,F9.3,A,F9.6,A)') ' ZPVE             ',
     &   ZPVE*6.27509541D+2,' kcal/mol     ',ZPVE,' au.'
      Write (6,'(A,13X,F15.6,A)') ' ZPVE corrected energy',
     &   ZPVE+Energy,' au.'
*
      Do i = 1, nUserPT
        Call Thermo_VibG(nFreq,Freq,UserT(i),UserP,TotalM,
     &                  nTR,nsRot,TRotA,TRotB,TRotC,iMult,Energy)
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nAtom)
      End
*
      Subroutine Thermo_VibG(nFreq,Freq,T,P,TotalM,
     &            nTR,nsRot,TRotA,TRotB,TRotC,iMult,Energy)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "constants.fh"
#include "constants2.fh"
      Real*8 Freq(nFreq), VibT(nFreq), T, P,Energy
*
      r_k  = 1.3806580d-23 ! Boltzmann / SI
      r_J2au= 1.0D-3 / CONV_AU_TO_KJ_
      rk = r_k * r_J2au ! Bolzmann constant in a.u./ K
      dNA  = CONST_AVOGADRO_ ! Avogadro
      dAU2kCal = 627.5095d0
      R_gas_kcal = 1.987216d-3
      q_e   = 1.0d0
      q_tr  = 1.0d0
      q_rot = 1.0d0
      q_vib_Tot = 1.0d0
      dS_e   = 0.0d0
      dS_tr  = 0.0d0
      dS_rot = 0.0d0
      dS_vib_Tot = 0.0d0
      dU_e   = 0.0d0
      dU_tr  = 0.0d0
      dU_rot = 0.0d0
      dU_vib_Tot = 0.0d0
      ZPVE       = 0.0d0
      q_TOT  = 1.0d0
      dS_TOT = 0.0d0
      dU_TOT = 0.0d0
      dH_TOT = 0.0d0
      dG_TOT = 0.0d0
*
* --- Electronic Contributions
*
*     Molecular Partition Function
      q_e    = iMult
*     Canonical Partition Function
      CQ_e = q_e
*     Entropy
      dS_e = R_gas_kcal * log(CQ_e)
*     Thermal
      dU_e = 0.0d0
*
* --- Translational Contributions
*
      If (T.GT.0.0d0) then
*       Molecular Partition Function (q/V in 1/m^3)
        dFact = 1.8793338d26        ! ((2*PI*k_B)/(h^2 *N_A*1000))^(3/2)
        dMT   = (TotalM/1.0d3) * T  ! m = PM/1000
        q_tr  = dFact * dMT
        q_tr  = q_tr * sqrt(dMT)
*       Entropy
        dVM   = 8.20575d-5 * T / P ! Gas Molar Volume / m**3
        dS_tr = 1.0d3 * R_gas_kcal *
     &  (log(dVM) + 1.5d0*(log(T)+log(TotalM/1.0d3)) + 18.605d0)
*       Thermal
        dU_tr = 1.5d0 * R_gas_kcal * T
      EndIf
*
* --- Rotational Contributions
*
      If (T.GT.0.0d0) then
*       Molecular Partition Function
        If (nTR.EQ.5) then
          q_rot = T / (nsRot*TRotC)
        Else If (nTR.EQ.3) Then
           q_rot = 1.0d0
        else
          q_rot = T*T*T
          q_rot = q_rot / TRotA
          q_rot = q_rot / TRotB
          q_rot = q_rot / TRotC
          q_rot = PI*q_rot
          q_rot = sqrt(q_rot)
          q_rot = q_rot/nsRot
        EndIf
      EndIf
*     Entropy
      If (nTR.EQ.5) then
        dS_rot = 1.0d3 * R_gas_kcal * ( log(q_rot) + 1.0d0 )
      Else If (nTR.EQ.3) Then
        dS_rot=0.0d0
      else
        dS_rot = 1.0d3 * R_gas_kcal * ( log(q_rot) + 1.5d0 )
      EndIf
*     Thermal
      If (nTR.EQ.5) then
        dU_rot = R_gas_kcal*T
      Else if (nTr.EQ.3) then
        dU_rot = 0.0d0
      else
        dU_rot = 1.5d0*R_gas_kcal*T
      EndIf
*
* --- Vibrational Contributions
*
      If (T.eq.Zero) then
         beta = 1.0D99
      else
         beta = One / (rk*T)
      EndIf
      Do i = 1, nFreq
        q_vib  = 1.0d0
        dU_vib = 0.0d0
        dS_vib = 0.0d0
        eta = Freq(i)
        VibT(i) = Freq(i)/rk
        If (eta.gt.Zero) Then
          ZPVE  = ZPVE  + eta/Two
          If (T.eq.Zero) Then
            q_vib  = One
            dU_vib = (eta/Two)
            dS_vib = Zero
          Else
*           Eq. (6-20)
            q_vib  = exp(-eta*beta/Two) / (One - exp(-eta*beta))
            dU_vib = (eta/Two) + (eta / (exp(eta*beta)-One))
            dS_vib = (eta*beta) / (exp(eta*beta) - One)
            dS_vib = dS_vib - log(One - exp(-eta*beta) )
          End If
          q_vib_Tot  =  q_vib_Tot *  q_vib
          dU_vib_Tot = dU_vib_Tot + dU_vib
          dS_vib_Tot = dS_vib_Tot + dS_vib
        End If
      End Do
      dU_vib_Tot = dU_vib_Tot * dAU2kCal
      dS_vib_Tot = dS_vib_Tot * R_gas_kcal * 1.0d3
*
      q_TOT  = q_e  * q_tr  * q_rot  * q_vib_Tot
      dS_TOT = dS_e + dS_tr + dS_rot + dS_vib_Tot
      dU_TOT = dU_e + dU_tr + dU_rot + dU_vib_Tot
      dH_TOT = dU_TOT + R_gas_kcal*T
      dG_TOT = dH_TOT - dS_TOT*T /1.0d3
*
* --- Print Results
*
      Write (6,*)
      Write (6,'(A)')
     &          ' *****************************************************'
      Write (6,'(A,F8.2,A,F7.2,A)') ' Temperature = ',T,
     &                              ' Kelvin, Pressure =',P,' atm'
      Write (6,'(A)')
     &          ' -----------------------------------------------------'
      Write (6,'(A)')
     &          ' Molecular Partition Function and Molar Entropy:'
      Write (6,'(A)')
     &          '                        q/V (M**-3)    S(kcal/mol*K)'
      Write (6,'(A,D17.6,F13.3)') ' Electronic       ',q_e,dS_e
      Write (6,'(A,D17.6,F13.3)') ' Translational    ',q_tr,dS_tr
      Write (6,'(A,D17.6,F13.3)') ' Rotational       ',q_rot,dS_rot
      Write (6,'(A,D17.6,F13.3)') ' Vibrational      ',q_vib_Tot,
     &                                                        dS_vib_Tot
      Write (6,'(A,D17.6,F13.3)') ' TOTAL            ',q_TOT,dS_TOT
*
      Write (6,*)
      Write (6,'(A)') ' Thermal contributions to INTERNAL ENERGY:'
      Write (6,'(A,F9.3,A,F9.6,A)') ' Electronic       ',
     &  dU_e      ,  ' kcal/mol     ',dU_e      /dAU2kCal,' au.'
      Write (6,'(A,F9.3,A,F9.6,A)') ' Translational    ',
     &  dU_tr     ,  ' kcal/mol     ',dU_tr     /dAU2kCal,' au.'
      Write (6,'(A,F9.3,A,F9.6,A)') ' Rotational       ',
     &  dU_rot    ,  ' kcal/mol     ',dU_rot    /dAU2kCal,' au.'
      Write (6,'(A,F9.3,A,F9.6,A)') ' Vibrational      ',
     &  dU_vib_Tot,  ' kcal/mol     ',dU_vib_Tot/dAU2kCal,' au.'
      Write (6,'(A,F9.3,A,F9.6,A)') ' TOTAL            ',
     &  dU_TOT    ,  ' kcal/mol     ',dU_TOT    /dAU2kCal,' au.'
      Write (6,*)
      Write (6,'(A)') ' Thermal contributions to'
      Write (6,'(A,F9.3,A,F9.6,A)') ' ENTHALPY         ',
     &  dH_TOT    ,  ' kcal/mol     ',dH_TOT    /dAU2kCal,' au.'
      Write (6,'(A,F9.3,A,F9.6,A)') ' GIBBS FREE ENERGY',
     &  dG_TOT    ,  ' kcal/mol     ',dG_TOT    /dAU2kCal,' au.'
      Write (6,*)
      Write (6,'(A)') ' Sum of energy and thermal contributions'
      Write (6,'(A,17X,F15.6,A)') ' INTERNAL ENERGY  ',
     &  Energy + dU_TOT    /dAU2kCal,' au.'
      Write (6,'(A,17X,F15.6,A)') ' ENTHALPY         ',
     &  Energy + dH_TOT    /dAU2kCal,' au.'
      Write (6,'(A,17X,F15.6,A)') ' GIBBS FREE ENERGY',
     &  Energy + dG_TOT    /dAU2kCal,' au.'
      Write (6,'(A)')
     &          ' -----------------------------------------------------'
*
      Return
      End
