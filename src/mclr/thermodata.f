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
       Subroutine ThermoData(in_Freq, in_nFreq)
*
       Use Constants, only: auTocm
       use Temperatures, only: DefTemp
       Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
*----- Compute thermodynamic data at different temperatures.
       Real*8 in_Freq(in_nFreq), Freq(MxAtom*3-6)

*
*----- Remove translational and rotational frequencies
*
       nFreq = 0
       Do i = 1, in_nFreq
          If (in_Freq(i).gt.20.0d0) Then
             nFreq = nFreq + 1
             Freq(nFreq)=in_Freq(i)
          End if
       End Do
*
*----- Is the system linear?
*
       nAtom=(nFreq+6)/3
       nTR=3*nAtom-nFreq           ! Number of trans and rot fg
*
       Do i = 1, nFreq
*         Convert frequecnies from cm-1 to hartree
C         Freq(i) = Freq(i) * 4.55633538D-06
          Freq(i) = Freq(i) / auTocm
       End Do
*
       Do i = 1, Size(DefTemp)
          Call Thermo_Vib(nFreq,Freq,DefTemp(i),nTR,i)
       End Do
       End
*
       Subroutine Thermo_Vib(nFreq,Freq,T,nTR,iter)
       Use Constants, only: auTokcalmol, auTokJ, cal_to_j, kBoltzmann,
     &                      rNAVO
       Implicit Real*8 (a-h,o-z)
       Real*8 Freq(nFreq), T
       Integer first
*
       First=1
       Zero=0.0D00
       One=1.0D00
       Two=2.0D00
C      r_k = 1.38065800D-23
C      r_J2au=2.29371049D+17 ! Convert joules to atomic units
       r_J2au=1.0D-3 / auTokJ
*      Bolzmann's constant in a.u./ K
       rk = kBoltzmann * r_J2au
*      Write (*,*) 'Bolzmann''s constant=',rK
       if(T.eq.Zero) then
          beta = 1.0D99
        else
          beta = One / (rk*T)
        endif
************************************************************************
*                                                                      *
*      The Canonical partition function Q for indistinguishable        *
*      molecules                                                       *
*                                                                      *
*      is Q = q^N /N!, where N is the number of indistinguishable      *
*                                                                      *
*      molecules and q is the function with                            *
*                                                                      *
*      q = q_e * q_v * q_t * q_r                                       *
*                                                                      *
*      N = n * N_A (N_A is Avogrado's number) and the molecular        *
*                                                                      *
*      partition function is defined as q_m=q/n.                       *
*                                                                      *
*      We approximate q as follows                                     *
*                                                                      *
*      electronic:                                                     *
*      q_e=g_e (g_e is the degeneracy factor, for a singlet            *
*               ground state g_e=1, for a triplet ground state         *
*               g_e to a good approximation is 3, etc).                *
*                                                                      *
*      vibrational:                                                    *
*      q_v=q_v(1)*q_v(2)*...                                           *
*                                                                      *
*            q_v(i)= exp(-e(i)*beta/2) / (1-exp(-e(i)*beta))           *
*                                                                      *
*      translational:                                                  *
*      q_t=V/L^3, where L=(h^2 beta/2 pi m)^(1/2)                      *
*                                                                      *
*      d lnq_t/d beta = - 3kT/2                                        *
*                                                                      *
*      rotational: (high temperature approximations!)                  *
*      for linear systems                                              *
*      q_r=1/(beta sigma hc B) sigma=symmetry number                   *
*                                                                      *
*      d lnq_r/d beta = - 2kT/2                                        *
*                                                                      *
*      for nonlinear systems                                           *
*      q_r=(1/sigma)(1/(beta hc))^(3/2) (pi/ABC)^(1/2)                 *
*                                                                      *
*      d lnq_r/d beta = - 3kT/2                                        *
*                                                                      *
*      Some useful relations!                                          *
*                                                                      *
*      The internal energy:                                            *
*      U-U(0)=-(d lnQ/d beta)_V=-N(d lnq/d beta)_V                     *
*                                                                      *
*      Gibbs free energy:                                              *
*      G-G(0)=-nRT ln(q_m/N_A)                                         *
*                                                                      *
*      The enthalpy:                                                   *
*      H-H(0)=G-G(0)+T(S-S(0))                                         *
*                                                                      *
*            =-(d lnQ/d beta)_V + kTV(d lnQ/d V)_T                     *
*                                                                      *
*            =-N(d lnq/d beta)_V + NkTV(d lnQ/d V)_T                   *
*                                                                      *
*            =U-U(0)+pV                                                *
*                                                                      *
*      p=(NkT/q)(dq/dV)_T                                              *
*                                                                      *
*      since the translational part is the only which depends on the   *
*      volume we have that                                             *
*                                                                      *
*      H-H(0)=U-U(0)+NkT                                               *
*                                                                      *
************************************************************************
*
       Write (6,*)
       Write (6,*)
       Write (6,'(A,F6.2,A)') ' Temperature = ',T,' kelvin'
       Write (6,'(A)') ' ---------------------------'
       Write (6,*)
*
*------Iterate over vibrations, and compute the molecular
*      partition functions of the vibrations
*
       q_vib_Tot = One
       ZPVE  = Zero
       U_vib_Tot    = Zero
       Do i = 1, nFreq
          eta = Freq(i)
          If (iter.eq.first) Then
          Write (6,*) ' Vibrational temperature =',eta/rk
          End if
          If (eta.gt.Zero) Then
             ZPVE  = ZPVE  + eta/Two
             If (T.eq.Zero) Then
                q_vib = Zero
                U_vib = (eta/Two)
             Else
*               Eq. (6-20)
                q_vib = exp(-eta*beta/Two) / (One - exp(-eta*beta))
                U_vib = (eta/Two) + (eta / (exp(eta*beta)-One))
             End If
             q_vib_Tot = q_vib_Tot * q_vib
             U_vib_Tot = U_vib_Tot + U_vib
          End If
       End Do
*
*----- Add translational and rotational energy, kT*nTR/2,
*      the classiscal limit is used.
*
*      U-U(0)=-(d lnQ/d beta)_V, U(0)=0
*
*      lnQ = N lnq - lnN!
*
*      U-U(0)=-N(d lnq/d beta)_V
*
C      r_J2kcalmol=1.43932522D+20 ! conversion J to kcal/mol
       r_J2kcalmol= rNAVO / (1.0D3 * cal_to_J)
       U_TR=DBLE(nTR)*(kBoltzmann*T*r_J2kcalmol/Two)
*
*
*----- G-G(0)=-kT lnQ + kTV(d lnQ /dV)_T, G(0)=0
*
       If (T.eq.Zero) Then
          DeltaG = Zero
       Else
          DeltaG = -log(q_vib_Tot) * rk * T
       End If
*
C      auTokcalmol = 6.27509541D+2 ! Conversion au to kcal/mol
       DeltaG = DeltaG * auTokcalmol
*
       DeltaU = U_vib_Tot     * auTokcalmol
       ZPVE   = ZPVE   * auTokcalmol
       Write (6,'(A,F6.2,A)') '         DeltaG =',DeltaG,' kcal/mol'
       Write (6,'(A,F6.2,A)') '           ZPVE =',  ZPVE,' kcal/mol'
       Write (6,'(A,F6.2,A)') '      TotDeltaU =',DeltaU,' kcal/mol'
       Write (6,'(A,F6.2,A)') ' TotDeltaU-ZPVE =',DeltaU-ZPVE,
     &  ' kcal/mol'
*
       If (T.gt.Zero) Then
          DeltaS = 1.0d+3 * (DeltaU - DeltaG)/T
       Else
          DeltaS = Zero
       End If
       Write (6,'(A,F8.4,A)') '      Entropy S =',DeltaS,' cal/(mol*K)'
       Write (6,'(A,F8.4,A)') '         U(T&R) =',U_TR,' kcal/mol'
       Write (6,'(A,F8.4,A)') '       Tot(temp)=',U_TR+DeltaU-ZPVE,
     &                         ' kcal/mol'
*
       Return
       End
