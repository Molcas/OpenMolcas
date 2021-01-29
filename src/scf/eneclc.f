************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2003, Valera Veryazov                                  *
************************************************************************
      SubRoutine EneClc(En1V,En2V,EnerV,Dens,OneHam,TwoHam,mBT,mDens,nD,
     &                  EDFT,nEDFT)
************************************************************************
*                                                                      *
* Purpose: Compute one- and two-electron energies                      *
*                                                                      *
* output:                                                              *
*   En1V    : one-electron energy (variational)                        *
*   En2V    : two-electron energy (variational)                        *
*   EnerV   : En1V + En2V                                              *
*                                                                      *
* called from: SCF_Energy                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* written by:                                                          *
* P.O. Widmark, M.P. Fuelscher and P. Borowski                         *
* University of Lund, Sweden, 1992                                     *
*                                                                      *
************************************************************************
#ifdef _FDE_
      Use SCF_Arrays, Only: Emb
#endif
      use OFembed, only: Do_OFemb
      Implicit Real*8 (a-h,o-z)
*
* Declaration of procedure parameters
*
      REAL*8 En1V,En2V,EnerV
      Real*8 Dens(mBT,nD,mDens), OneHam(mBT), TwoHam(mBT,nD,mDens),
     &       EDFT(nEDFT)
#include "real.fh"

#include "mxdm.fh"
#include "infscf.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
      COMMON  / OFembed_R / Rep_EN,Func_AB,Func_A,Func_B,Energy_NAD,
     &                      V_Nuc_AB,V_Nuc_BA,V_emb
*----------------------------------------------------------------------*
* Start                                                                *
*----------------------------------------------------------------------*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*define _DEBUGPRINT_
*
* Allocate memory for full Dens and TwoHam
*
c set to Zero for RHF
      En1V_ab=0.0D0
      En2V_ab=0.0D0

      iter_d=iter-iter0
*
      En1V  = DDot_(nBT,OneHam,1,Dens(1,1,iPsLst),1)
      If(iUHF.eq.1) Then
         En1V_ab  = DDot_(nBT,OneHam,1,Dens(1,2,iPsLst),1)
      End If
*
      E_DFT = EDFT(iter_d)
*
#ifdef _FDE_
      ! Embedding
      if (embPot) Eemb = DDot_(nBT*nD,Emb,1,Dens(1,1,iPsLst),1)
#endif
*
*     If just one electron make sure that the two-electron energy
*     is zero.
*
      nElec=0
      Do iSym = 1, nSym
         nElec = nElec
     &         + (2-iUHF)*nOcc(iSym,1)
     &         +    iUHF *nOcc(iSym,2)
      End Do
      En2V=0
      If ((nElec.le.1).and.(KSDFT.eq.'SCF')) Go To 999
*
      En2V  = DDot_(nBT,TwoHam(1,1,iPsLst),1,Dens(1,1,iPsLst),1)
      If (iUHF.eq.1) Then
         En2V_ab  = DDot_(nBT,TwoHam(1,2,iPsLst),1,Dens(1,2,iPsLst),1)
      End If
999   Continue
*
      If (Do_OFemb) Then
         If(iUHF.eq.1) Then ! equipartition
            En2V   = En2V    -Half*Rep_EN
            En2V_ab= En2V_ab -Half*Rep_EN
         Else
            En2V   = En2V -        Rep_EN
         End If
      End If
*
*     Note that the DFT energy can not be computed as a trace.
*
      If(iUHF.eq.1) Then
         Elst(iter,1)=En1V   +Half*En2V   +Half*PotNuc+Half*E_DFT
         Elst(iter,2)=En1V_ab+Half*En2V_ab+Half*PotNuc+Half*E_DFT
      Else
         Elst(iter,1)=En1V   +Half*En2V        +PotNuc+     E_DFT
      End If
*
      If(iUHF.eq.0) Then
         En2V  = Half*En2V
      Else
         En2V  = Half*(En2V+En2V_ab)
      End If
      En1V= (En1V+En1V_ab) + E_DFT
      EnerV = En1V + En2V + PotNuc
#ifdef __SUNPRO_F90
      If (iUHF.gt.3) Write (6,*) 'eneclc: Ene=',En1V,En1V_ab,En2V,EnerV
#endif
#ifdef _DEBUGPRINT_
      Write (6,*) 'eneclc: Ene=',En1V,En1V_ab,En2V,EnerV
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(14) = TimFld(14) + (Cpu2 - Cpu1)
*----------------------------------------------------------------------*
* Exit                                                                 *
*----------------------------------------------------------------------*
      Return
      End
