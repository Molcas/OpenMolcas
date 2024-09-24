!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine EneClc(En1V,En2V,EnerV)
!***********************************************************************
!                                                                      *
! Purpose: Compute one- and two-electron energies                      *
!                                                                      *
! output:                                                              *
!   En1V    : one-electron energy (variational)                        *
!   En2V    : two-electron energy (variational)                        *
!   EnerV   : En1V + En2V                                              *
!                                                                      *
!***********************************************************************
#ifdef _FDE_
      use Embedding_Global, only: Eemb, embInt, embPot
#endif
      use OFembed, only: Do_OFemb
      use OFembed, only: Rep_EN
      use Constants, only: Zero, Half
      use InfSCF, only: ipsLst, nD, Iter, nSym, KSDFT, PotNuc, ELst, nBT, nOcc, TimFld
      use SCF_Arrays, only: OneHam, TwoHam, Dens, EDFT
      Implicit None
!
! Declaration of procedure parameters
!
      Real*8 En1V,En2V,EnerV

      Real*8 :: En1V_AB, En2V_AB, E_DFT, CPU1, CPU2, Tim1, Tim2, Tim3
      Integer nElec, iSym
      Real*8, External :: DDot_

!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
!
! Allocate memory for full Dens and TwoHam
!
! set to Zero for RHF
      En1V_ab=Zero
      En2V_ab=Zero
!
      En1V  = DDot_(nBT,OneHam,1,Dens(1,1,iPsLst),1)
      If(nD==2) Then
         En1V_ab  = DDot_(nBT,OneHam,1,Dens(1,2,iPsLst),1)
      End If
!
      E_DFT = EDFT(iter)
!
#ifdef _FDE_
      ! Embedding
      if (embPot) Eemb = DDot_(nBT*nD,embInt,1,Dens(1,1,iPsLst),1)
#endif
!
!     If just one electron make sure that the two-electron energy
!     is zero.
!
      nElec=0
      If (nD==1) Then
         Do iSym = 1, nSym
            nElec = nElec + 2*nOcc(iSym,1)
         End Do
      Else
         Do iSym = 1, nSym
            nElec = nElec + nOcc(iSym,1) + nOcc(iSym,2)
         End Do
      End If
      En2V=0
      If ((nElec.le.1).and.(KSDFT.eq.'SCF')) Go To 999
!
      En2V  = DDot_(nBT,TwoHam(1,1,iPsLst),1,Dens(1,1,iPsLst),1)
      If (nD==2) Then
         En2V_ab  = DDot_(nBT,TwoHam(1,2,iPsLst),1,Dens(1,2,iPsLst),1)
      End If
999   Continue
!
      If (Do_OFemb) Then
         If(nD==2) Then ! equipartition
            En2V   = En2V    -Half*Rep_EN
            En2V_ab= En2V_ab -Half*Rep_EN
         Else
            En2V   = En2V -        Rep_EN
         End If
      End If
!
!     Note that the DFT energy cannot be computed as a trace.
!
#ifdef _DEBUGPRINT_
      If(nD==2) Then
         Write(6,*) 'EnerClc:',En1V,En2V,PotNuc,E_DFT
         Write(6,*) 'EnerClc:',En1V_ab,En2V_ab
      Else
         Write(6,*) 'EnerClc:',En1V,En2V,PotNuc,E_DFT
      End If
#endif
      If(nD==2) Then
         Elst(iter,1)=En1V   +Half*En2V   +Half*PotNuc+Half*E_DFT
         Elst(iter,2)=En1V_ab+Half*En2V_ab+Half*PotNuc+Half*E_DFT
      Else
         Elst(iter,1)=En1V   +Half*En2V        +PotNuc+     E_DFT
      End If
!
      If(nD==1) Then
         En2V  = Half*En2V
      Else
         En2V  = Half*(En2V+En2V_ab)
      End If
      En1V= (En1V+En1V_ab) + E_DFT
      EnerV = En1V + En2V + PotNuc
#ifdef _DEBUGPRINT_
      Write (6,*) 'EneClc: Ene=',En1V,En1V_ab,En2V,EnerV
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(14) = TimFld(14) + (Cpu2 - Cpu1)
!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*
      Return
      End SubRoutine EneClc
