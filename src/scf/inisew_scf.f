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
* Copyright (C) 1992, Roland Lindh                                     *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine IniSew_scf(DSCF,EThr,SIntTh,KSDFT)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             modified by M.Schuetz @teokem.lu.se, 1995                *
* input:      EThr,DThr,FThr,DltNTh: div. Threshold values for         *
*               of SCF WF. (only relevant for direct SCF, DSCF=TRUE)   *
* output:     SIntTh: computed cutoff for integrals (prescreening)     *
*                                                                      *
*                                                                      *
* Note :  the corresponding finalization subroutine is ClsSew          *
************************************************************************
      use Sizes_of_Seward, only: S
      use Gateway_Info, only: ThrInt, CutInt
      use OFembed, only: Do_OFemb
      use RICD_Info, only: Do_DCCD
      use InfSCF, only: nDisc, nCore
      Implicit Real*8 (A-H,O-Z)
      External EFP_On
#include "print.fh"
#include "addcorr.fh"
*
*
      Logical DSCF, RF_On, Langevin_On, PCM_On, EFP_On
      Character*(*) KSDFT
      Real*8 EThr
*
      If (DSCF.or.RF_On().or.Langevin_On().or.KSDFT.ne.'SCF'
     &                                    .or.Do_Addc
     &                                    .or.Do_Tw
     &                                    .or.Do_OFemb
     &                                    .or.EFP_On()) Then
         nDiff=0
         If (Langevin_On().and.S%iAngMx.eq.0) nDiff=1
         Call IniSew(DSCF.or.Langevin_On().or.PCM_On(),nDiff)
      End If
*
      If (Do_DCCD) Then
         nCore=0
         nDisc=0
      End If
      If (DSCF) Then
         CutInt=EThr*Min(1.0D-7,1.0D0/DBLE(S%nDim)**2)
         ThrInt=Cutint
         SIntTh=CutInt
      End If
      Return
      End
      Function Get_ThrInt()
      use Gateway_Info, only: ThrInt
      Implicit Real*8 (A-H,O-Z)
      Real*8 Get_ThrInt
      Get_ThrInt=ThrInt
      Return
      End
      Subroutine xSet_ThrInt(tmp)
      use Gateway_Info, only: ThrInt
      Implicit Real*8 (A-H,O-Z)
      ThrInt=tmp
      Return
      End
