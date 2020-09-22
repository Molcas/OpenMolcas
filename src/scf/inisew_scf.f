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
      SubRoutine IniSew_scf(DSCF,EThr,DThr,FThr,
     &                  DltNTh,SIntTh,KSDFT)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             modified by M.Schuetz @teokem.lu.se, 1995                *
* input:      EThr,DThr,FThr,DltNTh: div. Threshold values for         *
*               of SCF WF. (only relevant for direct SCF, DSCF=TRUE)   *
*             iPrLV(2*MxPrLv): int vector with routes and corres-      *
*               ponding print levels for SEWARD integral routines      *
* output:     SIntTh: computed cutoff for integrals (prescreening)     *
*                                                                      *
*                                                                      *
* Note :  the corresponding finalization subroutine is ClsSew          *
*                                                                      *
************************************************************************
      use Sizes_of_Seward, only: S
      use Real_Info, only: ThrInt
      Implicit Real*8 (A-H,O-Z)
      External EFP_On
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "iprlv.fh"
      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      Logical Do_Tw
      COMMON  / Tw_corr_L   / Do_Tw
      Logical Do_Addc
      COMMON  / ADDcorr_L   / Do_Addc
*
*
      Logical DSCF, RF_On, Langevin_On, PCM_On, EFP_On
      Character*(*) KSDFT
      Real*8 EThr,DThr,FThr,DltNTh
*
      If (DSCF.or.RF_On().or.Langevin_On().or.KSDFT.ne.'SCF'
     &                                    .or.Do_Addc
     &                                    .or.Do_Tw
     &                                    .or.Do_OFemb
     &                                    .or.EFP_On()) Then
         nDiff=0
         If (Langevin_On().and.iAngMx.eq.0) nDiff=1
         Call IniSew(DSCF.or.Langevin_On().or.PCM_On(),nDiff)
      End If
*
      If (DSCF) Then
C        CutInt=Min(EThr,DThr,FThr,DltNTh)*1.0d-5
         CutInt=EThr*Min(1.0D-7,1.0D0/DBLE(S%nDim)**2)
         Thrint=Cutint
         SIntTh=CutInt
      End If
      i=1
  100 If (iPrLV(i).le.0) Go To 110
       If (iPrLV(i).le.nRout) nPrint(iPrLV(i))=iPrLV(i+1)
       i=i+2
       Go To 100
  110 Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(DThr)
         Call Unused_real(FThr)
         Call Unused_real(DltNTh)
      End If
      End
      Function Get_ThrInt()
      use Real_Info, only: ThrInt
      Implicit Real*8 (A-H,O-Z)
      Real*8 Get_ThrInt
      Get_ThrInt=ThrInt
      Return
      End
      Subroutine xSet_ThrInt(tmp)
      use Real_Info, only: ThrInt
      Implicit Real*8 (A-H,O-Z)
      ThrInt=tmp
      Return
      End
