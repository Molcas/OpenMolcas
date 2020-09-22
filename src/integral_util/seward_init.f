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
* Copyright (C) 1990,2020, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Subroutine Seward_Init
************************************************************************
*                                                                      *
*     Object: to set data which is stored in common blocks             *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose      *
*             January 1990                                             *
************************************************************************
      use EFP_Module
      use k2_arrays
      use Basis_Info
      use Real_Info, only: ThrInt
      implicit real*8 (a-h,o-z)
      External Reduce_Prt
      Logical Reduce_Prt
#include "itmax.fh"
#include "info.fh"
#include "pstat.fh"
#include "print.fh"
#include "notab.fh"
#include "status.fh"
#include "twoswi.fh"
#include "rmat.fh"
#include "gam.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "relae.fh"
#include "FMM.fh"
#include "nac.fh"
#include "srint.fh"
      Logical lGENINT
      Character(LEN=180) Env
*                                                                      *
************************************************************************
*                                                                      *
C
C...  Parameters for srint
C
      shortrange=.False.
      isr_simulate=0
*
*     Initialize FMM.fh
*
      FMM_shortrange = .False.
      asymptotic_Rys = .False.
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Info
*
      Seward_Status=InActive
*
*     Info
*
      nMltpl=2
      iAngMx=-1
      nWel=0
      iRI_type=0
      Max_Center=15

      KVector(:)=Zero

      nRP=0
      iPAMcount=1

      IsChi=0
*
*-----LInfo
*
      NEMO=.False.
      Do_RI=.False.
      IRFLAG1=0
      DirInt=.False.
      EMFR  =.False.
      UnNorm=.False.
      lSchw=.True.
      Vlct=.True.
      lUPONLY=.False.
      lDOWNONLY=.False.
      lRel=.False.
      lAMFI=.False.
      lGENINT=.False.
      GIAO=.False.
      Cholesky=.False.
      Do_FckInt=.True.
      Do_GuessOrb=.True.
*
      Do_acCD_Basis=.True.
      Do_nacCD_Basis=.False.
      Skip_High_AC = .False.
      LDF=.False.
      LocalDF=.False.
*
      lRP=.False.
      Align_Only=.False.
      Do_Align=.True.
      VarT=.False.
      VarR=.False.
      FNMC=.False.
*
      Call GetEnvF('MOLCAS_NEW_DEFAULTS', Env)
      Call UpCase(Env)
      If (Env.eq.'YES') Then
         Do_RI=.True.
         iRI_Type=4
      End If
*
      Shake=-One
*
*     Flags to control build of FMM short-range integral components
*
      DoFMM = .False.
*
*-----RInfo
*
      RadMax=0.00d+00
*
      Thrshld_CD=1.0D-4
      Delta_RICD=0.0D0
      E1=0.0D0
      E2=0.0D0
      SadStep=0.1d0
      ChiI2=0.0D0
*
*     Flags to control build of FMM short-range integral components
*
      RPQMIN = 0.4d0
*
*-----PStat
*
      r1=0.d0
      r2=0.d0
      r3=0.d0
      r4=0.d0
      q1=0.d0
      q2=0.d0
      q3=0.d0
      q4=0.d0
      MaxReq=0
      MinXtr=0
      iTotal=0
      MaxMem=0
*
*-----Print
*
      iPL=iPrintLevel(-1)
      If (iPL.eq.2) Then
         iPL=5
      Else If (iPL.eq.3) Then
         iPL=6
      Else If (iPL.eq.4) Then
         iPL=7
      Else If (iPL.eq.5) Then
         iPL=49  ! 99 would be just too much
      End If
      Do i = 1, nRout
         nPrint(i)=iPL
      End Do
      If ((Reduce_Prt().and.iPL.lt.6).or.iPL.eq.0) Then
         Show=.False.
      Else
         Show=.True.
      End If
*
*-----NoTab
*
      NoTab=.false.
*
      NDDO=.False.
*
      GT_Status=InActive
      T_Status =InActive
      PP_Status=InActive
      k2_Status=InActive
      RctFld_Status=InActive
      ERI_Status=InActive
      Indexation_Status=Inactive
      XMem_Status=Inactive
      NQ_Status=Inactive
      Seward_Status=Active
*
      Call Set_Binom
      Call Set_CanInd
*
*     Set some default value for RMAT type integration
*
*---- rmat.fh
*
      RmatR     = 10.0D0
      Epsabs    = 10.D-10
      Epsrel    = 1.D-14
      qCoul     = 0.0D0
      Epsq      = 1.D-8
      bParm     = 0.0D0
      dipol(1)  = 0.0D0
      dipol(2)  = 0.0D0
      dipol(3)  = 0.0D0
      Dipol1    = 0.0D0
      keyr      = 6
      Quadpack  = .True.
      nagint    = .False.
      testint   = .False.
      RMat_On   = .False.
*
*---- gam.fh
*
      lgamma = 9
*
*     relae.fh
*
      irelae=-1
*
      Call DCR_Init()
*
      Call Set_Basis_Mode('Valence')
*
      Call Mk_TriInd()
*
      Call CovRadT_Init()
*
*     nac.fh
*
      isNAC = .False.
      isCSF = .False.
*
*     EFP stuff
*
      lEFP=.False.
      nEFP_fragments=0
*
      Return
      End
