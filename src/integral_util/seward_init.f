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
      Character*180 Env
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
c     data iPhase/ 1, 1, 1,   -1, 1, 1,   1,-1, 1,  -1,-1, 1,
c    &             1, 1,-1,   -1, 1,-1,   1,-1,-1,  -1,-1,-1/
      Seward_Status=InActive
      do 10 j=0,4
      do 11 i=1,3
      iPhase(i,j)=1
11    continue
10    continue
      iPhase(1,1)=-1
      iPhase(2,2)=-1
      iPhase(1,3)=-1
      iPhase(2,3)=-1
      iPhase(3,4)=-1
      do 20 j=5,7
      do 21 i=1,3
      iPhase(i,j)=-1
21    continue
20    continue
      iPhase(2,5)=1
      iPhase(1,6)=1
*
*     Info
*
      MemHid=1
      nMltpl=2
      m2Max=0
      iAngMx=-1
      nWel=0
      iRI_type=0
      jMax = 5
      nTtl=0
      Max_Center=15

      iChCar(:)=0
      KVector(:)=Zero
      lOffAO(1:Mxdbsc)=0
      kOffAO(1:Mxdbsc,0:MxAng)=0
      iAOtSO(1:MxAO,0:7)=-999999999
      MaxBas(0:MxAng)=0
      MaxPrm(0:MxAng)=0
      IrrCmp(1:MxUnq)=0
      NrInt(-20:9)=0
      iSkip(0:7)=0

      nOrdEF=-1
      nDMS=0
      nRP=0
      iPack=0
      iSquar=0
      iWRopt=0
      iPAMcount=1
      Do i=1, Mxdbsc
         ECP(i)      =.false.
         AuxCnttp(i) =.false.
         FragCnttp(i)=.false.
         FockOp(i)   = .False.
         IsMM(i)         =0
         Parent_iCnttp(i)=0
         mdciCnttp(i)    =0
         ExpNuc(i)       =-One
         w_mGauss(i)     = One
         aCD_Thr(i)      = One
         fmass(i)        = One
      End Do
      inttot=0
      nOrd_XF = 1
      iOrdFm=0
      iXPolType=0
      IsChi=0
      MolWgh=2
*
      Do i=1,MxShll
         AuxShell(i) =.False.
      End Do
*
*-----LInfo
*
      NEMO=.False.
      Do_RI=.False.
      Primitive_Pass=.True.
      DKroll=.False.
      LDKroll=.False.
      IRFLAG1=0
      BSS   =.False.
      Onenly=.False.
      DirInt=.False.
      Expert=.False.
      EMFR  =.False.
      Petite=.False.
      lSOInt=.True.
      UnNorm=.False.
      lSchw=.True.
      Test=.False.
      Vlct=.True.
      lOAM=.False.
      lUPONLY=.False.
      lDOWNONLY=.False.
      lOMQ=.False.
      lDMS=.False.
      lRel=.False.
      SW_FileOrb='INPORB'
      Prprt=.False.
      Short=.True.
*--sdong, Apr. 2018--*
      ifallorb=.False.
*--sdong end---------*
      lECP=.False.
      lAux=.False.
      lPAM2=.False.
      Dist=.False.
      lXF=.False.
      lPP=.False.
      lAMP=.False.
      lAMFI=.False.
      lGENINT=.False.
      Nuclear_Model=Point_Charge
      force_part_c=.False.
      force_part_p=.False.
      GIAO=.False.
      Cholesky=.False.
      lFAIEMP=.False.
      Do_FckInt=.True.
      Do_GuessOrb=.True.
      Do_acCD_Basis=.True.
      Do_nacCD_Basis=.False.
      Skip_High_AC = .False.
      LDF=.False.
      LocalDF=.False.
      lRP=.False.
      Align_Only=.False.
      Do_Align=.True.
      Align_Weights='MASS'
      Do_Numerical_Gradients=.False.
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
      Sum=0.00d+00
      Sumsq=0.00d+00
      SumAbs=0.00d+00
      RadMax=0.00d+00
      AccMch=1.d-15
*
*     Integral thresholds
*
      ThrInt=1.d-14
      CutInt=1.d-16
*
*     Two-electron integral packing threshold
*
      PkAcc=ThrInt
*
      Rtrnc = Three
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
      Thrs=1.d-6
*
*-----CInfo
*
      Bline=' '
      Do i = 1, 10
         Title(i)=' '
      End Do
      Do i = 1, Mxdbsc
         Bsl    (i) = ' '
         Bsl_Old(i) = ' '
      End Do
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
      Call iPrmt_Init()
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
