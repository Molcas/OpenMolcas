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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetInc()
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Initialize all entries in LDF include files.
C
C     localdf.fh
C     localdf_print.fh
C     localdf_bas.fh
C     localdf_int.fh
C     localdf_int2.fh
C     localdf_int3.fh
C     ldf_atom_info.fh
C     ldf_atom_pair_info.fh
C     ldf_cio.fh
C     ldf_a2ap.fh
C     ldf_integral_prescreening_info.fh
C     ldf_qdiag.fh
C     ldf_oneel.fh
C     ldf_charge_constraint_info.fh
C     ldf_atomiclabels.fh
C
      Implicit None
#include "localdf.fh"
#include "localdf_print.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "ldf_atom_info.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_cio.fh"
#include "ldf_a2ap.fh"
#include "ldf_integral_prescreening_info.fh"
#include "ldf_qdiag.fh"
#include "ldf_oneel.fh"
#include "ldf_charge_constraint_info.fh"
#include "ldf_atomiclabels.fh"

C     localdf.fh
C     ===========

      LDF2=.False.
      VerifyFit=.False.
      CheckPairIntegrals=.False.
      CheckOverlapIntegrals=.False.
      Thr_Prescreen=-1.0d9
      Thr_Accuracy=-1.0d9
      LDF_Run_Mode=0
      LDF_Constraint=-1
      WriteUnconstrainedC=.False.
      UseUniqueAtomPairs=.False.

C     localdf_print.fh
C     =================

      iPrint=0

C     localdf_bas.fh
C     ===============

      nBas_Valence=0
      nBas_Auxiliary=0
      nShell_Valence=0
      nShell_Auxiliary=0
      ip_iSOShl=0
      l_iSOShl=0
      ip_iShlSO=0
      l_iShlSO=0
      ip_nBasSh=0
      l_nBasSh=0

C     localdf_int.fh
C     ===============

      SHA=0
      SHB=0
      SHC=0
      SHD=0
      SPAB=0
      SPCD=0
      ip_IndxG=0
      l_IndxG_1=0
      l_IndxG_2=0
      ip_IndxG2=0
      l_IndxG2_1=0
      l_IndxG2_2=0
      ip_2CList=0
      l_2CList_1=0
      l_2CList_2=0
      ip_iOff=0
      l_iOff=0
      nRow_G=0
      nRow_uvJ=0
      iOffuv=0

C     localdf_int2.fh
C     ================

      Call LDF_Set_localdf_int2()

C     localdf_int3.fh
C     ================

      Call LDF_Set_localdf_int3()

C     ldf_atom_info.fh
C     =================

      LDF_AtomInfo_Status=LDF_AtomInfo_Unset
      NumberOfAtoms=0
      ip_Coord=0
      l_Coord=0
      ip_A_Unique=0
      l_A_Unique=0
      ip_A_Shells=0
      l_A_Shells=0
      ip_A_AuxShells=0
      l_A_AuxShells=0

C     ldf_atom_pair_info.fh
C     ======================

      LDF_AtomPairInfo_Status=LDF_AtomPairInfo_Unset
      NumberOfAtomPairs=0
      ip_AP_Atoms=0
      l_AP_Atoms=0
      ip_AP_Unique=0
      l_AP_Unique=0
      ip_AP_Diag=0
      l_AP_Diag=0
      ip_AP_DiagBak=0
      l_AP_DiagBak=0
      ip_AP_1CLinDep=0
      l_AP_1CLinDep=0
      ip_AP_2CFunctions=0
      l_AP_2CFunctions=0
      ip_AP_DiskC=0
      l_AP_DiskC=0

C     ldf_cio.fh.
C     ============

      Lu_LDFC=0
      LastAtomPair=0
      ip_LDFC_Buffer=0
      l_LDFC_Buffer=0
      ip_LDFC_Blocks=0
      l_LDFC_Blocks=0

C     ldf_a2ap.fh.
C     =============

      ip_A2AP=0
      l_A2AP=0

C     ldf_integral_prescreening_info.fh.
C     ===================================

      ip_GDiag_1C=0
      l_GDiag_1C=0
      ip_GDiag_1C_Mx=0
      l_GDiag_1C_Mx=0
      ip_GDiag_1C_Sm=0
      l_GDiag_1C_Sm=0
      ip_GDiag_2C=0
      l_GDiag_2C=0
      ip_GDiag_2C_Mx=0
      l_GDiag_2C_Mx=0
      ip_GDiag_2C_Sm=0
      l_GDiag_2C_Sm=0
      ip_IDiag=0
      l_IDiag=0
      ip_IDiag_Mx=0
      l_IDiag_Mx=0
      ip_IDiag_Sm=0
      l_IDiag_Sm=0

C     ldf_qdiag.fh.
C     ==============

      ip_AP_QDiag=0
      l_AP_QDiag=0

C     ldf_oneel.fh.
C     ==============

      OperatorLabel='IS_UNSET'
      nComp=0
      nIC=0
      Call iZero(iStabO,8)
      nStabO=0
      ip_lOper=0
      l_lOper=0
      ip_kOper=0
      l_kOper=0
      ip_CCoor=0
      l_CCoor=0
      ip_xZeta=0
      l_xZeta=0
      ip_xZI=0
      l_xZI=0
      ip_xKappa=0
      l_xKappa=0
      ip_xPCoor=0
      l_xPCoor=0
      rHrmt=-9.9d9

C     ldf_charge_constraint_info.fh.
C     ===============================

      ChargeConstraintSet=.False.
      ip_CC_AuxIntVec=0
      ip_CC_Overlap=0
      l_CC_Overlap=0
      ip_CC_lambda=0
      l_CC_lambda=0

C     ldf_atomiclabels.fh.
C     =====================

      AtomicLabelsSet=.False.
      ip_AtomicLabels=0
      l_AtomicLabels=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Set_localdf_int2()
      Implicit None
#include "localdf_int2.fh"

      SHA=0
      SHB=0
      SHC=0
      SHD=0

      SPAB=0
      SPCD=0

      ip_AB_IndxG=0
      l_AB_IndxG_1=0
      l_AB_IndxG_2=0
      ip_AB_IndxG2=0
      l_AB_IndxG2_1=0
      l_AB_IndxG2_2=0
      ip_AB_2CList=0
      l_AB_2CList_1=0
      l_AB_2CList_2=0

      ip_CD_IndxG=0
      l_CD_IndxG_1=0
      l_CD_IndxG_2=0
      ip_CD_IndxG2=0
      l_CD_IndxG2_1=0
      l_CD_IndxG2_2=0
      ip_CD_2CList=0
      l_CD_2CList_1=0
      l_CD_2CList_2=0

      nAB=0
      nCD=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Set_localdf_int3()
      Implicit None
#include "localdf_int3.fh"

      SHA=0
      SHB=0
      SHC=0
      SHD=0

      nRow=0
      iRow0=0
      iCol0=0

      End
