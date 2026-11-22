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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

! The interface of the SGUGA library: re-exports the public entities of
! its modules (SGUGA_Core, SGUGA_Tables, SGUGA_Dets, SGUGA_Epq,
! SGUGA_Setup, SGUGA_MatEl, SGUGA_Triplet, SGUGA_Cnf, SGUGA_Checks).
!#define _TEST_CSFDET_
module SGUGA

use SGUGA_Core
use SGUGA_Tables
use SGUGA_Dets
use SGUGA_Epq
use SGUGA_Setup
use SGUGA_MatEl
use SGUGA_Triplet
use SGUGA_Cnf
use SGUGA_Checks

implicit none
private

public :: SGStruct, CIStruct, EXStruct, SG_Free, SG_Init, SG_Init_Light, SG_Epq_Psi, MkCOT, MkSgNum, SG_NUM, SG_PHASE
public :: SG_EpqEqp_Diag, sg_csfdet, generate_detspin, generate_twos, StepVec, SG_TRACI
public :: SG_HPsi_Proc, SG_Test_EpqErs_Ref
public :: SG_Mk_SDREO
#ifdef _TEST_CSFDET_
public :: test_sg_csfdet
#endif
public :: nPack
Public:: SGS, CIS, EXS, State_is_used
public :: SG_Checks
public :: SG_New_Tag, SG_No_Reuse
public :: SG_CSF_Occ, SG_MatEl_H, SG_MatEl_Init, SG_Level_Map, SG_Level_Ints
public :: SG_Map_CSFs, SG_Num_Setup
public :: SGCnfList, SG_Mk_Cnf_List, SG_Free_Cnf_List, SG_Cnf_Occ, SG_Select_Explicit, SG_Block_H, SG_Explicit_Block
public :: SG_Cnf_Avg_Diag
public :: CSF_Det_Expand, SG_CSF_StepVec, SG_MatEl_Tpq, SG_Tpq_Psi, SG_Tpq_Psi_Ref, SG_Xpq_Psi

end module SGUGA
