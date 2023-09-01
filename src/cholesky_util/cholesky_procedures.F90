!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! This module contains procedures that need an interface
module Cholesky_procedures

private
public :: Build_Mp2Dens, Cho_CGM_InfVec, Cho_GetLQ, Cho_P_GetLQ, Cho_VecBuf_Subtr, Cho_X_GetIP_InfVec

contains

#define _IN_MODULE_
#include "build_mp2dens.F90"
#include "cho_cgm_infvec.F90"
#include "cho_getlq.F90"
#include "cho_p_getlq.F90"
#include "cho_vecbuf_getlq.F90"
#include "cho_vecbuf_subtr.F90"
#include "cho_x_getip_infvec.F90"

end module Cholesky_procedures
