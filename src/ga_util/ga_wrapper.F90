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

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

module GA_Wrapper

implicit none
private

! This module offers a wrapper around the GlobalArrays include files.
! Other routines should use this module instead of the include files.

#include "global.fh"
#include "mafdecls.fh"

! External procedures should forgo TKR checks
external :: GA_Brdcst

public :: DBL_MB, MT_BYTE, MT_DBL, MT_INT
public :: GA_Brdcst, GA_Create, GA_Create_Irreg, GA_DDot, GA_Destroy, GA_Duplicate, GA_NNodes, GA_NodeId, GA_Read_Inc

end module GA_Wrapper

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GA_Wrapper)

#endif
