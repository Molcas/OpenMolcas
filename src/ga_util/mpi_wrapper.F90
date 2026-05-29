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

module MPI_Wrapper

! This module offers a wrapper around the MPI module.
! Other routines should use this module instead of the bare MPI

use MPI, only: MPI_ADDRESS_KIND, MPI_COMM_WORLD, MPI_INTEGER, MPI_INTEGER4, MPI_INTEGER8, MPI_LOGICAL, MPI_REAL8

implicit none
private

! Not all implementations provide usable interfaces in the old-style MPI module,
! and not all implementations provide an usable MPI_F08 module
! External procedures should forgo TKR checks
external :: MPI_AllGather, MPI_AllGatherV, MPI_AllToAll, MPI_AllToAllV, MPI_Bcast

public :: MPI_ADDRESS_KIND, MPI_COMM_WORLD, MPI_INTEGER, MPI_INTEGER4, MPI_INTEGER8, MPI_LOGICAL, MPI_REAL8
public :: MPI_AllGatherV, MPI_AllToAll, MPI_AllToAllV, MPI_Bcast

end module MPI_Wrapper

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(MPI_Wrapper)

#endif
