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
#ifdef _EFP_

! File from the EFP submodule
#include "efp.f90"

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"

module EFP
contains
dummy_empty_procedure(efp_mod)
end module EFP

#endif
