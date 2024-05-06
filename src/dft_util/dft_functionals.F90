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

module DFT_Functionals

use Definitions, only: iwp

implicit none
private

abstract interface
  subroutine DFT_FUNCTIONAL( &
#                           define _CALLING_
#                           include "dft_functional.fh"
                           )
    import :: iwp
#   include "dft_functional.fh"
  end subroutine DFT_FUNCTIONAL
end interface

public :: DFT_FUNCTIONAL, NDSD_Ts, NucAtt, Overlap

contains

#define _IN_MODULE_
#include "overlap.F90"
#include "nucatt.F90"
#include "ndsd_ts.F90"

end module DFT_Functionals
