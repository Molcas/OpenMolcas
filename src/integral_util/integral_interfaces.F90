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

module Integral_Interfaces

! Dummy modules to get correct order of compilation
use IOBuf, only:
use Real_Spherical, only:
use k2_arrays, only:
use rmat, only:
use define_af, only:
use Property_Label, only:
use Basis_Info, only:
use Center_Info, only:
use Gateway_Info, only:
use Sizes_of_Seward, only:
use Symmetry_Info, only:
use SOAO_Info, only:
use Dens_Stuff, only:
use stdalloc, only:
use Constants, only:

use iSD_data, only: nSD
use Index_Functions, only: nTri_Elem1
use k2_structure, only: k2_type
use Definitions, only: wp, iwp

implicit none
private

abstract interface
  subroutine int_kernel( &
#                       define _CALLING_
#                       include "int_interface.fh"
                       )
    import :: nTri_Elem1, wp, iwp
#   include "int_interface.fh"
  end subroutine int_kernel

  subroutine int_mem( &
#                    define _CALLING_
#                    include "mem_interface.fh"
                    )
    import :: iwp
#   include "mem_interface.fh"
  end subroutine int_mem

  subroutine int_wrout( &
#                      define _CALLING_
#                      include "int_wrout_interface.fh"
                      )
    import :: wp, iwp
#   include "int_wrout_interface.fh"
  end subroutine int_wrout

  subroutine prm_kernel( &
#                       define _CALLING_
#                       include "prm_interface.fh"
                       )
    import :: nTri_Elem1, wp, iwp
#   include "prm_interface.fh"
  end subroutine prm_kernel

  subroutine twoel_kernel( &
#                         define _CALLING_
#                         include "twoel_interface.fh"
                         )
    import :: wp, iwp, k2_type, nSD
#   include "twoel_interface.fh"
  end subroutine twoel_kernel
end interface

! Intel 13 compiler only works with "public" here
procedure(int_wrout), public, pointer :: Int_postprocess => null()

public :: DeDe_SCF, Drv2El_dscf, FckAcc, int_kernel, int_mem, int_wrout, OneEl_ij, OneEl_Integrals, prm_kernel, twoel_kernel

contains

! Subroutines that need an explicit interface (target or allocatable arguments)
#define _IN_MODULE_
#include "dede_scf.F90"
#include "drv2el_dscf.F90"
#include "fckacc.F90"
#include "oneel_ij.F90"
#include "oneel_integrals.F90"

end module Integral_Interfaces
