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
use setup,only:
use Real_Spherical,only:
use iSD_data,only:
use k2_arrays,only:
use rmat,only:
use define_af,only:
use Property_Label,only:
use ri_glob,only:

private

public :: DeDe_SCF, int_kernel, int_mem, int_wrout, OneEl_ij, OneEl_Inner, OneEl_Integrals, Integral_WrOut2, Integral_RI_3, &
          Integral_RICD, Integral_RI_2, Integral_WrOut_Cho, No_Routine, Integral_WrOut_Cho_Diag

#define _FIXED_FORMAT_
abstract interface
  subroutine int_kernel( &
#                       define _CALLING_
#                       include "int_interface.fh"
                       )
    use Index_Functions, only: nTri_Elem1
#   include "int_interface.fh"
  end subroutine int_kernel

  subroutine int_mem( &
#                    define _CALLING_
#                    include "mem_interface.fh"
                    )
#   include "mem_interface.fh"
  end subroutine int_mem

  subroutine int_wrout( &
#                      define _CALLING_
#                      include "int_wrout_interface.fh"
                      )
#   include "int_wrout_interface.fh"
  end subroutine int_wrout
end interface

procedure(int_wrout), pointer, public :: Int_postprocess => null()

contains

#define _IN_MODULE_
#include "oneel_ij.F90"
#include "oneel_inner.F90"
#include "oneel_integrals.F90"
#include "dede_scf.F90"
#include "integral_wrout2.F90"
#include "no_routine.F90"
#include "../ri_util/integral_ri_3.F90"
#include "../ri_util/integral_ri_2.F90"
#include "../ri_util/integral_ricd.F90"
#include "../cholesky_util/integral_wrout_cho.F90"
#include "../cholesky_util/integral_wrout_cho_diag.F90"

end module Integral_Interfaces
