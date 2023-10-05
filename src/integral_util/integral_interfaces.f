************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Module Integral_Interfaces

      ! Dummy modules to get correct order of compilation
      Use Real_Spherical, only:
      Use iSD_data, only:
      Use k2_arrays, only:

      Private

      Public :: DeDe_SCF, int_kernel, int_mem, int_wrout, OneEl_ij,
     &          OneEl_Inner, OneEl_Integrals

#define _FIXED_FORMAT_
      Abstract interface
        Subroutine int_kernel(
#                             define _CALLING_
#                             include "int_interface.fh"
     &                       )
        Use Index_Functions, only: nTri_Elem1
#       include "int_interface.fh"
        End Subroutine int_kernel

        Subroutine int_mem(
#                          define _CALLING_
#                          include "mem_interface.fh"
     &                    )
#       include "mem_interface.fh"
        End Subroutine int_mem

        Subroutine int_wrout(
#                            define _CALLING_
#                            include "int_wrout_interface.fh"
     &                      )
#       include "int_wrout_interface.fh"
        End Subroutine int_wrout
      End interface

      Contains

#define _IN_MODULE_
#include "oneel_ij.f90"
#include "oneel_inner.f90"
#include "oneel_integrals.f"
#include "dede_scf.f"

      End Module Integral_Interfaces
