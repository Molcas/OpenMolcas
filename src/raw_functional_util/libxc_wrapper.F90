! ************************************************************************
! * This file is part of OpenMolcas.                                     *
! *                                                                      *
! * OpenMolcas is free software; you can redistribute it and/or modify   *
! * it under the terms of the GNU Lesser General Public License, v. 2.1. *
! * OpenMolcas is distributed in the hope that it will be useful, but it *
! * is provided "as is" and without any express or implied warranties.   *
! * For more details see the full text of the license in the file        *
! * LICENSE or in <http://www.gnu.org/licenses/>.                        *
! *                                                                      *
! * Copyright (C) 2022 Molecular Sciences Software Institute             *
! ************************************************************************

module libxc_wrapper
  use, intrinsic :: iso_c_binding
  implicit none
  interface
     ! Query libxc if the functional is a mixed functional
     integer(c_int) function libxc_num_aux_funcs(id) bind(C)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: id
     end function libxc_num_aux_funcs

     ! Query libxc for the id numbers and weights of the mixture, e.g.
     ! B3LYP = 0.08 Slater + 0.72 Becke'88 + 0.19 VWN5 + 0.81 LYP
     subroutine libxc_get_aux_funcs(id, ids, weights) bind(C)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: id
       integer(c_int), pointer :: ids(:)
       real(c_double), pointer :: weights(:)
     end subroutine libxc_get_aux_funcs
  end interface
end module libxc_wrapper
