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

! Global variables of the CASPT2 module
! TODO: move here all variables from CASPT2 common blocks defined in caspt2.fh
module caspt2_global

  use definitions, only: iwp,wp

  real(kind=wp)     :: ipea_shift = 0.0_wp
  real(kind=wp)     :: imag_shift = 0.0_wp
  real(kind=wp)     :: real_shift = 0.0_wp

  ! sigma-p regularization
  real(kind=wp)     :: sigma_p_epsilon  = 0.0_wp
  integer(kind=iwp) :: sigma_p_exponent = 2_iwp

end module caspt2_global