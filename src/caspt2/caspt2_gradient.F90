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

! Variables related to CASPT2 gradients
! TODO: move here everything that is in caspt2_grad.fh
module caspt2_gradient

  use definitions, only: iwp,wp

  ! gradients and NAC switches
  logical(kind=iwp) :: do_grad = .false.
  logical(kind=iwp) :: do_nac  = .false.
  logical(kind=iwp) :: do_csf  = .false. ! CSF term in deriv. coup.
  integer(kind=iwp) :: iRoot1  = 0_iwp
  integer(kind=iwp) :: iRoot2  = 0_iwp
  integer(kind=iwp) :: nStpGrd = 1_iwp

end module caspt2_gradient
