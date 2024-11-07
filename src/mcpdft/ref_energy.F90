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

subroutine ref_energy(mcscf_energy,nroots)
  use definitions,only:iwp,wp
  use mcpdft_input,only:mcpdft_options
  use mspdft,only:heff
  implicit none

  integer(kind=iwp),intent(in) :: nroots
  real(kind=wp),dimension(nroots),intent(out) :: mcscf_energy

  integer(kind=iwp) :: root

  if(mcpdft_options%mspdft) then
    do root = 1,nroots
      mcscf_energy(root) = heff(root,root)
    enddo
  endif

endsubroutine
