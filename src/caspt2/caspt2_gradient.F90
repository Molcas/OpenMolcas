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
  logical(kind=iwp) :: do_grad         = .false.
  logical(kind=iwp) :: do_nac          = .false.
  logical(kind=iwp) :: do_csf          = .false. ! CSF term in deriv. coup.
  integer(kind=iwp) :: iRoot1          = 0_iwp
  integer(kind=iwp) :: iRoot2          = 0_iwp
  integer(kind=iwp) :: nStpGrd         = 1_iwp

  ! for removing the weired loop
  integer(kind=iwp) :: iStpGrd         = 1_iwp
  integer(kind=iwp) :: LUGRAD          = 0_iwp

  ! for IPEA
  logical(kind=iwp) :: do_lindep       = .false.
  logical(kind=iwp) :: if_invar        = .true. ! active invariance
  integer(kind=iwp) :: LUSTD           = 0_iwp
  integer(kind=iwp) :: IDSAVGRD        = 0_iwp
  integer(kind=iwp) :: idBoriMat(8,13) = 0_iwp
  real(kind=wp)     :: ConvInvar       = 0.0_wp

  ! whether CASPT2 energy is invariant wrt rotations among inactive
  ! and secondary orbitals
  logical(kind=iwp) :: if_invaria      = .true.

  ! natural <-> quasi-canonical transformation of frozen orbitals
  real(kind=wp),allocatable :: TraFro(:)

  ! number of CI vectors per batch in mkfg3.f and derfg3.f
  integer(kind=iwp) :: nbuf1_grad      = 0_iwp

end module caspt2_gradient
