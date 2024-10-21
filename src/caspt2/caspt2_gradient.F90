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
module caspt2_gradient

  use definitions, only: iwp,wp

  ! some gradient stuff
  integer(kind=iwp) :: iVecL           = 7_iwp ! Solution of the Lambda equation
  ! iVecG (G is probably gradient stuff) is used in
  ! caspt2_res.f to temporarily store residual vectors in solving the lambda equation
  ! sigder.f and clagx.f to temporarily store derivatives of overlap
  integer(kind=iwp) :: iVecG           = 8_iwp
  integer(kind=iwp) :: idSDMat(8,13)   = 0_iwp  ! offset of overlap derivative; can be defined with 11
  logical(kind=iwp) :: if_SSDM         = .false. ! State-dependent DM is used in Fock or not
  ! The state for which derivatives of the Lagrangian is computed.
  ! This is equivalent to jState
  integer(kind=iwp) :: jStLag          = 0_iwp

  ! unit numbers
  integer(kind=iwp) :: LuPT2           = 0_iwp
  integer(kind=iwp) :: LuGAMMA         = 0_iwp
  integer(kind=iwp) :: LuCMOPT2        = 0_iwp
  integer(kind=iwp) :: LuSTD           = 0_iwp
  integer(kind=iwp) :: LuAPT2          = 0_iwp
  integer(kind=iwp) :: LuPT2GRD        = 0_iwp

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
  integer(kind=iwp) :: IDSAVGRD        = 0_iwp
  integer(kind=iwp) :: idBoriMat(8,13) = 0_iwp
  real(kind=wp)     :: ConvInvar       = 0.0_wp

  ! whether PT2 energy is invariant wrt rotations among inactive
  ! and secondary orbitals
  logical(kind=iwp) :: if_invaria      = .true.

  ! some derivatives of Lagrangian etc.
  real(kind=wp),allocatable :: CLag(:,:)
  real(kind=wp),allocatable :: CLagFull(:,:)
  real(kind=wp),allocatable :: OLag(:)
  real(kind=wp),allocatable :: OLagFull(:)
  real(kind=wp),allocatable :: SLag(:,:)
  real(kind=wp),allocatable :: WLag(:)
  integer(kind=iwp) :: nCLag           = 0_iwp
  integer(kind=iwp) :: nOLag           = 0_iwp
  integer(kind=iwp) :: nSLag           = 0_iwp
  integer(kind=iwp) :: nWLag           = 0_iwp

  ! some correlated density matrices
  real(kind=wp),allocatable :: DPT2_tot(:)
  real(kind=wp),allocatable :: DPT2C_tot(:)
  real(kind=wp),allocatable :: DPT2_AO_tot(:)
  real(kind=wp),allocatable :: DPT2C_AO_tot(:)
  real(kind=wp),allocatable :: DPT2Canti_tot(:)

  ! Fock-related matrices
  real(kind=wp),allocatable :: FIMO_all(:)
  real(kind=wp),allocatable :: FIFA_all(:)
  real(kind=wp),allocatable :: FIFASA_all(:)

  ! natural <-> quasi-canonical transformation of frozen orbitals
  real(kind=wp),allocatable :: TraFro(:)

  ! derivative of the weight factor for XDW-CASPT2
  real(kind=wp),allocatable :: OMGDER(:,:)

  ! number of CI vectors per batch in mkfg3.f and derfg3.f
  integer(kind=iwp) :: nbuf1_grad      = 0_iwp

end module caspt2_gradient
