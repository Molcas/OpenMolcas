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

  private

  real(kind=wp)     , public:: ipea_shift = 0.0_wp
  real(kind=wp)     , public:: imag_shift = 0.0_wp
  real(kind=wp)     , public:: real_shift = 0.0_wp

  ! sigma-p regularization
  real(kind=wp)     , public:: sigma_p_epsilon  = 0.0_wp
  integer(kind=iwp) , public:: sigma_p_exponent = 2_iwp


  ! some gradient stuff
  integer(kind=iwp) , public:: iVecL           = 7_iwp ! Solution of the Lambda equation
  ! iVecG (G is probably gradient stuff) is used in
  ! caspt2_res.f to temporarily store residual vectors in solving the lambda equation
  ! sigder.f and clagx.f to temporarily store derivatives of overlap
  integer(kind=iwp) , public:: iVecG           = 8_iwp
  integer(kind=iwp) , public:: idSDMat(8,13)   = 0_iwp  ! offset of overlap derivative; can be defined with 11
  logical(kind=iwp) , public:: if_SSDM         = .false. ! State-dependent DM is used in Fock or not
  ! The state for which derivatives of the Lagrangian is computed.
  ! This is equivalent to jState
  integer(kind=iwp) , public:: jStLag          = 0_iwp

  ! unit numbers
  integer(kind=iwp) , public:: LuPT2           = 0_iwp
  integer(kind=iwp) , public:: LuGAMMA         = 0_iwp
  integer(kind=iwp) , public:: LuCMOPT2        = 0_iwp
  integer(kind=iwp) , public:: LuSTD           = 0_iwp
  integer(kind=iwp) , public:: LuAPT2          = 0_iwp
  integer(kind=iwp) , public:: LuPT2GRD        = 0_iwp

  ! gradients and NAC switches
  logical(kind=iwp) , public:: do_grad         = .false.
  logical(kind=iwp) , public:: do_nac          = .false.
  logical(kind=iwp) , public:: do_csf          = .false. ! CSF term in deriv. coup.
  integer(kind=iwp) , public:: iRoot1          = 0_iwp
  integer(kind=iwp) , public:: iRoot2          = 0_iwp
  integer(kind=iwp) , public:: nStpGrd         = 1_iwp


  ! for removing the weired loop
  integer(kind=iwp) , public:: iStpGrd         = 1_iwp
  integer(kind=iwp) , public:: LUGRAD          = 0_iwp

  ! for IPEA
  logical(kind=iwp) , public:: do_lindep       = .false.
  logical(kind=iwp) , public:: if_invar        = .true. ! active invariance
  integer(kind=iwp) , public:: IDSAVGRD        = 0_iwp
  integer(kind=iwp) , public:: idBoriMat(8,13) = 0_iwp
  real(kind=wp)     , public:: ConvInvar       = 0.0_wp

  ! whether PT2 energy is invariant wrt rotations among inactive
  ! and secondary orbitals
  logical(kind=iwp) , public:: if_invaria      = .true.

  ! some derivatives of Lagrangian etc.
  real(kind=wp),allocatable , public:: CLag(:,:)
  real(kind=wp),allocatable , public:: CLagFull(:,:)
  real(kind=wp),allocatable , public:: OLag(:)
  real(kind=wp),allocatable , public:: OLagFull(:)
  real(kind=wp),allocatable , public:: SLag(:,:)
  real(kind=wp),allocatable , public:: WLag(:)
  integer(kind=iwp) , public:: nCLag           = 0_iwp
  integer(kind=iwp) , public:: nOLag           = 0_iwp
  integer(kind=iwp) , public:: nSLag           = 0_iwp
  integer(kind=iwp) , public:: nWLag           = 0_iwp

  ! some correlated density matrices
  real(kind=wp),allocatable , public:: DPT2_tot(:)
  real(kind=wp),allocatable , public:: DPT2C_tot(:)
  real(kind=wp),allocatable , public:: DPT2_AO_tot(:)
  real(kind=wp),allocatable , public:: DPT2C_AO_tot(:)
  real(kind=wp),allocatable , public:: DPT2Canti_tot(:)

  ! Fock-related matrices
  real(kind=wp),allocatable , public:: FIMO_all(:)
  real(kind=wp),allocatable , public:: FIFA_all(:)
  real(kind=wp),allocatable , public:: FIFASA_all(:)

  ! natural <-> quasi-canonical transformation of frozen orbitals
  real(kind=wp),allocatable , public:: TraFro(:)

  ! derivative of the weight factor for XDW-CASPT2
  real(kind=wp),allocatable , public:: OMGDER(:,:)

  ! number of CI vectors per batch in mkfg3.f and derfg3.f
  integer(kind=iwp) , public:: nbuf1_grad      = 0_iwp

end module caspt2_global
