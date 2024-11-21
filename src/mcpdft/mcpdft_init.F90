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
!  mcpdft_init
!
!> @brief
!>   Initialize variables in commons, and set default values.
!>   Determine whether orbital files should be read, etc.
!>
!> @details
!> Sets values in common blocks in timers.fh and the modules
!> rasscf_global.F90 and general_data.F90
!***********************************************************************

subroutine mcpdft_init()
  use constants,only:zero,one
  use definitions,only:iwp
  Use Fock_util_global,only:DoCholesky
  Use Cholesky,only:ChFracMem
  use mcpdft_output,only:set_print_level
  use general_data,only:ispin,nactel,nelec3,nhole1,stsym,nfro,nish,nash,nrs1,nrs2,nrs3,nssh,ndel,nbas
  use rasscf_global,only:iroot,weight,DFTFOCK,ExFac,IPT2,iTRIM,lROOTS,NonEq,NROOTS,TITLE,iXSym,iTRI

  implicit none

#include "timers.fh"

  integer(kind=iwp) :: i
!----------------------------------------------------------------------*

! Set print levels, and adjust them if needed:
  call set_print_level()

! Cholesky-related settings:
  Call DecideOnCholesky(DoCholesky)

#ifdef _MOLCAS_MPP_
  ChFracMem = 0.3d0
#else
  ChFracMem = zero
#endif

! Default title line:
  TITLE(1) = '(No title given)'

! number of roots required in CI
  NROOTS = 1
! number of roots actually used in CI-DAVIDSON
  LROOTS = 1
! sequence numbers for roots in CI counted from lowest energy.
  iRoot = 0
  IROOT(1) = 1
! weights used for average energy calculations
  WEIGHT = zero
  WEIGHT(1) = one

! Default value for type of CASSCF (used for DFT)
  DFTFOCK = "ROKS"
  ExFac = zero

! default spin value (singlet)
  ISPIN = 1
! default symmetry
  STSYM = 1
! default number of active electrons
  NACTEL = 0
! default maximum number of holes in RAS1
  NHOLE1 = 0
! default maximum number of electrons in RAS3
  NELEC3 = 0
! This run will not be the start for a CASPT2 calculation
  IPT2 = 0

! These keys will activate the calculation of the high
! frequency contribution to the reaction field
! ???
! This key controls if a non-equilibrium reaction field
! calculation is performed.
  NonEq = .False.

! set default values for orbitals
  nFro(:) = 0
  nIsh(:) = 0
  nAsh(:) = 0
  nRs1(:) = 0
  NRS2(:) = 0
  NRS3(:) = 0
  NSSH(:) = 0
  NDEL(:) = 0
  NBAS(:) = 0
  ixsym(:) = 0

!     Auxiliary vector ITRI(I)=I*(I-1)/2
  do I = 2,ITRIM
    ITRI(I) = ITRI(I-1)+I-1
  enddo

! Initialize Timing Variables
  Ebel_3 = zero
  Eterna_3 = zero
  Rado_3 = zero
  Rolex_3 = zero
  Omega_3 = zero
  Tissot_3 = zero
  Piaget_3 = zero
  Candino_3 = zero
  Fortis_3 = Zero
  Zenith_3 = zero
  Gucci_3 = zero
  Alfex_3 = zero
  WTC_3 = zero
  Longines_3 = zero
  Oris_2 = zero
  Movado_2 = Zero

END
