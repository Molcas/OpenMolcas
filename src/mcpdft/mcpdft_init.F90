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
!> Sets values in common blocks in rasscf.fh, general.fh, timers.fh
!***********************************************************************

subroutine mcpdft_init()
  Use Fock_util_global,only:DoCholesky
  Use Cholesky,only:ChFracMem
  Use KSDFT_Info,Only:CoefR,CoefX
  use mcpdft_output,only:set_print_level

  implicit none

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "timers.fh"

  integer i
!----------------------------------------------------------------------*

! Set print levels, and adjust them if needed:
  call set_print_level()

! Cholesky-related settings:
  Call DecideOnCholesky(DoCholesky)

#ifdef _MOLCAS_MPP_
  ChFracMem = 0.3d0
#else
  ChFracMem = 0.0d0
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
  WEIGHT = 0.0d0
  WEIGHT(1) = 1.0D0
! iteration energies
  ENER = 0.0D0
! prethr: energy threshold for printout of orbitals
  prethr = 0.15d0

! Default value for type of CASSCF (used for DFT)
  DFTFOCK = "ROKS"
  ExFac = 0.0d0

! Initialize KSDFT coefficients (S Dong, 2018)
  CoefR = 1.0D0
  CoefX = 1.0D0

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
  DO I = 1,mxSym
    NFRO(I) = 0
    NISH(I) = 0
    NASH(I) = 0
    NRS1(I) = 0
    NRS2(I) = 0
    NRS3(I) = 0
    NSSH(I) = 0
    NDEL(I) = 0
    NBAS(I) = 0
  ENDDO

! initialize occupation numbers for GAS
  NGAS = 3
  NGSSH = 0
  IGSOCCX = 0
  do I = 1,mxOrb
    IXSYM(I) = 0
  enddo

!     Auxiliary vector ITRI(I)=I*(I-1)/2
  do I = 2,ITRIM
    ITRI(I) = ITRI(I-1)+I-1
  enddo

! Initialize Timing Variables
  Ebel_3 = 0.0d0
  Eterna_3 = 0.0d0
  Rado_3 = 0.0d0
  Rolex_3 = 0.0d0
  Omega_3 = 0.0d0
  Tissot_3 = 0.0d0
  Piaget_3 = 0.0d0
  Candino_3 = 0.0d0
  Fortis_3 = 0.0d0
  Zenith_3 = 0.0d0
  Gucci_3 = 0.0d0
  Alfex_3 = 0.0d0
  WTC_3 = 0.0d0
  Longines_3 = 0.0d0
  Oris_2 = 0.0d0
  Movado_2 = 0.0d0

END
