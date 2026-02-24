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
!> Sets values in the modules timers, rasscf_global and general_data.
!***********************************************************************

subroutine mcpdft_init()

use Fock_util_global, only: DoCholesky
use Cholesky, only: ChFracMem
use timers, only: TimeAoMo, TimeCIOpt, TimeDavid, TimeDens, TimeFock, TimeHCSCE, TimeHDiag, TimeHSel, TimeInput, TimeOrb, &
                  TimePage, TimeRelax, TimeSigma, TimeTotal, TimeTrans, TimeWfn
use mcpdft_output, only: set_print_level
use general_data, only: ispin, nactel, nash, nbas, ndel, nelec3, nfro, nhole1, nish, nrs1, nrs2, nrs3, nssh, stsym
use rasscf_global, only: DFTFOCK, ExFac, IPT2, iroot, iXSym, lROOTS, NonEq, NROOTS, TITLE, weight
use Constants, only: Zero, One
#ifdef _MOLCAS_MPP_
use Definitions, only: wp
#endif

implicit none

!----------------------------------------------------------------------*

! Set print levels, and adjust them if needed:
call set_print_level()

! Cholesky-related settings:
call DecideOnCholesky(DoCholesky)

#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Zero
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
WEIGHT = Zero
WEIGHT(1) = One

! Default value for type of CASSCF (used for DFT)
DFTFOCK = 'ROKS'
ExFac = Zero

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
NonEq = .false.

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

! Initialize Timing Variables
TimeTotal = Zero
TimeInput = Zero
TimeWfn = Zero
TimeDens = Zero
TimeSigma = Zero
TimeHSel = Zero
TimeHDiag = Zero
TimeFock = Zero
TimeAoMo = Zero
TimeTrans = Zero
TimeCIOpt = Zero
TimeOrb = Zero
TimeDavid = Zero
TimePage = Zero
TimeHCSCE = Zero
TimeRelax = Zero

end subroutine mcpdft_init
