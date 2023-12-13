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

module Fock_util_global

use Definitions, only: wp, iwp

implicit none
private

!=======================================================================
! Contents of old include files (programs that used them):
!-----------------------------------------------------------------------
! - choscf.fh (SCF)
!   DECO
!-----------------------------------------------------------------------
! - chounit.fh (SCF)
!   Lunit
!-----------------------------------------------------------------------
! - chodensity.fh (SCF, RASSCF, MCPDFT, CASPT2)
!   DensityCheck
!-----------------------------------------------------------------------
! - choscreen.fh (SCF, RASSCF, MCPDFT, LOCALISATION)
!   Estimate, Update
!-----------------------------------------------------------------------
! - chlcas.fh (RASSCF, MCPDFT)
!   ALGO, DoCholesky
!-----------------------------------------------------------------------
! - cholk.fh (RASSCF, MCPDFT)
!   Deco, dmpk, DoLocK Nscreen
!-----------------------------------------------------------------------
! - chotodo.fh (RASSCF, MCPDFT)
!   DoActive
!-----------------------------------------------------------------------
! - choras.fh (CASPT2)
!   ALGO, DECO, REORD
!-----------------------------------------------------------------------
! - chomclr.fh (MCLR)
!   Deco, dmpk, Estimate, NScreen, Update
!-----------------------------------------------------------------------
! - lkscreen.fh (RASSI)
!   Deco, Estimate, PseudoChoMOs, Update
!-----------------------------------------------------------------------
! - cho_jobs.fh (RASSI)
!   Fake_CMO2
!=======================================================================

integer(kinD=iwp) :: ALGO = 1, Lunit(8) = -1, Nscreen = 10
real(kind=wp) :: dmpk = 0.1_wp
logical(kind=iwp) :: Deco = .true., DensityCheck = .false., DoActive = .true., DoCholesky = .false., DoLocK = .true., &
                     Estimate = .false., Fake_CMO2 = .false., PseudoChoMOs = .false., REORD = .false., Update = .true.

public :: ALGO, Deco, DensityCheck, dmpk, DoActive, DoCholesky, DoLocK, Estimate, Fake_CMO2, Lunit, Nscreen, PseudoChoMOs, REORD, &
          Update

end module Fock_util_global
