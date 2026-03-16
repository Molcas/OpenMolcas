!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2018, Jesper Norell                                    *
!               2018, Joel Creutzberg                                  *
!               2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

! Subroutine to correctly bunch together spin-free Dyson orbitals
! and pass them to the molden_dysorb interface for .molden export

! IFG: Added DysOrb export, not sure it's correct

subroutine WRITEDYS(DYSAMPS,SFDYS,NZ,ENERGY)

use Cntrl, only: NSTATE, DYSEXPSF
use Symmetry_Info, only: nSym => nIrrep
use rassi_data, only: NBASF
use Constants, only: Zero
use Definitions, only: wp

implicit none
integer NZ
real*8 DYSAMPS(NSTATE,NSTATE)
real*8 SFDYS(NZ,NSTATE,NSTATE)
real*8 ENERGY(NSTATE)
integer DYSCIND
integer NDUM
integer ORBNUM
real*8 DYSEN(NSTATE)
real*8 AMPS(NSTATE)
real*8 CMO(NZ*NSTATE)
character(len=30) Filename
character(len=80) TITLE
integer ISTATE, JSTATE, LUNIT
integer, external :: IsFreeUnit

!+++  J. Creutzberg, J. Norell  - 2018 (.molden export )

! In principle we should here take into account the diagonalization
! matrix from RASSCF. This will only rarely be needed for regular
! Dyson calculations so we will leave it out for now.

do JSTATE=1,DYSEXPSF

  ! For each initial state JSTATE up to DYSEXPSF we will gather all the obtained Dysorbs
  ! and export to a shared .molden file
  DYSCIND = 0 ! Orbital coeff. index
  ORBNUM = 0 ! Dysorb index for given JSTATE
  CMO = Zero ! Orbital coefficients
  DYSEN = Zero ! Orbital energies
  AMPS = Zero ! Transition amplitudes (shown as occupations)

  do ISTATE=JSTATE+1,NSTATE

    if (DYSAMPS(JSTATE,ISTATE) > 1.0e-5_wp) then
      do NDUM=1,NZ
        DYSCIND = DYSCIND+1
        CMO(DYSCIND) = SFDYS(NDUM,ISTATE,JSTATE)
      end do
      ORBNUM = ORBNUM+1
      DYSEN(ORBNUM) = ENERGY(ISTATE)-ENERGY(JSTATE)
      AMPS(ORBNUM) = DYSAMPS(JSTATE,ISTATE)*DYSAMPS(JSTATE,ISTATE)
    end if

  end do ! ISTATE

  ! If at least one orbital was found, export it/them
  if (ORBNUM > 0) then
    write(filename,'(A,I0)') 'MD_DYS.SF.',JSTATE
    call Molden_DysOrb(filename,DYSEN,AMPS,CMO,ORBNUM,NZ)

    write(filename,'(A,I0)') 'DYSORB.SF.',JSTATE
    LUNIT = IsFreeUnit(50)
    write(TITLE,'(A,I0)') '* Spin-free Dyson orbitals for state ',JSTATE
    call WRVEC_DYSON(filename,LUNIT,NSYM,NBASF,ORBNUM,CMO,AMPS,DYSEN,trim(TITLE),NZ)
    close(LUNIT)
  end if

end do ! JSTATE

end subroutine WRITEDYS
