************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2018, Jesper Norell                                    *
*               2018, Joel Creutzberg                                  *
*               2023, Ignacio Fdez. Galvan                             *
************************************************************************

!     Subroutine to correctly bunch together spin-free Dyson orbitals
!     and pass them to the molden_dysorb interface for .molden export

!     IFG: Added DysOrb export, not sure it's correct

      SUBROUTINE WRITEDYS(DYSAMPS,SFDYS,NZ,ENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"

      INTEGER   NZ,ORBNUM
      INTEGER   DYSCIND
      INTEGER   NDUM

      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SFDYS(NZ,NSTATE,NSTATE)
      DIMENSION ENERGY(NSTATE)
      DIMENSION DYSEN(NSTATE)
      DIMENSION AMPS(NSTATE)
      DIMENSION CMO(NZ*NSTATE)

      Character*30 Filename
      Character*80 TITLE

!+++  J. Creutzberg, J. Norell  - 2018 (.molden export )

! In principle we should here take into account the diagonalization
! matrix from RASSCF. This will only rarely be needed for regular
! Dyson calculations so we will leave it out for now.

      DO JSTATE=1,DYSEXPSF

!     For each initial state JSTATE up to DYSEXPSF we will gather all the obtained Dysorbs
!     and export to a shared .molden file
         DYSCIND=0 ! Orbital coeff. index
         ORBNUM=0 ! Dysorb index for given JSTATE
         CMO=0.0D0 ! Orbital coefficients
         DYSEN=0.0D0 ! Orbital energies
         AMPS=0.0D0 ! Transition amplitudes (shown as occupations)

         DO ISTATE=JSTATE+1,NSTATE

         IF (DYSAMPS(JSTATE,ISTATE).GT.1.0D-5) THEN
          DO NDUM=1,NZ
             DYSCIND=DYSCIND+1
             CMO(DYSCIND)=SFDYS(NDUM,ISTATE,JSTATE)
          END DO
          ORBNUM=ORBNUM+1
          DYSEN(ORBNUM)=ENERGY(ISTATE)-ENERGY(JSTATE)
          AMPS(ORBNUM)=DYSAMPS(JSTATE,ISTATE)*DYSAMPS(JSTATE,ISTATE)
         END IF

       END DO ! ISTATE

! If at least one orbital was found, export it/them
        IF(ORBNUM.GT.0) THEN
         Write(filename,'(A,I0)') 'MD_DYS.SF.',JSTATE
         Call Molden_DysOrb(filename,DYSEN,AMPS,CMO,ORBNUM,NZ)

         Write(filename,'(A,I0)') 'DYSORB.SF.',JSTATE
         LUNIT=IsFreeUnit(50)
         Write(TITLE,'(A,I0)') '* Spin-free Dyson orbitals for state ',
     &                         JSTATE
         Call WRVEC_DYSON(filename,LUNIT,NSYM,NBASF,ORBNUM,CMO,AMPS,
     &                    DYSEN,Trim(TITLE),NZ)
         Close(LUNIT)
        END IF

      END DO ! JSTATE

      RETURN
      END
