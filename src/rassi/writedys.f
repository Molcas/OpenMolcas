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
************************************************************************

!     Subroutine to correctly bunch together spin-free Dyson orbitals
!     and pass them to the interface for .DysOrb and .molden file
!     creation.
!     Heavily based on SODYSORB subroutine.

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
      DIMENSION DYSEN(NZ)
      DIMENSION AMPS(NZ)
      DIMENSION CMO(NZ*NZ)

!+++  J. Creutzberg, J. Norell  - 2018 (.DysOrb and .molden export )

! In principle we should here take into account the diagonalization
! matrix from RASSCF. This will only rarely be needed for regular
! Dyson calculations so we will leave it out for now.

      IF (NSYM.GT.1) THEN
       WRITE(6,*)""
       WRITE(6,*)"! Molden export of Dyson orbitals is "//
     & "currently not supported for calculations with symmetry !"
       RETURN
      END IF

      EN_IND=1
      CMO_IND=1
      DYSCIND=0
      ORBNUM=0
      DO JSTATE=1,DYSEXPSF
         DO ISTATE=JSTATE+1,NSTATE

!     For each initial state JSTATE we will gather all the obtained Dysorbs
!     and export to a shared .DysOrb file and .molden file if
!     requested
!     Each file can however only contain NZ numbe of orbitals, so we
!     might have to split into several files IFILE
          IF ( ISTATE.EQ.(JSTATE+1) ) THEN
              IFILE=1
              DYSCIND=0 ! Orbital coeff. index
              ORBNUM=0 ! Dysorb index for given JSTATE
              CMO=0.0D0
              DYSEN=0.0D0
              AMPS=0.0D0
          END IF

         IF (DYSAMPS(JSTATE,ISTATE).GT.1.0D-6) THEN
          DO NDUM=1,NZ
             DYSCIND=DYSCIND+1
             CMO(DYSCIND)=SFDYS(NDUM,ISTATE,JSTATE)
          END DO
          ORBNUM=ORBNUM+1
          DYSEN(ORBNUM)=ENERGY(ISTATE)-ENERGY(JSTATE)
          AMPS(ORBNUM)=DYSAMPS(JSTATE,ISTATE)*DYSAMPS(JSTATE,ISTATE)
         END IF

! Write the Dysorbs from JSTATE to .DysOrb and .molden file
! (Enough to fill one file)
         IF(ORBNUM.EQ.NZ) THEN
          Call Dys_Interf(0,JSTATE,IFILE,NZ,CMO,
     &        DYSEN,AMPS)
          IFILE=IFILE+1
          SODYSCIND=0 ! Orbital coeff. index
          ORBNUM=0 ! Dysorb index for given JSTATE
          SODYSCMO=0.0D0
          DYSEN=0.0D0
          AMPS=0.0D0
         END IF

       END DO ! ISTATE

! Write the Dysorbs from JSTATE to .DysOrb and .molden file
! (All remaining, if any)
        IF(ORBNUM.GT.0) THEN
        Call Dys_Interf(0,JSTATE,IFILE,NZ,CMO,
     &        DYSEN,AMPS)
        END IF
! +++


      END DO ! JSTATE

      RETURN
      END



