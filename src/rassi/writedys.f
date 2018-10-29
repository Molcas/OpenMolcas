************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
!     +++ J. Creutzberg, J. Norell  - 2018
!     Subroutine to correctly bunch together spin-free Dyson orbitals
!     and pass them to the interface for .DysOrb and .molden file
!     creation.
!     Heavily based on SODYSORB subroutine.

      SUBROUTINE WRITEDYS(DYSAMPS,SFDYS,NZ,NSTATE2,ENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"
      CHARACTER*16 ROUTINE

      INTEGER   NZ,ORBNUM
      INTEGER   IDISK
      INTEGER   DYSCIND,ENIND
      INTEGER   INDJ,INDI,SFI,SFJ,ZI,ZJ,NSZZ,NDUM

      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SFDYS(NZ,NSTATE,NSTATE)
      DIMENSION ENERGY(NSTATE)
      DIMENSION DYSEN(NZ)
      DIMENSION AMPS(NZ)
      DIMENSION CMO(NZ*NZ)

C Read in all the previously saved SF Dyson orbitals in the
C atomic basis from disk

! In principle we should here take into account the diagonalization
! matrix from RASSCF. This will only rarely be needed for regular
! Dyson calculations so we will leave it out for now.

      EN_IND=1
      CMO_IND=1
      DO JSTATE=1,NSTATE

       IF (JSTATE.GT.DYSEXPSF) THEN
        EXIT
       END IF

         DO ISTATE=JSTATE+1,NSTATE
            IF (DYSAMPS(JSTATE,ISTATE).GT.1.0D-6) THEN

             IDISK=IWORK(LIDDYS+(ISTATE-1)*NSTATE+JSTATE-1)
             CALL DDAFILE(LUDYS,2,SFDYS(:,JSTATE,ISTATE),NZ,IDISK)
             ! Loops over states are performed triangularly, but
             ! permutation of degenerate states in the SO part
             ! might 'escape' this, therefore we fill out the
             ! full matrix to be safe.
             IDISK=IWORK(LIDDYS+(ISTATE-1)*NSTATE+JSTATE-1)
             CALL DDAFILE(LUDYS,2,SFDYS(:,ISTATE,JSTATE),NZ,IDISK)
            END IF

!+++  J. Creutzberg, J. Norell  - 2018 (.DysOrb and .molden export )
!     For each initial state JSTATE we will gather all the obtained Dysorbs
!     and export to a shared .DysOrb file and .molden file if
!     requested
          IF ( ISTATE.EQ.(JSTATE+1) ) THEN
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

       END DO ! ISTATE

! Write the Dysorbs from JSTATE to .DysOrb and .molden file
         Call Dys_Interf(.FALSE.,JSTATE,NZ,CMO,
     &        DYSEN,AMPS)
! +++


      END DO ! JSTATE

      RETURN
      END



