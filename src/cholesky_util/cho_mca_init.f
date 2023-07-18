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
      SUBROUTINE CHO_MCA_INIT(SKIP_PRESCREEN)
!
!     Purpose: initialization of Cholesky decomposition in MOLCAS.
!
      use index_arrays, only: iSO2Sh
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iSP2F, iShlSO,  &
     &                  iShP2RS, iShP2Q
      use stdalloc
      Implicit Real*8 (a-h,o-z)
      LOGICAL SKIP_PRESCREEN
#include "cholesky.fh"
#include "choorb.fh"

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_INIT')

!     INTEGER  CHO_ISAO
!     EXTERNAL CHO_ISAO
!     INTEGER  CHO_ISAOSH
!     EXTERNAL CHO_ISAOSH

!     Check that the number of shells is within limits.
!     -------------------------------------------------

      IF (NSHELL .LT. 1) THEN
         WRITE(LUPRI,*) 'NSHELL out of bounds: ',NSHELL
         CALL CHO_QUIT('NSHELL out of bounds in '//SECNAM,102)
      END IF

!     Compute total #shell pair.
!     --------------------------

      NNSHL_TOT = NSHELL*(NSHELL + 1)/2
      IF (NNSHL_TOT .LT. 1) THEN
         WRITE(LUPRI,*)                                                 &
     &   'NNSHL_TOT=NSHELL*(NSHELL+1)/2 is non-positive: ',             &
     &   NNSHL_TOT
         WRITE(LUPRI,*) 'Integer overflow ?'
         CALL CHO_QUIT('NNSHL_TOT out of bounds in '//SECNAM,102)
      END IF

!     Compute contributing #shell pair by prescreening (if requested).
!     iSP2F(k): returns global shell pair of contributing shell pair k.
!               (Allocated and defined in CHO_DIAPS.)
!     -----------------------------------------------------------------

      IF (SKIP_PRESCREEN) THEN
         IF (NNSHL.LT.1 .OR. NNSHL.GT.NNSHL_TOT) THEN
            WRITE(LUPRI,*) SECNAM,': flag SKIP_PRESCREEN is ',          &
     &                     SKIP_PRESCREEN
            WRITE(LUPRI,*) 'NNSHL is out-of-bounds: ',NNSHL
            WRITE(LUPRI,*) 'Condition: 0 < NNSHL < ',NNSHL_TOT
            CALL CHO_QUIT('Initialization error in '//SECNAM,102)
         END IF
         IF (SIZE(iSP2F) .NE. NNSHL) THEN
            WRITE(LUPRI,*) SECNAM,': flag SKIP_PRESCREEN is ',          &
     &                     SKIP_PRESCREEN
            WRITE(LUPRI,*) 'NNSHL is: ',NNSHL
            WRITE(LUPRI,*) 'SIZE(iSP2F) must be equal to NNSHL, ',      &
     &                     'SIZE(iSP2F)= ',SIZE(iSP2F)
            CALL CHO_QUIT('Initialization error in '//SECNAM,102)
         END IF
      ELSE
         CALL CHO_DIASP()
      END IF

!     Get the number of symmetries.
!     -----------------------------

      CALL GET_ISCALAR('nSym',NSYM)  ! Get # irreps.
      IF ((NSYM.LT.1) .OR. (NSYM.GT.8)) THEN
         WRITE(LUPRI,*) 'NSYM out of bounds: ',NSYM
         CALL CHO_QUIT('NSYM out of bounds in '//SECNAM,102)
      END IF

!     NBAS(ISYM): # basis functions (SOs) in symmetry ISYM
!     IBAS(ISYM): offset to basis functions in symmetry ISYM
!     NBAST     : total number of basis functions
!     ------------------------------------------------------

      CALL GET_IARRAY('nBas',NBAS,NSYM)
      IBAS(1) = 0
      NBAST   = NBAS(1)
      DO ISYM = 2,NSYM
         IBAS(ISYM) = NBAST
         NBAST = NBAST + NBAS(ISYM)
      END DO
      IF (NBAST .LT. 1) THEN
         WRITE(LUPRI,*) 'NBAST out of bounds: ',NBAST
         CALL CHO_QUIT('NBAST out of bounds in '//SECNAM,102)
      END IF

!     Allocate shell based index arrays.
!     ----------------------------------

      Call mma_allocate(iBasSh,nSym,nShell,Label='iBasSh')
      Call mma_allocate(nBasSh,nSym,nShell,Label='nBasSh')
      Call mma_allocate(nBstSh,nShell,Label='nBstSh')

!     ISOSHL(I): shell to which SO I belongs
!     --------------------------------------

      Call mma_allocate(iSOShl,NBAST,Label='iSOShl')
      DO ISYM = 1,NSYM
         DO IA = 1,NBAS(ISYM)
            I = IBAS(ISYM) + IA
            iSOShl(I) = ISO2SH(I)
         END DO
      END DO

!     NBASSH(ISYM,ISHL): total dimension of shell ISHL, sym. ISYM
!     NBSTSH(ISHL): total dimension of shell ISHL
!     MXORSH      : max. shell dimension
!     -----------------------------------------------------------

      CALL CHO_SETSH(IBASSH,NBASSH,NBSTSH,                              &
     &               IBAS,NBAS,ISOSHL,NSYM,NSHELL,NBAST)

      MXORSH = NBSTSH(1)
      DO ISHL = 2,NSHELL
         MXORSH = MAX(MXORSH,NBSTSH(ISHL))
      END DO

!     MX2SH: max. dimension of contributing shell pair.
!     -------------------------------------------------

      MX2SH = -1
      DO IJSHL = 1,NNSHL
         IJ = iSP2F(IJSHL)
         CALL CHO_INVPCK(IJ,I,J,.TRUE.)
         IF (I .EQ. J) THEN
            NUMIJ = NBSTSH(I)*(NBSTSH(I) + 1)/2
         ELSE
            NUMIJ = NBSTSH(I)*NBSTSH(J)
         END IF
         MX2SH = MAX(MX2SH,NUMIJ)
      END DO
      IF (MX2SH .LT. 1) THEN
         WRITE(LUPRI,*) 'Max. shell pair dimension non-positive: ',     &
     &                  MX2SH
         CALL CHO_QUIT('Initialization problem in '//SECNAM,102)
      END IF

!     If needed, allocate memory for extracting qualified columns
!     directly in reduced set from Seward.
!     -----------------------------------------------------------

      IF (IFCSEW .EQ. 2) THEN
         Call mma_allocate(iShP2RS,2,Mx2Sh,Label='iShP2RS')
         Call mma_allocate(iShP2Q ,2,Mx2Sh,Label='iShP2Q ')
      END IF

!     ISHLSO(I): index of SO I within its shell
!     -----------------------------------------

      Call mma_allocate(iShlSO,nBasT,Label='iShlSO')
      CALL CHO_SETSH2(iShlSO,iSOShl,NBSTSH,NBAST,NSHELL)

      END
