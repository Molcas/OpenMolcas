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
      SUBROUTINE CHO_RDRSTC(IFAIL)
C
C     Purpose: read decomposition restart info and store in common
C              block. If IFAIL != 0 on exit, some error occurred and,
C              most likely, some of the restart info is not
C              defined/initialized.
C
C     NB!!!! the restart files MUST be open on entry...
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_RDRSTC')

      CHARACTER*6 FRST
      PARAMETER (FRST = 'CHORST')

      PARAMETER (LSCR = 8)
      REAL*8  DSCR(LSCR)
      INTEGER JSCR(LSCR)

      INFRED(I)=IWORK(ip_INFRED-1+I)

C     Set return code.
C     ----------------

      IFAIL = 0

C     Read molecular info.
C     --------------------

      IADR = 0

      IOPT = 2
      NRD  = 4
      CALL IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
      XNSYM   = JSCR(1)
      XNSHELL = JSCR(2)
      XNNSHL  = JSCR(3)
      IF (XNSYM.LT.1 .OR. XNSYM.GT.8) THEN
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': #irreps from restart file: ',XNSYM
         IFAIL = 1
         GO TO 100
      ELSE
         IOPT = 2
         CALL IDAFILE(LURST,IOPT,XNBAS,XNSYM,IADR)
      END IF

C     Read decomposition configuration info.
C     --------------------------------------

      IOPT = 2
      NRD  = 2
      CALL IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
      IF (JSCR(1) .EQ. 0) THEN
         XSCDIAG = .FALSE.
      ELSE IF (JSCR(1) .EQ. 1) THEN
         XSCDIAG = .TRUE.
      ELSE
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': integer flag for screening not recognized:',JSCR(1)
         IFAIL = 2
         GO TO 100
      END IF
      XCHO_ADRVEC = JSCR(2)

      IOPT = 2
      NRD  = 8
      CALL DDAFILE(LURST,IOPT,DSCR,NRD,IADR)
      XTHRCOM  = DSCR(1)
      XTHRDIAG = DSCR(2)
      XDAMP(1) = DSCR(3)
      XDAMP(2) = DSCR(4)
      XSPAN    = DSCR(5)
      XTHRNEG  = DSCR(6)
      XWARNEG  = DSCR(7)
      XTOONEG  = DSCR(8)

C     Read decomposition info.
C     ------------------------

      IOPT = 2
      NRD  = 1
      CALL IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
      XNPASS = JSCR(1)
      IF (XNPASS.LT.1 .OR. XNPASS.GT.MAXRED) THEN
         WRITE(LUPRI,'(A,A,I10)')
     &   SECNAM,': #reduced sets in restart:',XNPASS
         IFAIL = 3
         GO TO 100
      ELSE
         IOPT = 2
         CALL IDAFILE(LURST,IOPT,IWORK(ip_INFRED),XNPASS,IADR)
         IF (INFRED(1) .NE. 0) THEN
            WRITE(LUPRI,'(A,A,I10)')
     &      SECNAM,': disk address of 1st reduced set:',IWORK(ip_INFRED)
            IFAIL = 4
            GO TO 100
         END IF
         LREST = MAXRED - XNPASS
         IF (LREST .GT. 0) CALL CHO_IZERO(IWORK(ip_INFRED+XNPASS),LREST)
      END IF
      DO ISYM = 1,NSYM
         IOPT = 2
         NRD  = 1
         CALL IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
         NUMCHO(ISYM) = JSCR(1)
         IF (NUMCHO(ISYM).LT.0 .OR. NUMCHO(ISYM).GT.MAXVEC) THEN
            WRITE(LUPRI,'(A,A,I2,A,I10)')
     &      SECNAM,': #Cholesky vectors (sym.',ISYM,'): ',NUMCHO(ISYM)
            IFAIL = 5
            GO TO 100
         ELSE IF (NUMCHO(ISYM) .EQ. 0) THEN
            NDIM = MAXVEC*INFVEC_N2
            KOFF = ip_INFVEC + NDIM*(ISYM-1)
            CALL CHO_IZERO(IWORK(KOFF),NDIM)
         ELSE
            DO J = 1,INFVEC_N2
               IOPT = 2
               KOFF = ip_INFVEC + MAXVEC*INFVEC_N2*(ISYM-1)
     &              + MAXVEC*(J-1)
               CALL IDAFILE(LURST,IOPT,IWORK(KOFF),NUMCHO(ISYM),IADR)
               LREST = MAXVEC - NUMCHO(ISYM)
               IF (LREST .GT. 0) THEN
                  KOFF = ip_INFVEC + MAXVEC*INFVEC_N2*(ISYM-1)
     &                 + MAXVEC*(J-1) + NUMCHO(ISYM)
                  CALL CHO_IZERO(IWORK(KOFF),LREST)
               END IF
            END DO
         END IF
      END DO

  100 IF (IFAIL .NE. 0) THEN  ! failures jump to this point
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,': refusing to read more restart info!'
      END IF

      END
