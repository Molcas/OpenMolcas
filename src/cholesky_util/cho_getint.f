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
      SUBROUTINE CHO_GETINT(DIAG,DIASH,ISYSH,LSTQSP,NPOTSH,ICOUNT)
C
C     Purpose: get qualified integral columns for Cholesky decomposition.
C
C     DIASH(ij): max. diagonal in shell pair i,j
C     NPOTSH   : the number of shell pairs that can be qualified.
C
      use ChoArr, only: iSP2F
#include "implicit.fh"
      DIMENSION DIAG(*), DIASH(*)
      INTEGER   ISYSH(*)
      INTEGER   LSTQSP(NPOTSH)
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_GETINT')

      PARAMETER (ZERO = 0.0D0)

      LOGICAL DODECO, FULL, SYNC, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MEMQ(1)

      INTEGER  CHO_ISUMELM
      EXTERNAL CHO_ISUMELM

      INTMAP(I)=IWORK(ip_INTMAP-1+I)

C-tbp: some debugging...
c     IF (LOCDBG) THEN
c        DO LEVEL = 1,3
c           CALL CHO_MCA_INT_1_DBG(DIAG,LEVEL)
c        END DO
c        call cho_quit(SECNAM//' end of test',100)
c     END IF

C     Initializations.
C     ----------------

      CALL CHO_IZERO(NQUAL,NSYM)
      ICOUNT = 0
      IF (MXSHPR .GT. 0) THEN
         MCOUNT = MIN(NPOTSH,MXSHPR)
      ELSE
         MCOUNT = NPOTSH
      END IF
      DODECO = .FALSE.

      MXDIM = NNBSTR(1,2)
      DO ISYM = 2,NSYM
         MXDIM = MAX(MXDIM,NNBSTR(ISYM,2))
      END DO
      CALL CHO_MEM('GetMax','GETM','REAL',KDUM,LMAX)
      XMMQ = DBLE(N1_QUAL)*DBLE(LMAX)/DBLE(N2_QUAL)
      MEMQ(1) = INT(XMMQ)
      CALL CHO_GAIGOP(MEMQ,1,'min')
      IF (MEMQ(1) .LT. MXDIM) THEN
         WRITE(LUPRI,*) SECNAM,': memory split error!'
         WRITE(LUPRI,*) 'Memory for storing qualified columns: ',MEMQ(1)
         WRITE(LUPRI,*) 'Minimal memory needed to store one column: ',
     &                  MXDIM
         WRITE(LUPRI,*) 'Total memory available: ',LMAX
         WRITE(LUPRI,*) 'Memory split is ',N1_QUAL,'/',N2_QUAL,
     &                  ' for qualified columns.'
         WRITE(LUPRI,*) 'Change memory split in input file...'
         CALL CHO_QUIT('Memory split error in '//SECNAM,101)
      END IF

C     Shell pair qualification loop.
C     ------------------------------

      DO WHILE ((.NOT.DODECO) .AND. (ICOUNT.LT.MCOUNT))

C        Update shell pair counter.
C        --------------------------

         ICOUNT = ICOUNT + 1

C        Get shell pair corresponding to largest diagonal.
C        -------------------------------------------------

         CALL CHO_P_GETMAXSHL(DIASH,SMAX,ISHLAB)
         CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
         ISYMAB = ISYSH(ISHLAB)

         IF ((SMAX.EQ.ZERO) .OR. (ABS(SMAX).LT.DIAMIN(ISYMAB))) THEN

C           Diagonal too small to be qualified for decomposition.
C           -----------------------------------------------------

            IF (ICOUNT .EQ. 1) THEN
               WRITE(LUPRI,*) SECNAM,': no integrals calculated; ',
     &                        'unable to proceed to decomposition!'
               WRITE(LUPRI,*) 'Max. abs. diagonal for shell pair ',
     &                        ISHLA,', ',ISHLB,': ',ABS(SMAX)
               WRITE(LUPRI,*) 'Max. abs. diagonal allowed: ',
     &                        DIAMIN(ISYMAB),' (sym. ',ISYMAB,')'
               CALL CHO_QUIT('Severe error in '//SECNAM,104)
            ELSE
               ICOUNT = ICOUNT - 1
               NSEL   = CHO_ISUMELM(NQUAL,NSYM)
               DODECO = NSEL .GT. 0
            END IF

         ELSE

C           Qualify diagonals within this shell pair.
C           -----------------------------------------

            SYNC = .FALSE.
            FULL = .FALSE.
            CALL CHO_P_QUALIFY(DIAG,SYNC,ISHLAB,ISYMAB,MEMQ(1),FULL)

C           Calculate integral columns; get qualified ones stored in
C           current reduced set; write these to disk on temporary
C           file(s).
C           --------------------------------------------------------

            NSEL   = CHO_ISUMELM(NQUAL,NSYM)
            NCOLAB = NSEL - CHO_ISUMELM(IOFFQ,NSYM)

            IF (NCOLAB .GT. 0) THEN

               IWORK(ip_INTMAP-1+ISHLAB) = INTMAP(ISHLAB) + 1
               IF (IPRINT .GE. INF_IN2) THEN
                  WRITE(LUPRI,'(/,A,I5,1X,I5,A,I9,A)')
     &            'Calculating shell pair (**|',ISHLA,ISHLB,
     &            '):',NCOLAB,' columns have been qualified'
                  WRITE(LUPRI,'(80A)') ('=',i=1,77)
                  WRITE(LUPRI,'(A,I12)')
     &            'Number of calculations so far for this shell pair: ',
     &            INTMAP(ISHLAB)
               END IF

               LSTQSP(ICOUNT) = ISHLAB
               CALL CHO_MCA_CALCINT(ISHLAB)

C              Enough integral columns for proceeding to decomposition?
C              --------------------------------------------------------

               DODECO = FULL .OR. NSEL.GE.MINQUAL

            ELSE IF (NCOLAB .EQ. 0) THEN

               IF (NSEL .LT. 1) THEN
                  WRITE(LUPRI,*) SECNAM,': logical error: ',
     &                                  'unable to qualify diagonals'
                  WRITE(LUPRI,*) SECNAM,': NCOLAB = ',NCOLAB
                  WRITE(LUPRI,*) SECNAM,': NSEL   = ',NSEL
                  CALL CHO_QUIT('[0] Logical error in '//SECNAM,104)
               ELSE
                  ICOUNT = ICOUNT - 1
                  DODECO = .TRUE.
               END IF

            ELSE

               WRITE(LUPRI,*) SECNAM,': logical error: ',
     &                               'unable to qualify diagonals'
               WRITE(LUPRI,*) SECNAM,': NCOLAB = ',NCOLAB
               WRITE(LUPRI,*) SECNAM,': NSEL   = ',NSEL
               CALL CHO_QUIT('[1] Logical error in '//SECNAM,104)

            END IF

         END IF

      END DO

C     Test loop exit (we may have calculated all possible integrals, yet
C     NSEL < MINQUAL or allowed memory may have been used).
C     ------------------------------------------------------------------

      IF (.NOT. DODECO) THEN
         NSEL = CHO_ISUMELM(NQUAL,NSYM)
         IF (NSEL .LT. 1) THEN
            WRITE(LUPRI,*) SECNAM,': logical error: ',
     &                            'unable to qualify diagonals'
            WRITE(LUPRI,*) SECNAM,': Flag DODECO is ',DODECO
            WRITE(LUPRI,*) SECNAM,': NSEL    = ',NSEL
            WRITE(LUPRI,*) SECNAM,': ICOUNT  = ',ICOUNT
            WRITE(LUPRI,*) SECNAM,': MCOUNT  = ',MCOUNT
            WRITE(LUPRI,*) SECNAM,': NPOTSH  = ',NPOTSH
            WRITE(LUPRI,*) SECNAM,': MINQUAL = ',MINQUAL
            CALL CHO_QUIT('[2] Logical error in '//SECNAM,103)
         ELSE
            DODECO = .TRUE.
         END IF
      END IF

C     Set indices for local qualified (parallel runs).
C     ------------------------------------------------

      CALL CHO_P_SETLQ()

      END
