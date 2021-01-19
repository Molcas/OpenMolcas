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
      SUBROUTINE CHO_QUALIFY(DIAG,ISHLAB,ISYMAX,MEM,FULL)
C
C     Purpose: qualify diagonal elements for decomposition in
C              current reduced set. ISYMAX is the symmetry block
C              to which the largest diagonal belongs.
C              MEM is the (total!) max. allowed
C              memory for storing the qualified columns.
C              If no more columns can be qualified on exit,
C              FULL=.true. is returned.
C
      use ChoArr, only: iSP2F
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRed
#include "implicit.fh"
      DIMENSION DIAG(*)
      LOGICAL   FULL
#include "cholesky.fh"
#include "choptr.fh"

      INTEGER  CHO_IDOT
      EXTERNAL CHO_IDOT

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_QUALIFY')

      LOGICAL LOCDBG
#if defined (_DEBUGPRINT_)
      PARAMETER (LOCDBG = .TRUE.)
#else
      PARAMETER (LOCDBG = .FALSE.)
#endif

C     Copy counter to offset array.
C     -----------------------------

      CALL ICOPY(NSYM,NQUAL,1,IOFFQ,1)

C     Check memory.
C     -------------

      MEM0 = CHO_IDOT(NSYM,NQUAL,1,NNBSTR(1,2),1)
      LEFT = MEM - MEM0
      IF (IALQUA .EQ. 0) THEN
         MINM = NNBSTR(1,2)
         DO ISYM = 1,NSYM
            MINM = MAX(MINM,NNBSTR(ISYM,2))
         END DO
      ELSE
         MINM = NNBSTR(ISYMAX,2)
      END IF
      FULL = LEFT .LT. MINM
      IF (FULL) RETURN

C     Qualify.
C     --------

      IF (IALQUA .EQ. 0) THEN  ! qualify until full (dalton style)
         DO ISYM = 1,NSYM
            CALL CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
         END DO
      ELSE IF (IALQUA .EQ. 1) THEN  ! qualify until full
         CALL CHO_QUALIFY_1(DIAG,ISYMAX,ISHLAB,MEM,MEM0,LEFT)
         DO ISYM = 1,ISYMAX-1
            CALL CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
         END DO
         DO ISYM = ISYMAX+1,NSYM
            CALL CHO_QUALIFY_1(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
         END DO
      ELSE ! qualify until full, then largest
         CALL CHO_QUALIFY_2(DIAG,ISYMAX,ISHLAB,MEM,MEM0,LEFT)
         DO ISYM = 1,ISYMAX-1
            CALL CHO_QUALIFY_2(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
         END DO
         DO ISYM = ISYMAX+1,NSYM
            CALL CHO_QUALIFY_2(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
         END DO
      END IF

C     Set FULL flag:
C     FULL=.true. if a) not enough memory to qualify another column of
C     any symmetry, or b) MAXQUAL reached in any symmetry.
C     ----------------------------------------------------------------

      NEED = CHO_IDOT(NSYM,NQUAL,1,NNBSTR(1,2),1)
      IF (NEED.LT.1 .OR. NEED.GT.MEM) THEN
         CALL CHO_QUIT('Logical error (2) in '//SECNAM,104)
      ELSE
         LEFT = MEM - NEED
         FULL = .FALSE.
         ISYM = 0
         DO WHILE (ISYM.LT.NSYM .AND. .NOT.FULL)
            ISYM = ISYM + 1
            IF (NQUAL(ISYM).LT.IOFFQ(ISYM) .OR. NQUAL(ISYM).LT.0 .OR.
     &          NQUAL(ISYM).GT.MAXQUAL) THEN
               CALL CHO_QUIT('Logical error (3) in '//SECNAM,104)
            ELSE
               FULL = NQUAL(ISYM) .EQ. MAXQUAL
            END IF
            IF (NNBSTR(ISYM,2) .GT. 0) THEN
               FULL = FULL .OR. LEFT.LT.NNBSTR(ISYM,2)
            END IF
         END DO
      END IF

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': qualified diagonals from shell-pair ',
     &                  ISHLA,ISHLB,':'
         WRITE(LUPRI,*) 'Qualification algorithm: ',IALQUA
         WRITE(LUPRI,*) 'Total memory for qualification: ',MEM,
     &                  '  Memory left: ',LEFT
         DO ISYM = 1,NSYM
            NUM = NQUAL(ISYM) - IOFFQ(ISYM)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Sym.,dimension,#qualified,threshold: ',
     &                     ISYM,NNBSTRSH(ISYM,ISHLAB,2),NUM,DIAMIN(ISYM)
            IF (NNBSTRSH(ISYM,ISHLAB,2) .GT. 0) THEN
               I1 = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2) + 1
               I2 = I1 + NNBSTRSH(ISYM,ISHLAB,2) - 1
               WRITE(LUPRI,*) 'Diagonal (current reduced set):'
               WRITE(LUPRI,'(5F15.8)') (DIAG(INDRED(I,2)), I=I1,I2)
               K1 = IOFFQ(ISYM) + 1
               K2 = NQUAL(ISYM)
               WRITE(LUPRI,*) 'Qualified diagonals:'
               WRITE(LUPRI,'(5F15.8)') (DIAG(INDRED(IQUAB(K,ISYM),2)),
     &                                  K = K1,K2)
            END IF
         END DO
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
      END IF

      END
