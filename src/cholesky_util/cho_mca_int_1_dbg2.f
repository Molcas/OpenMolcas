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
      SUBROUTINE CHO_MCA_INT_1_DBG2()
C
C     Purpose: test symmetry of integral matrix.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*18 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_INT_1_DBG2')

      LOGICAL PRTINT
      PARAMETER (PRTINT = .FALSE.)

      PARAMETER (THR = 1.0D-14)

      INTEGER  CHO_ISAOSH
      EXTERNAL CHO_ISAOSH

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
      NBSTSH(I)=IWORK(ip_NBSTSH-1+I)
      ISP2F(I)=IWORK(ip_iSP2F-1+I)

      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) SECNAM,': testing integral matrix symmetry'
      WRITE(LUPRI,*)

C     Force computation of full shell quadruple.
C     ------------------------------------------

      IF (IFCSEW .NE. 1) THEN
         WRITE(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',
     &                  IFCSEW,' to 1.'
         IFCSEW = 1
      END IF

      LINTT = 2*MX2SH*MX2SH
      CALL GETMEM('Int1.dbg2.1','ALLO','REAL',KINT1,LINTT)
      CALL GETMEM('Int1.dbg2.2','MAX ','REAL',KSEW,LSEW)
      CALL XSETMEM_INTS(LSEW)

      NTST = 0
      NERR = 0
      DO ISHLAB = 1,NNSHL

         CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
         IF (ISHLB .EQ. ISHLA) THEN
            NUMAB = NBSTSH(ISHLA)*(NBSTSH(ISHLB) + 1)/2
         ELSE
            NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
         END IF

         DO ISHLCD = 1,ISHLAB

            CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
            IF (ISHLD .EQ. ISHLC) THEN
               NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLD) + 1)/2
            ELSE
               NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
            END IF
            LINT  = NUMCD*NUMAB
            KINT2 = KINT1 + LINT

C           Calculate integrals (CD|AB).
C           ----------------------------

            CALL CHO_DZERO(WORK(KINT1),LINT)
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,WORK(KINT1),LINT,
     &                         PRTINT)

C           Calculate integrals (AB|CD).
C           ----------------------------

            CALL CHO_DZERO(WORK(KINT2),LINT)
            CALL CHO_MCA_INT_1(ISHLAB,ISHLCD,WORK(KINT2),LINT,
     &                         PRTINT)

C           Compare.
C           --------

            ITST = 0
            IERR = 0
            CALL CHO_MCA_INT1_1_DBG2_CMP(WORK(KINT1),WORK(KINT2),
     &                                  NUMCD,NUMAB,ERRMIN,ICDMN,IABMN,
     &                                  ERRMAX,ICDMX,IABMX,ITST,IERR,
     &                                  THR,.FALSE.,LUPRI)
            NTST = NTST + ITST
            NERR = NERR + IERR

            WRITE(LUPRI,*) '#sym. errors for ',
     &      '(',ISHLC,ISHLD,'|',ISHLA,ISHLB,'): ',IERR,
     &      ' #tests: ',ITST
            IF (IERR .NE. 0) THEN
c              WRITE(LUPRI,*) '    Here is the shell quadruple in INT1:'
c              CALL CHO_OUTPUT(WORK(KINT1),1,NUMCD,1,NUMAB,NUMCD,NUMAB,
c    &                         1,LUPRI)
c              WRITE(LUPRI,*) '    And the shell quadruple in INT2:'
c              CALL CHO_OUTPUT(WORK(KINT1),1,NUMAB,1,NUMCD,NUMAB,NUMCD,
c    &                         1,LUPRI)
               DO IB = 1,NBSTSH(ISHLB)
                  ISYMB = CHO_ISAOSH(IB,ISHLB)
                  DO IA = 1,NBSTSH(ISHLA)
                     ISYMA  = CHO_ISAOSH(IA,ISHLA)
                     ISYMAB = MULD2H(ISYMA,ISYMB)
                     IF (ISHLB .EQ. ISHLA) THEN
                        IAB = ITRI(IA,IB)
                     ELSE
                        IAB = NBSTSH(ISHLA)*(IB - 1) + IA
                     END IF
                     DO ID = 1,NBSTSH(ISHLD)
                        ISYMD = CHO_ISAOSH(ID,ISHLD)
                        DO IC = 1,NBSTSH(ISHLC)
                           ISYMC = CHO_ISAOSH(IC,ISHLC)
                           ISYMCD = MULD2H(ISYMC,ISYMD)
                           IF (ISHLC .EQ. ISHLD) THEN
                              ICD = ITRI(IC,ID)
                           ELSE
                              ICD = NBSTSH(ISHLC)*(ID - 1) + IC
                           END IF
                           IABCD = NUMCD*(IAB - 1) + ICD
                           TST   = ABS(WORK(KINT1+IABCD-1))
                           IF ((TST.GT.0.0D0).AND.(ISYMCD.NE.ISYMAB))
     &                     THEN
                              WRITE(LUPRI,*) 'Symmetry break!!'
                              WRITE(LUPRI,*) 'element ',ICD,IAB,
     &                                       ' is non-zero: ',
     &                                       WORK(KINT1+IABCD-1)
                              WRITE(LUPRI,*) 'Symmetry is: ',
     &                                       MULD2H(ISYMCD,ISYMAB)
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END IF

         END DO
      END DO

      CALL XRLSMEM_INTS
      CALL GETMEM('Int1.flsh','FLUSH','REAL',KINT1,LINTT)
      CALL GETMEM('Int1.free','FREE','REAL',KINT1,LINTT)

      WRITE(LUPRI,*) '***END OF ',SECNAM,': #tests: ',NTST,
     &               ' #errors: ',NERR

      END
