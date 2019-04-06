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
      SUBROUTINE TRANSFOCK(TORB,F,IDIR)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      DIMENSION TORB(*),F(*)
* Purpose: given an orbital transformation array and some
* one-electron matrix in storage format as e.g. HONE, FIFA,
* transform the matrix to use the new orbital basis.
      CALL QENTER('TRANSFOCK')

      NT=0
      NOMX=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NI+NR1+NR2+NR3+NS
        NOMX=MAX(NOMX,NO)
        NT=NT+NI**2+NR1**2+NR2**2+NR3**2+NS**2
      END DO
      CALL GETMEM('FSQ','ALLO','REAL',LFSQ,NOMX**2)
      CALL GETMEM('TSQ','ALLO','REAL',LTSQ,NOMX**2)
      CALL GETMEM('TMP','ALLO','REAL',LTMP,NOMX**2)
      IJOFF=0
      ITOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NI+NR1+NR2+NR3+NS
        IF (NO.eq.0) GOTO 99
* Copy the matrices to square storage: first fill with zeroes.
        CALL DCOPY_(NO**2,[0.0D0],0,WORK(LTSQ),1)
* Copy inactive TORB block to TSQ
        IOFF=0
        DO I=1,NI
         DO J=1,NI
          WORK(LTSQ-1+I+NO*(J-1))=TORB(ITOFF+I+NI*(J-1))
         END DO
        END DO
* Copy ras1 T block to TSQ, and so on..
        ITOFF=ITOFF+NI**2
        IOFF=IOFF+NI
        DO I=1,NR1
         II=IOFF+I
         DO J=1,NR1
          JJ=IOFF+J
          WORK(LTSQ-1+II+NO*(JJ-1))=TORB(ITOFF+I+NR1*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR1**2
        IOFF=IOFF+NR1
        DO I=1,NR2
         II=IOFF+I
         DO J=1,NR2
          JJ=IOFF+J
          WORK(LTSQ-1+II+NO*(JJ-1))=TORB(ITOFF+I+NR2*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR2**2
        IOFF=IOFF+NR2
        DO I=1,NR3
         II=IOFF+I
         DO J=1,NR3
          JJ=IOFF+J
          WORK(LTSQ-1+II+NO*(JJ-1))=TORB(ITOFF+I+NR3*(J-1))
         END DO
        END DO
*--- Finally, the secondary orbitals (non-deleted, virtual).
        ITOFF=ITOFF+NR3**2
        IOFF=IOFF+NR3
        DO I=1,NS
         II=IOFF+I
         DO J=1,NS
          JJ=IOFF+J
          WORK(LTSQ-1+II+NO*(JJ-1))=TORB(ITOFF+I+NS*(J-1))
         END DO
        END DO
        ITOFF=ITOFF+NS**2
* Now transfer the Fock matrix block to square storage:
        IJ=0
        DO I=1,NO
         DO J=1,I
          IJ=IJ+1
          WORK(LFSQ-1+J+NO*(I-1))=F(IJOFF+IJ)
          WORK(LFSQ-1+I+NO*(J-1))=F(IJOFF+IJ)
         END DO
        END DO
       IF (IDIR.GE.0) THEN
* Transform, first do FSQ*TSQ -> TMP...
        CALL DGEMM_('N','N',NO,NO,NO,1.0D0,WORK(LFSQ),NO,WORK(LTSQ),NO,
     &              0.0D0,WORK(LTMP),NO)
* ... and then do TSQ(transpose)*TMP -> FSQ.
        CALL DGEMM_('T','N',NO,NO,NO,1.0D0,WORK(LTSQ),NO,WORK(LTMP),NO,
     &              0.0D0,WORK(LFSQ),NO)
       ELSE
* Or inverse transformation
        CALL DGEMM_('N','T',NO,NO,NO,1.0D0,WORK(LFSQ),NO,WORK(LTSQ),NO,
     &              0.0D0,WORK(LTMP),NO)
        CALL DGEMM_('N','N',NO,NO,NO,1.0D0,WORK(LTSQ),NO,WORK(LTMP),NO,
     &              0.0D0,WORK(LFSQ),NO)
       END IF
* Transfer FSQ values back to F, in triangular storage.
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         F(IJOFF+IJ)=WORK(LFSQ-1+I+NO*(J-1))
        END DO
       END DO
       IJOFF=IJOFF+(NO*(NO+1))/2
* and repeat, using next symmetry block.
  99   CONTINUE
      END DO
      CALL GETMEM('FSQ','FREE','REAL',LFSQ,NOMX**2)
      CALL GETMEM('TSQ','FREE','REAL',LTSQ,NOMX**2)
      CALL GETMEM('TMP','FREE','REAL',LTMP,NOMX**2)

      CALL QEXIT('TRANSFOCK')

      RETURN
      END
