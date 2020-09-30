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
      SUBROUTINE TRANSDREF(TORB,DREF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      DIMENSION TORB(*),DREF(NDREF)
* Purpose: given an orbital transformation array
* transform the DREF array (blocked triangular, active)

      NT=0
      NAMX=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NI+NR1+NR2+NR3+NS
        NA=NR1+NR2+NR3
        NAMX=MAX(NAMX,NA)
        NT=NT+NI**2+NR1**2+NR2**2+NR3**2+NS**2
      END DO
      CALL GETMEM('DSQ','ALLO','REAL',LDSQ,NAMX**2)
      CALL GETMEM('TSQ','ALLO','REAL',LTSQ,NAMX**2)
      CALL GETMEM('TMP','ALLO','REAL',LTMP,NAMX**2)
      IDOFF=0
      ITOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NA=NR1+NR2+NR3
        NS=NSSH(ISYM)
        NO=NI+NA+NS
        IF (NO.eq.0) GOTO 99
* Copy the matrices to square storage: first fill with zeroes.
        CALL DCOPY_(NA**2,[0.0D0],0,WORK(LTSQ),1)
        IOFF=0
        ITOFF=ITOFF+NI**2
        DO I=1,NR1
         II=IOFF+I
         DO J=1,NR1
          JJ=IOFF+J
          WORK(LTSQ-1+II+NA*(JJ-1))=TORB(ITOFF+I+NR1*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR1**2
        IOFF=IOFF+NR1
        DO I=1,NR2
         II=IOFF+I
         DO J=1,NR2
          JJ=IOFF+J
          WORK(LTSQ-1+II+NA*(JJ-1))=TORB(ITOFF+I+NR2*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR2**2
        IOFF=IOFF+NR2
        DO I=1,NR3
         II=IOFF+I
         DO J=1,NR3
          JJ=IOFF+J
          WORK(LTSQ-1+II+NA*(JJ-1))=TORB(ITOFF+I+NR3*(J-1))
         END DO
        END DO
*--- Finally, the secondary orbitals (non-deleted, virtual).
        ITOFF=ITOFF+NR3**2
        ITOFF=ITOFF+NS**2
* Now transfer the DREF matrix block to square storage:
        IJ=0
        DO I=1,NA
         DO J=1,I
          IJ=IJ+1
          WORK(LDSQ-1+J+NA*(I-1))=DREF(IDOFF+IJ)
          WORK(LDSQ-1+I+NA*(J-1))=DREF(IDOFF+IJ)
         END DO
        END DO
* Transform, first do DSQ*TSQ -> TMP...
       CALL DGEMM_('N','N',NA,NA,NA,1.0D0,WORK(LDSQ),NA,WORK(LTSQ),NA,
     &              0.0D0,WORK(LTMP),NA)
* ... and then do TSQ(transpose)*TMP -> DSQ...
       CALL DGEMM_('T','N',NA,NA,NA,1.0D0,WORK(LTSQ),NA,WORK(LTMP),NA,
     &              0.0D0,WORK(LDSQ),NA)
* Transfer DSQ values back to D, in triangular storage.
       IJ=0
       DO I=1,NA
        DO J=1,I
         IJ=IJ+1
         DREF(IDOFF+IJ)=WORK(LDSQ-1+I+NA*(J-1))
        END DO
       END DO
       IDOFF=IDOFF+(NA*(NA+1))/2
* and repeat, using next symmetry block.
  99   CONTINUE
      END DO
      CALL GETMEM('DSQ','FREE','REAL',LDSQ,NAMX**2)
      CALL GETMEM('TSQ','FREE','REAL',LTSQ,NAMX**2)
      CALL GETMEM('TMP','FREE','REAL',LTMP,NAMX**2)


      RETURN
      END
