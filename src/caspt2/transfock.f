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
      SUBROUTINE TRANSFOCK(TORB,NTORB,F,NF,IDIR)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nOMx, nIsh, nRas1, nRas2, nRas3, nSsh,
     &                         nSym
      IMPLICIT None
      integer(kind=iwp), intent(in):: NTORB, NF, IDIR
      real(kind=wp), intent(in)::  TORB(NTORB)
      real(kind=wp), intent(inout):: F(NF)

      real(kind=wp), ALLOCATABLE:: FSQ(:), TSQ(:), TMP(:)
      integer(kind=iwp) NT, NI, NR1, NR2, NR3, NS, NO, IJOFF,
     &                  ITOFF, I, J, II, JJ, IJ, ISYM, IOFF
* Purpose: given an orbital transformation array and some
* one-electron matrix in storage format as e.g. HONE, FIFA,
* transform the matrix to use the new orbital basis.

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
      CALL mma_allocate(FSQ,NOMX**2,LABEL='FSQ')
      CALL mma_allocate(TSQ,NOMX**2,LABEL='TSQ')
      CALL mma_allocate(TMP,NOMX**2,LABEL='TMP')
      IJOFF=0
      ITOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NI+NR1+NR2+NR3+NS
        IF (NO.eq.0) Cycle
* Copy the matrices to square storage: first fill with zeroes.
        TSQ(1:NO**2)=Zero
* Copy inactive TORB block to TSQ
        IOFF=0
        DO I=1,NI
         DO J=1,NI
          TSQ(I+NO*(J-1))=TORB(ITOFF+I+NI*(J-1))
         END DO
        END DO
* Copy ras1 T block to TSQ, and so on..
        ITOFF=ITOFF+NI**2
        IOFF=IOFF+NI
        DO I=1,NR1
         II=IOFF+I
         DO J=1,NR1
          JJ=IOFF+J
          TSQ(II+NO*(JJ-1))=TORB(ITOFF+I+NR1*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR1**2
        IOFF=IOFF+NR1
        DO I=1,NR2
         II=IOFF+I
         DO J=1,NR2
          JJ=IOFF+J
          TSQ(II+NO*(JJ-1))=TORB(ITOFF+I+NR2*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR2**2
        IOFF=IOFF+NR2
        DO I=1,NR3
         II=IOFF+I
         DO J=1,NR3
          JJ=IOFF+J
          TSQ(II+NO*(JJ-1))=TORB(ITOFF+I+NR3*(J-1))
         END DO
        END DO
*--- Finally, the secondary orbitals (non-deleted, virtual).
        ITOFF=ITOFF+NR3**2
        IOFF=IOFF+NR3
        DO I=1,NS
         II=IOFF+I
         DO J=1,NS
          JJ=IOFF+J
          TSQ(II+NO*(JJ-1))=TORB(ITOFF+I+NS*(J-1))
         END DO
        END DO
        ITOFF=ITOFF+NS**2
* Now transfer the Fock matrix block to square storage:
        IJ=0
        DO I=1,NO
         DO J=1,I
          IJ=IJ+1
          FSQ(J+NO*(I-1))=F(IJOFF+IJ)
          FSQ(I+NO*(J-1))=F(IJOFF+IJ)
         END DO
        END DO
       IF (IDIR.GE.0) THEN
* Transform, first do FSQ*TSQ -> TMP...
        CALL DGEMM_('N','N',NO,NO,NO,One,FSQ,NO,TSQ,NO,
     &              Zero,TMP,NO)
* ... and then do TSQ(transpose)*TMP -> FSQ.
        CALL DGEMM_('T','N',NO,NO,NO,One,TSQ,NO,TMP,NO,
     &              Zero,FSQ,NO)
       ELSE
* Or inverse transformation
        CALL DGEMM_('N','T',NO,NO,NO,One,FSQ,NO,TSQ,NO,
     &              Zero,TMP,NO)
        CALL DGEMM_('N','N',NO,NO,NO,One,TSQ,NO,TMP,NO,
     &              Zero,FSQ,NO)
       END IF
* Transfer FSQ values back to F, in triangular storage.
       IJ=0
       DO I=1,NO
        DO J=1,I
         IJ=IJ+1
         F(IJOFF+IJ)=FSQ(I+NO*(J-1))
        END DO
       END DO
       IJOFF=IJOFF+(NO*(NO+1))/2
* and repeat, using next symmetry block.
      END DO
      CALL mma_deallocate(FSQ)
      CALL mma_deallocate(TSQ)
      CALL mma_deallocate(TMP)

      END SUBROUTINE TRANSFOCK
