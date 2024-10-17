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
      SUBROUTINE DIAFCK(NO,FOCK,IOSTA,IOEND,TSCT,NB,CMO1SCT,CMO2SCT)
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer NO,IOSTA,IOEND,NB
      Real*8 FOCK(NO,NO)
      Real*8 TSCT(IOSTA:IOEND,IOSTA:IOEND)
      Real*8 CMO1SCT(NB,*)
      Real*8 CMO2SCT(NB,*)

      Real*8, Allocatable:: TMP(:)

* Use orbital rotations on a subblock of active orbitals
* in order to diagonalize that block of a Fock matrix, and
* in addition apply a counterrotation to the CI array that
* keeps the wave function invariant to the rotation.

* FOCK is the Fock matrix of this symmetry, in square form.
* to be tranformed. TSCT is a corresponding section of the
* transformation matrix, CMO1SCT is a section of the CMO array (in)
* and CMO2SCT the same, transformed.

* Size of orbital section to process:
      NSCT=IOEND+1-IOSTA
* Local scratch area:
      CALL mma_allocate(TMP,NO*NSCT,Label='TMP')
* Put part of FOCK to be diagonalized into triangular scratch area:
      IJ=0
      DO I=IOSTA,IOEND
       DO J=IOSTA,I
        IJ=IJ+1
        TMP(IJ)=FOCK(I,J)
       END DO
      END DO
* Put unit matrix into TSCT:
      CALL DCOPY_(NSCT**2,[0.0D0],0,TSCT,1)
      CALL DCOPY_(NSCT,[1.0D0],0,TSCT,NSCT+1)
C Diagonalize, and order for best submatrix condition:
      CALL Jacob(TMP,TSCT,NSCT,NSCT)
      DO I=IOSTA,IOEND
       JMX=I
       VMX=ABS(TSCT(I,JMX))
       DO J=I+1,IOEND
        V=ABS(TSCT(I,J))
        IF(V.GT.VMX) THEN
          JMX=J
          VMX=V
        END IF
       END DO
       IF(JMX.GT.I) THEN
        DO K=IOSTA,IOEND
         SWAP=TSCT(K,I)
         TSCT(K,I)=TSCT(K,JMX)
         TSCT(K,JMX)=SWAP
        END DO
       END IF
       IF(TSCT(I,I).LT.0.0D0) THEN
        DO K=IOSTA,IOEND
         TSCT(K,I)=-TSCT(K,I)
        END DO
       END IF
      END DO
C Transform the Fock matrix:
      CALL DGEMM_('N','N',NO,NSCT,NSCT,1.0D0,FOCK(1,IOSTA),NO,TSCT,NSCT,
     &              0.0D0,TMP,NO)
      CALL DCOPY_(NO*NSCT,TMP,1,FOCK(1,IOSTA),1)
      DO I=IOSTA,IOEND
       DO J=1,NO
        FOCK(I,J)=TMP(J+NO*(I-IOSTA))
       END DO
      END DO
      CALL DGEMM_('T','N',NSCT,NSCT,NSCT,1.0D0,TSCT,NSCT,
     &         TMP(IOSTA),NO,0.0D0,FOCK(IOSTA,IOSTA),NO)

C Then transform the CMO coeffs:
      CALL DGEMM_('N','N',NB,NSCT,NSCT,1.0D0,CMO1SCT,NB,TSCT,NSCT,
     &              0.0D0,CMO2SCT,NB)

      CALL mma_deallocate(TMP)
      RETURN
      END
