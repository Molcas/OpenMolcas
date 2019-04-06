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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FOCK(NO,NO)
      DIMENSION TSCT(IOSTA:IOEND,IOSTA:IOEND)
      DIMENSION CMO1SCT(NB,*)
      DIMENSION CMO2SCT(NB,*)
#include "WrkSpc.fh"
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
      CALL GETMEM('TMP','ALLO','REAL',LTMP,NO*NSCT)
* Put part of FOCK to be diagonalized into triangular scratch area:
      IJ=0
      DO I=IOSTA,IOEND
       DO J=IOSTA,I
        IJ=IJ+1
        WORK(LTMP-1+IJ)=FOCK(I,J)
       END DO
      END DO
* Put unit matrix into TSCT:
      CALL DCOPY_(NSCT**2,[0.0D0],0,TSCT,1)
      CALL DCOPY_(NSCT,[1.0D0],0,TSCT,NSCT+1)
C Diagonalize, and order for best submatrix condition:
      CALL Jacob(WORK(LTMP),TSCT,NSCT,NSCT)
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
     &              0.0D0,WORK(LTMP),NO)
      CALL DCOPY_(NO*NSCT,WORK(LTMP),1,FOCK(1,IOSTA),1)
      DO I=IOSTA,IOEND
       DO J=1,NO
        FOCK(I,J)=WORK(LTMP-1+J+NO*(I-IOSTA))
       END DO
      END DO
      CALL DGEMM_('T','N',NSCT,NSCT,NSCT,1.0D0,TSCT,NSCT,
     &         WORK(LTMP-1+IOSTA),NO,0.0D0,FOCK(IOSTA,IOSTA),NO)

C Then transform the CMO coeffs:
      CALL DGEMM_('N','N',NB,NSCT,NSCT,1.0D0,CMO1SCT,NB,TSCT,NSCT,
     &              0.0D0,CMO2SCT,NB)

      CALL GETMEM('TMP','FREE','REAL',LTMP,NO*NSCT)
      RETURN
      END
