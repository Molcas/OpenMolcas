************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
      SUBROUTINE ZTRNSF_MASKED(N,UR,UI,AR,AI,IJ,IST,INUM,JST,JNUM)
      IMPLICIT NONE
      INTEGER :: N,INUM,JNUM
      REAL*8 :: UR(N,N),UI(N,N)
      REAL*8 :: AR(N,N),AI(N,N)
      INTEGER :: IJ(4),IST(INUM),JST(JNUM)
#include "real.fh"
#include "stdalloc.fh"
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MR,MI,VR,VI,TR,TI
      INTEGER :: I,J,II,JJ,NI,NJ

      NI=IJ(2)-IJ(1)+1
      NJ=IJ(4)-IJ(3)+1

      CALL mma_allocate(MR,INUM,JNUM,LABEL='MR')
      CALL mma_allocate(MI,INUM,JNUM,LABEL='MI')
      DO J=1,JNUM
        DO I=1,INUM
          MR(I,J)=AR(IST(I),JST(J))
          MI(I,J)=AI(IST(I),JST(J))
        END DO
      END DO

      CALL mma_allocate(VR,JNUM,NJ,LABEL='VR')
      CALL mma_allocate(VI,JNUM,NJ,LABEL='VI')
      DO J=1,NJ
        JJ=IJ(3)+J-1
        DO I=1,JNUM
          VR(I,J)=UR(JST(I),JJ)
          VI(I,J)=UI(JST(I),JJ)
        END DO
      END DO

      CALL mma_allocate(TR,INUM,NJ,LABEL='TR')
      CALL mma_allocate(TI,INUM,NJ,LABEL='TI')
      CALL DGEMM_('N','N',INUM,NJ,JNUM, One,MR,INUM,VR,JNUM,
     &                                 Zero,TR,INUM)
      CALL DGEMM_('N','N',INUM,NJ,JNUM,-One,MI,INUM,VI,JNUM,
     &                                  One,TR,INUM)
      CALL DGEMM_('N','N',INUM,NJ,JNUM, One,MR,INUM,VI,JNUM,
     &                                 Zero,TI,INUM)
      CALL DGEMM_('N','N',INUM,NJ,JNUM, One,MI,INUM,VR,JNUM,
     &                                  One,TI,INUM)

      CALL mma_deallocate(VR)
      CALL mma_deallocate(VI)
      CALL mma_allocate(VR,INUM,NI,LABEL='VR')
      CALL mma_allocate(VI,INUM,NI,LABEL='VI')
      DO J=1,NI
        JJ=IJ(1)+J-1
        DO I=1,INUM
          VR(I,J)=UR(IST(I),JJ)
          VI(I,J)=UI(IST(I),JJ)
        END DO
      END DO

      CALL mma_deallocate(MR)
      CALL mma_deallocate(MI)
      CALL mma_allocate(MR,NI,NJ,LABEL='MR')
      CALL mma_allocate(MI,NI,NJ,LABEL='MI')
      CALL DGEMM_('T','N',NI,NJ,INUM, One,VR,INUM,TR,INUM,Zero,MR,NI)
      CALL DGEMM_('T','N',NI,NJ,INUM, One,VI,INUM,TI,INUM, One,MR,NI)
      CALL DGEMM_('T','N',NI,NJ,INUM, One,VR,INUM,TI,INUM,Zero,MI,NI)
      CALL DGEMM_('T','N',NI,NJ,INUM,-One,VI,INUM,TR,INUM, One,MI,NI)

      CALL DCOPY_(N*N,[Zero],0,AR,1)
      CALL DCOPY_(N*N,[Zero],0,AI,1)
      DO J=1,NJ
        JJ=IJ(3)+J-1
        DO I=1,NI
          II=IJ(1)+I-1
          AR(II,JJ)=MR(I,J)
          AI(II,JJ)=MI(I,J)
        END DO
      END DO

      CALL mma_deallocate(TR)
      CALL mma_deallocate(TI)
      CALL mma_deallocate(VR)
      CALL mma_deallocate(VI)
      CALL mma_deallocate(MR)
      CALL mma_deallocate(MI)

      RETURN
      END
