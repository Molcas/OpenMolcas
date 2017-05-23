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
!SVC: DGA can't distribute stripes for some small dimensions so we have
!to explicitly use the irregular versions. These is a wrapper routine to
!create horizontal (H) and vertical (V) stripes. The dimensions are
!divided evenly, with the remainder spread over the leading stripes.
#ifdef _MOLCAS_MPP_
      SUBROUTINE GA_CREATE_STRIPED (ORI,NROW,NCOL,LABEL,LG_M)
      IMPLICIT NONE
      CHARACTER :: ORI
      CHARACTER(LEN=*) :: LABEL
      INTEGER :: NROW, NCOL, LG_M
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL :: BSTAT
      INTEGER :: NPROCS
      INTEGER :: NBLOCK1, NBLOCK2
      INTEGER, ALLOCATABLE :: MAP1(:), MAP2(:)
      INTEGER :: NDIM, NBASE, NREST, IOFF, I

      NPROCS=GA_NNODES()

      NBLOCK1=1
      ALLOCATE(MAP1(NBLOCK1))
      MAP1(1)=1

      NDIM=0
      IF (ORI.EQ.'H') THEN
        NDIM=NROW
      ELSE IF (ORI.EQ.'V') THEN
        NDIM=NCOL
      END IF
      NBLOCK2=MIN(NDIM,NPROCS)
      NBASE=NDIM/NPROCS
      NREST=MOD(NDIM,NPROCS)
      ALLOCATE(MAP2(NBLOCK2))
      IOFF=1
      DO I=1,NBLOCK2
        MAP2(I)=IOFF
        IF (I.LE.NREST) THEN
          IOFF=IOFF+NBASE+1
        ELSE
          IOFF=IOFF+NBASE
        END IF
      END DO

      BSTAT=.FALSE.
      IF (ORI.EQ.'H') THEN
        BSTAT = GA_CREATE_IRREG (MT_DBL,NROW,NCOL,LABEL,
     &                     MAP2,NBLOCK2,MAP1,NBLOCK1,LG_M)
      ELSE IF (ORI.EQ.'V') THEN
        BSTAT = GA_CREATE_IRREG (MT_DBL,NROW,NCOL,LABEL,
     &                     MAP1,NBLOCK1,MAP2,NBLOCK2,LG_M)
      END IF

      DEALLOCATE(MAP1,MAP2)

      IF (.NOT.bStat) THEN
        WRITE(6,*) 'GA_CREATE_HS: could not create array, abort'
        CALL AbEnd()
      END IF
      END
#elif defined (NAGFOR)
c Some compilers do not like empty files
      SUBROUTINE EMPTY_GA_CREATE_STRIPED ()
      END
#endif
