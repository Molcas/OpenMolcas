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
      SUBROUTINE MLTR1 (IMLTOP,LST1,X,F,Y)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: MyRank, nProcs, Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),F(*),Y(*)
      DIMENSION LST1(4,NLST1)
#include "sigma.fh"

C Given a lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
C four entries called L1,L2,L3,L4 for short, for a given
C item, and with V=VAL1(L4),
C compute the Rank-1 updates, if IMLTOP=0,
C X(L1,p,q):= Add V*F(L2,p)*Y(L3,q), p=1..LEN1, q=1..LEN2
C else, the conjugate expression
C Y(L3,q):= Add V*F(L2,p)*X(L1,p,q), p=1..LEN1, q=1..LEN2
C or for IMLTOP=2, compute
C     F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
C However, strides etc can vary from case to case: The indices
C may appear in 1st, 2nd or 3rd position, Fortran-address-wise,
C thus an independent stride is to be given for each of the
C indices in the above 'formal' index ordering, as follows:
C The formal Y(p,q,r) is accessed as
C Y(1+INCX1*(p-1)+INCX2*(q-1)+INCX3*(r-1)), etc.

CSVC: determine outer loop properties
#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        ILST1_IOFF=MYRANK+1
        ILST1_SKIP=NPROCS
      ELSE
        ILST1_IOFF=1
        ILST1_SKIP=1
      ENDIF
#else
      ILST1_IOFF=1
      ILST1_SKIP=1
#endif

      IF(IMLTOP.EQ.0) THEN
        DO ILST=ILST1_IOFF,NLST1,ILST1_SKIP
          L1=LST1(1,ILST)
          L2=LST1(2,ILST)
          L3=LST1(3,ILST)
          L4=LST1(4,ILST)
          V=VAL1(L4)
          IX=INCX1*(L1-1)+1
          IF=INCF1*(L2-1)+1
          IY=INCY1*(L3-1)+1
          DO I=1,LEN1
            A=V*F(IF)
            CALL DAXPY_(LEN2,A,Y(IY),INCY2,X(IX),INCX3)
            IX=IX+INCX2
            IF=IF+INCF2
          END DO
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        DO ILST=ILST1_IOFF,NLST1,ILST1_SKIP
          L1=LST1(1,ILST)
          L2=LST1(2,ILST)
          L3=LST1(3,ILST)
          L4=LST1(4,ILST)
          V=VAL1(L4)
          IX=INCX1*(L1-1)+1
          IF=INCF1*(L2-1)+1
          IY=INCY1*(L3-1)+1
          DO I=1,LEN2
            Y(IY)=Y(IY)+V*DDOT_(LEN1,F(IF),INCF2,X(IX),INCX2)
            IX=IX+INCX3
            IY=IY+INCY2
          END DO
        END DO
      ELSE
        DO ILST=ILST1_IOFF,NLST1,ILST1_SKIP
          L1=LST1(1,ILST)
          L2=LST1(2,ILST)
          L3=LST1(3,ILST)
          L4=LST1(4,ILST)
          V=VAL1(L4)
          IX=INCX1*(L1-1)+1
          IF=INCF1*(L2-1)+1
          IY=INCY1*(L3-1)+1
C     F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
          DO I=1,LEN2
            A=V*Y(IY)
            CALL DAXPY_(LEN1,A,X(IX),INCX2,F(IF),INCF2)
            IX=IX+INCX3
            IY=IY+INCY2
          END DO
        END DO
      END IF

      NFR1 =NFR1 +2*NLST1*LEN1*LEN2

      RETURN
      END
      SUBROUTINE PMLTR1 (KOD,IMLTOP,LST1,
     &                   lg_X,NAS1,NIS1,JXOFF,
     &                   F,NFI,NFJ,
     &                   lg_Y,NAS2,NIS2)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      DIMENSION F(NFI,NFJ)
      DIMENSION LST1(4,NLST1)
#include "sigma.fh"

#ifdef _MOLCAS_MPP_
C SVC: Determine the index ranges of the local chunks of lg_X and lg_Y.
C The boundaries and leading dimension are stored in a common block for
C access inside the lower-level routines.
C For now, only case H is handled as a distributed array, which is
C always the Y array.
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
*     CALL GA_Distribution (lg_X,myRank,iXLo,iXHi,jXLo,jXHi)
*     IF (iXLo.NE.0.AND.jXLo.NE.0) THEN
*       CALL GA_Access (lg_X,iXLo,iXHi,jXLo,jXHi,mX,LDX)
*     END IF
        CALL GA_Distribution (lg_Y,myRank,iYLo,iYHi,jYLo,jYHi)
        IF (iYLo.NE.0.AND.jYLo.NE.0) THEN
          CALL GA_Access (lg_Y,iYLo,iYHi,jYLo,jYHi,mY,LDY)
          IF (KOD.EQ.17 .OR. KOD.EQ.18) THEN
            CALL MLTR1_EH(IMLTOP,LST1,
     &                  WORK(lg_X),NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  DBL_MB(mY),NAS2,jYLo,jYHi)
          ELSE IF (KOD.EQ.21 .OR. KOD.EQ.22) THEN
            CALL MLTR1_GH(IMLTOP,LST1,
     &                  WORK(lg_X),NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  DBL_MB(mY),NAS2,jYLo,jYHi)
          END IF
          CALL GA_Release_Update (lg_Y,iYLo,iYHi,jYLo,jYHi)
        END IF
        CALL GA_Sync()
      ELSE
        IF (KOD.EQ.17 .OR. KOD.EQ.18) THEN
          CALL MLTR1_EH(IMLTOP,LST1,
     &                  WORK(lg_X),NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  WORK(lg_Y),NAS2,1,NIS2)
        ELSE IF (KOD.EQ.21 .OR. KOD.EQ.22) THEN
          CALL MLTR1_GH(IMLTOP,LST1,
     &                  WORK(lg_X),NAS1,NIS1,JXOFF,
     &                  F,NFI,NFJ,
     &                  WORK(lg_Y),NAS2,1,NIS2)
        END IF
      END IF
#else
      IF (KOD.EQ.17 .OR. KOD.EQ.18) THEN
        CALL MLTR1_EH(IMLTOP,LST1,
     &                WORK(lg_X),NAS1,NIS1,JXOFF,
     &                F,NFI,NFJ,
     &                WORK(lg_Y),NAS2,1,NIS2)
      ELSE IF (KOD.EQ.21 .OR. KOD.EQ.22) THEN
        CALL MLTR1_GH(IMLTOP,LST1,
     &                WORK(lg_X),NAS1,NIS1,JXOFF,
     &                F,NFI,NFJ,
     &                WORK(lg_Y),NAS2,1,NIS2)
      END IF
#endif
      RETURN
      END

      SUBROUTINE MLTR1_EH (IMLTOP,LST1,
     &                     X,NAS1,NIS1,JXOFF,
     &                     F,NFT,NFA,
     &                     Y,NAS2,jYLo,jYHi)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(NAS1,NIS1),F(NFT,NFA),Y(NAS2,jYHi-jYLo+1)
      DIMENSION LST1(4,NLST1)
#include "sigma.fh"

C this routine is adapted to use chunks of a distributed array (lg_Y)
C for the H case.  The chunks span all rows (NAS2) and columns jYlo to
C jYHi. The array Y points to the beginning of a chunk.

      IF(IMLTOP.EQ.0) THEN
        NA=INCX3/NAS1
        DO ILST=1,NLST1
        L1=LST1(1,ILST)
        L2=LST1(2,ILST)
        L3=LST1(3,ILST)
        L4=LST1(4,ILST)
        V=VAL1(L4)
        JX=JXOFF+L1
        DO I=1,NAS1
        A=V*F(I,L2)
C X(L1,p,q):= Add V*F(L2,p)*Y(L3,q), p=1..LEN1, q=1..LEN2
        CALL DAXPY_(jYHi-jYLo+1,A,Y(L3,1),NAS2,
     &                 X(I,JX+NA*(jYLo-1)),INCX3)
        END DO
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        NA=INCX3/NAS1
        DO ILST=1,NLST1
        L1=LST1(1,ILST)
        L2=LST1(2,ILST)
        L3=LST1(3,ILST)
        L4=LST1(4,ILST)
        V=VAL1(L4)
        JX=JXOFF+L1
        DO J=jYLo,jYHi
C Y(L3,q):= Add V*F(L2,p)*X(L1,p,q), p=1..LEN1, q=1..LEN2
        Y(L3,J-jYLo+1)=Y(L3,J-jYLo+1)+
     &       V*DDOT_(NAS1,F(1,L2),1,X(1,JX+NA*(J-1)),1)
        END DO
        END DO
      ELSE
        NI=INCX3/NAS1
        DO ILST=1,NLST1
        L1=LST1(1,ILST)
        L2=LST1(2,ILST)
        L3=LST1(3,ILST)
        L4=LST1(4,ILST)
        V=VAL1(L4)
C F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
        JX=JXOFF+L1
        DO J=jYLo,jYHi
        A=V*Y(L3,J-jYLo+1)
        CALL DAXPY_(NAS1,A,X(1,JX+(J-1)*NI),1,F(1,L2),1)
        END DO
        END DO
      END IF

*     NFR1 =NFR1 +2*NLST1*LEN1*LEN2

      RETURN
      END
      SUBROUTINE MLTR1_GH (IMLTOP,LST1,
     &                     X,NAS1,NIS1,JXOFF,
     &                     F,NFT,NFI,
     &                     Y,NAS2,jYLo,jYHi)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(NAS1,NIS1),F(NFT,NFI),Y(NAS2,jYHi-jYLo+1)
      DIMENSION LST1(4,NLST1)
#include "sigma.fh"

C this routine is adapted to use chunks of a distributed array (lg_Y)
C for the H case.  The chunks span all rows (NAS2) and columns jYlo to
C jYHi. The array Y points to the beginning of a chunk.

      IF(IMLTOP.EQ.0) THEN
        DO ILST=1,NLST1
        L1=LST1(1,ILST)
        L2=LST1(2,ILST)
        L3=LST1(3,ILST)
        L4=LST1(4,ILST)
        V=VAL1(L4)
        IF (L3.GE.jYLo .AND. L3.LE.jYHi) THEN
          JX=JXOFF+L1
          JY=L3-jYLo+1
          DO I=1,NAS1
          A=V*F(I,L2)
C X(L1,p,q):= Add V*F(L2,p)*Y(L3,q), p=1..LEN1, q=1..LEN2
          CALL DAXPY_(NAS2,A,Y(1,JY),1,X(I,JX),INCX3)
          END DO
        END IF
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        NI=INCX3/NAS1
        DO ILST=1,NLST1
        L1=LST1(1,ILST)
        L2=LST1(2,ILST)
        L3=LST1(3,ILST)
        L4=LST1(4,ILST)
        V=VAL1(L4)
        IF (L3.GE.jYLo .AND. L3.LE.jYHi) THEN
          JX=JXOFF+L1
          JY=L3-jYLo+1
          DO I=1,NAS2
C Y(L3,q):= Add V*F(L2,p)*X(L1,p,q), p=1..LEN1, q=1..LEN2
          Y(I,JY)=Y(I,JY)+
     &       V*DDOT_(NAS1,F(1,L2),1,X(1,JX+NI*(I-1)),1)
          END DO
        END IF
        END DO
      ELSE
        NI=INCX3/NAS1
        DO ILST=1,NLST1
        L1=LST1(1,ILST)
        L2=LST1(2,ILST)
        L3=LST1(3,ILST)
        L4=LST1(4,ILST)
        V=VAL1(L4)
        IF (L3.GE.jYLo .AND. L3.LE.jYHi) THEN
          JX=JXOFF+L1
          JY=L3-jYLo+1
C F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
          DO I=1,NAS2
          A=V*Y(I,JY)
          CALL DAXPY_(NAS1,A,X(1,JX+NI*(I-1)),1,F(1,L2),1)
          END DO
        END IF
        END DO
      END IF

*     NFR1 =NFR1 +2*NLST1*LEN1*LEN2

      RETURN
      END
