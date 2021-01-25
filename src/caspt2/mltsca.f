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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MLTSCA(IMLTOP,LST1,LST2,X,F,Y)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: MyRank, nProcs, Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),F(*),Y(*)
      DIMENSION LST1(4,NLST1), LST2(4,NLST2)
#include "sigma.fh"

C Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
C four entries called L11,L12,L13,L14 for short, for a given
C item, and with V1=VAL1(L14), and similar for the other list,
C compute, for IMLTOP=0 or 1 respectively,
C     X(L11,L21) := Add V1*V2*F(L12,L22)*Y(L13,L23)
C  or Y(L13,L23) := Add V1*V2*F(L12,L22)*X(L11,L21)
C or for IMLTOP=2, compute
C     F(L12,L22) := Add V1*V2*X(L11,L21)*Y(L13,L23)

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
        DO ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
          L11=LST1(1,ILST1)
          L12=LST1(2,ILST1)
          L13=LST1(3,ILST1)
          L14=LST1(4,ILST1)
          V1=VAL1(L14)
          DO ILST2=1,NLST2
            L21=LST2(1,ILST2)
            L22=LST2(2,ILST2)
            L23=LST2(3,ILST2)
            L24=LST2(4,ILST2)
            V2=VAL2(L24)
            IX=1+INCX1*(L11-1)+INCX2*(L21-1)
            IF=1+INCF1*(L12-1)+INCF2*(L22-1)
            IY=1+INCY1*(L13-1)+INCY2*(L23-1)
            X(IX)=X(IX)+V1*V2*F(IF)*Y(IY)
          END DO
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        DO ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
          L11=LST1(1,ILST1)
          L12=LST1(2,ILST1)
          L13=LST1(3,ILST1)
          L14=LST1(4,ILST1)
          V1=VAL1(L14)
          DO ILST2=1,NLST2
            L21=LST2(1,ILST2)
            L22=LST2(2,ILST2)
            L23=LST2(3,ILST2)
            L24=LST2(4,ILST2)
            V2=VAL2(L24)
            IX=1+INCX1*(L11-1)+INCX2*(L21-1)
            IF=1+INCF1*(L12-1)+INCF2*(L22-1)
            IY=1+INCY1*(L13-1)+INCY2*(L23-1)
            Y(IY)=Y(IY)+V1*V2*F(IF)*X(IX)
          END DO
        END DO
      ELSE
        DO ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
          L11=LST1(1,ILST1)
          L12=LST1(2,ILST1)
          L13=LST1(3,ILST1)
          L14=LST1(4,ILST1)
          V1=VAL1(L14)
          DO ILST2=1,NLST2
            L21=LST2(1,ILST2)
            L22=LST2(2,ILST2)
            L23=LST2(3,ILST2)
            L24=LST2(4,ILST2)
            V2=VAL2(L24)
            IX=1+INCX1*(L11-1)+INCX2*(L21-1)
            IF=1+INCF1*(L12-1)+INCF2*(L22-1)
            IY=1+INCY1*(L13-1)+INCY2*(L23-1)
            F(IF)=F(IF)+V1*V2*X(IX)*Y(IY)
          END DO
        END DO
      END IF

      NFSCA=NFSCA+4*NLST1*NLST2
      RETURN
      END
      SUBROUTINE PMLTSCA(KOD,IMLTOP,LST1,LST2,
     &                   X,NXI,NXA,F,NFI,NFA,
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
      DIMENSION X(NXI,NXA),F(NFI,NFA)
      DIMENSION LST1(4,NLST1), LST2(4,NLST2)
#include "sigma.fh"

#ifdef _MOLCAS_MPP_
C SVC: Determine the index ranges of the local chunks of lg_Y.
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
          IF (KOD.EQ.23 .OR. KOD.EQ.24) THEN
            CALL MLTSCA_DH(IMLTOP,LST1,LST2,
     &                   X,NXI,NXA,F,NFI,NFA,
     &                   DBL_MB(mY),NAS2,jYLo,jYHi)
          ELSE
            WRITE(6,*) 'PMLTSCA: not supposed to be here'
            CALL AbEnd()
          END IF
          CALL GA_Release_Update (lg_Y,iYLo,iYHi,jYLo,jYHi)
        END IF
        CALL GA_Sync()
      ELSE
        IF (KOD.EQ.23 .OR. KOD.EQ.24) THEN
          CALL MLTSCA_DH(IMLTOP,LST1,LST2,
     &                   X,NXI,NXA,F,NFI,NFA,
     &                   WORK(lg_Y),NAS2,1,NIS2)
        ELSE
          WRITE(6,*) 'PMLTSCA: not supposed to be here'
          CALL AbEnd()
        END IF
      END IF
#else
      IF (KOD.EQ.23 .OR. KOD.EQ.24) THEN
        CALL MLTSCA_DH(IMLTOP,LST1,LST2,
     &                 X,NXI,NXA,F,NFI,NFA,
     &                 WORK(lg_Y),NAS2,1,NIS2)
      ELSE
        WRITE(6,*) 'PMLTSCA: not supposed to be here'
        CALL AbEnd()
      END IF
#endif
      RETURN
      END
      SUBROUTINE MLTSCA_DH(IMLTOP,LST1,LST2,
     &                     X,NXI,NXA,F,NFI,NFA,
     &                     Y,NAS2,jYLo,jYHi)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(NXI,NXA),F(NFI,NFA),Y(NAS2,jYHi-jYLo+1)
      DIMENSION LST1(4,NLST1), LST2(4,NLST2)
#include "sigma.fh"

C Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
C four entries called L11,L12,L13,L14 for short, for a given
C item, and with V1=VAL1(L14), and similar for the other list,
C compute, for IMLTOP=0 or 1 respectively,
C     X(L11,L21) := Add V1*V2*F(L12,L22)*Y(L13,L23)
C  or Y(L13,L23) := Add V1*V2*F(L12,L22)*X(L11,L21)
C or for IMLTOP=2, compute
C     F(L12,L22) := Add V1*V2*X(L11,L21)*Y(L13,L23)

      IF(IMLTOP.EQ.0) THEN
        DO ILST1=1,NLST1
        L11=LST1(1,ILST1)
        L12=LST1(2,ILST1)
        L13=LST1(3,ILST1)
        L14=LST1(4,ILST1)
        V1=VAL1(L14)
        IF (L13.GE.jYLo .AND. L13.LE.jYHi) THEN
          JY=L13-jYLo+1
          DO ILST2=1,NLST2
          L21=LST2(1,ILST2)
          L22=LST2(2,ILST2)
          L23=LST2(3,ILST2)
          L24=LST2(4,ILST2)
          V2=VAL2(L24)
          X(L11,L21)=X(L11,L21)+V1*V2*F(L12,L22)*Y(L23,JY)
          END DO
        END IF
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        DO ILST1=1,NLST1
        L11=LST1(1,ILST1)
        L12=LST1(2,ILST1)
        L13=LST1(3,ILST1)
        L14=LST1(4,ILST1)
        V1=VAL1(L14)
        IF (L13.GE.jYLo .AND. L13.LE.jYHi) THEN
          JY=L13-jYLo+1
          DO ILST2=1,NLST2
          L21=LST2(1,ILST2)
          L22=LST2(2,ILST2)
          L23=LST2(3,ILST2)
          L24=LST2(4,ILST2)
          V2=VAL2(L24)
          Y(L23,JY)=Y(L23,JY)+V1*V2*F(L12,L22)*X(L11,L21)
          END DO
        END IF
        END DO
      ELSE
        DO ILST1=1,NLST1
        L11=LST1(1,ILST1)
        L12=LST1(2,ILST1)
        L13=LST1(3,ILST1)
        L14=LST1(4,ILST1)
        V1=VAL1(L14)
        IF (L13.GE.jYLo .AND. L13.LE.jYHi) THEN
          JY=L13-jYLo+1
          DO ILST2=1,NLST2
          L21=LST2(1,ILST2)
          L22=LST2(2,ILST2)
          L23=LST2(3,ILST2)
          L24=LST2(4,ILST2)
          V2=VAL2(L24)
          F(L12,L22)=F(L12,L22)+V1*V2*X(L11,L21)*Y(L23,JY)
          END DO
        END IF
        END DO
      END IF

*     NFSCA=NFSCA+4*NLST1*NLST2
      RETURN
      END
