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
      SUBROUTINE MLTMV (IMLTOP,LST1,X,F,Y)
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
C compute, for IMLTOP=0 or 1, respectively.
C    X(L1,i) := Add V*F(L2,a)*Y(L3,i,a), i=1..LEN1, a=1..LEN2
C or Y(L3,i,a):= Add V*F(L2,a)*X(L1,i), i=1..LEN1, a=1..LEN2
C or for IMLTOP=2, compute
C     F(L2,a) := Add V*X(L1,i)*Y(L3,i,a)
C However, strides etc can vary from case to case: The indices
C may appear in 1st, 2nd or 3rd position, Fortran-address-wise,
C thus an independent stride is to be given for each of the
C indices in the above 'formal' index ordering, as follows:
C The formal X(p,q,r) is accessed as
C X(1+INCX1*(p-1)+INCX2*(q-1)+INCX3*(r-1)), etc.

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
          L1=LST1(1,ILST1)
          L2=LST1(2,ILST1)
          L3=LST1(3,ILST1)
          L4=LST1(4,ILST1)
          V=VAL1(L4)
          IX=INCX1*(L1-1)+1
          IF=INCF1*(L2-1)+1
          IY=INCY1*(L3-1)+1
C    X(L1,i) := Add V*F(L2,a)*Y(L3,i,a), i=1..LEN1, a=1..LEN2
          DO I=1,LEN1
            X(IX)=X(IX)+V*DDOT_(LEN2,F(IF),INCF2,Y(IY),INCY3)
            IX=IX+INCX2
            IY=IY+INCY2
          END DO
        END DO
      ELSE IF(IMLTOP.EQ.1) THEN
        DO ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
          L1=LST1(1,ILST1)
          L2=LST1(2,ILST1)
          L3=LST1(3,ILST1)
          L4=LST1(4,ILST1)
          V=VAL1(L4)
          IX=INCX1*(L1-1)+1
          IF=INCF1*(L2-1)+1
          IY=INCY1*(L3-1)+1
C or Y(L3,i,a):= Add V*F(L2,a)*X(L1,i), i=1..LEN1, a=1..LEN2
          DO I=1,LEN2
            A=V*F(IF)
            CALL DAXPY_(LEN1,A,X(IX),INCX2,Y(IY),INCY2)
            IY=IY+INCY3
            IF=IF+INCF2
          END DO
        END DO
      ELSE
        DO ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
          L1=LST1(1,ILST1)
          L2=LST1(2,ILST1)
          L3=LST1(3,ILST1)
          L4=LST1(4,ILST1)
          V=VAL1(L4)
          IX=INCX1*(L1-1)+1
          IF=INCF1*(L2-1)+1
          IY=INCY1*(L3-1)+1
C     F(L2,a) := Add V*X(L1,i)*Y(L3,i,a)
          DO I=1,LEN1
            A=V*X(IX)
            CALL DAXPY_(LEN2,A,Y(IY),INCY3,F(IF),INCF2)
            IX=IX+INCX2
            IY=IY+INCY2
          END DO
        END DO
      END IF

      NFMV =NFMV +2*NLST1*LEN1*LEN2

      RETURN
      END
