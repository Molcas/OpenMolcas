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
      SUBROUTINE MLTDXP(IMLTOP,LST1,LST2,X,F,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),F(*),Y(*)
      DIMENSION LST1(4,NLST1), LST2(4,NLST2)
#include "sigma.fh"
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#endif

C Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
C four entries called L11,L12,L13,L14 for short, for a given
C item, and with V1=VAL1(L14), and similar for the other list,
C compute, for IMLTOP=0 or 1 resp.,
C     X(L11,L21,a) := Add V1*V2*F(L12,L22)*Y(L13,L23,a)
C or  Y(L13,L23,a) := Add V1*V2*F(L12,L22)*X(L11,L21,a), a=1,LEN1
C or for IMLTOP=2, compute
C     F(L12,L22) := Add V1*V2*X(L11,L21,a)*Y(L13,L23,a)
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
            IY=1+INCY1*(L13-1)+INCY2*(L23-1)
            IF=1+INCF1*(L12-1)+INCF2*(L22-1)
            A=V1*V2*F(IF)
            CALL DAXPY_(LEN1,A,Y(IY),INCY3,X(IX),INCX3)
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
            IY=1+INCY1*(L13-1)+INCY2*(L23-1)
            IF=1+INCF1*(L12-1)+INCF2*(L22-1)
            A=V1*V2*F(IF)
            CALL DAXPY_(LEN1,A,X(IX),INCX3,Y(IY),INCY3)
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
            IY=1+INCY1*(L13-1)+INCY2*(L23-1)
            IF=1+INCF1*(L12-1)+INCF2*(L22-1)
            A=V1*V2*F(IF)
            V=V1*V2
            F(IF)=F(IF)+V*DDOT_(LEN1,X(IX),INCX3,Y(IY),INCY3)
          END DO
        END DO
      END IF
      NFDXP=NFDXP+2*NLST1*NLST2*LEN1
      RETURN
      END
