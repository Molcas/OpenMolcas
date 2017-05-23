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
      SUBROUTINE LOOP70(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,BUF,
     *      IBUF,A,B,F,FSEC,IPOF,IPOA,IPOB,
     * MYL,NYL,INDA,INDB,INMY,INNY,IFTB,IFTA,FACS,
     * IAB,CPL,CPLA, NVIRA,NVIRC,NVIRB)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
      DIMENSION INTSYM(*),INDX(*),C(*),S(*),
     *          ABIJ(*),AIBJ(*),AJBI(*),
     *          BUF(NBITM3),IBUF(NBITM3+2),
     *          A(*),B(*),F(*),FSEC(*)
      DIMENSION IPOF(9),IPOA(9),IPOB(9)


      DO 70 IASYM=1,NSYM
      IAB=IPOF(IASYM+1)-IPOF(IASYM)
      IF(IAB.EQ.0)GO TO 70
      ICSYM=MUL(MYL,IASYM)
      IBSYM=MUL(NYL,ICSYM)
      IF(INDA.EQ.INDB.AND.IBSYM.GT.IASYM)GO TO 70
      NVIRC=NVIR(ICSYM)
      IF(NVIRC.EQ.0)GO TO 70
      NVIRA=NVIR(IASYM)
      NVIRB=NVIR(IBSYM)
      IF(ICSYM.GE.IASYM)GO TO 31
      IF(ICSYM.GE.IBSYM)GO TO 32
C     CASE 1, IASYM > ICSYM AND IBSYM > ICSYM
      IPF=IPOF(IASYM)+1
      CALL DYAX(IAB,CPL,AIBJ(IPF),1,F,1)
      CALL DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
      IF(INDA.EQ.INDB)CALL SETZZ(F,NVIRA)
C     CALL VCLR(A,1,NBC)
C     CALL FMMM(C(INMY+IPOA(IASYM)),F,A,NVIRC,NVIRB,NVIRA)
C     CALL DAXPY_(NBC,FACS,A,1,S(INNY+IPOB(IBSYM)),1)
      CALL DGEMM_('N','N',NVIRC,NVIRB,NVIRA,
     *           FACS,C(INMY+IPOA(IASYM)),NVIRC,
     *           F,NVIRA,1.0D00,S(INNY+IPOB(IBSYM)),NVIRC)
      IF(INDA.EQ.INDB)GO TO 70
      IPF=IPOF(IBSYM)+1
      CALL VCLR(F,1,IAB)
      CALL DYAX(IAB,CPL,AJBI(IPF),1,F,1)
      CALL DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
      CALL DGEMM_('N','N',NVIRC,NVIRA,NVIRB,
     *           FACS,C(INNY+IPOB(IBSYM)),NVIRC,
     *           F,NVIRB,1.0D00,S(INMY+IPOA(IASYM)),NVIRC)
      GO TO 70
C     CASE 2, IASYM > ICSYM AND ICSYM > OR = IBSYM
32    IPF=IPOF(IBSYM)+1
      CALL DYAX(IAB,CPL,AJBI(IPF),1,F,1)
      CALL DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
C
      IF(NYL.EQ.1) THEN
       CALL DGEMM_('N','T',NVIRB,NVIRC,NVIRA,
     *            FACS,F,NVIRB,C(INMY+IPOA(IASYM)),NVIRC,
     *            0.D0,A,NVIRB)
       IF(IFTB.EQ.1) THEN
        CALL TRADD(A,S(INNY+IPOB(ICSYM)),NVIRB)
        CALL SQUARN(C(INNY+IPOB(IBSYM)),A,NVIRB)
       ELSE
        CALL SIADD(A,S(INNY+IPOB(ICSYM)),NVIRB)
        CALL SQUAR(C(INNY+IPOB(IBSYM)),A,NVIRB)
       ENDIF
       CALL DGEMM_('N','N',NVIRC,NVIRA,NVIRB,
     *            FACS,A,NVIRC,F,NVIRB,
     *            1.D0,S(INMY+IPOA(IASYM)),NVIRC)
      ELSE
       FACSX=FACS
       IF(IFTB.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','T',NVIRB,NVIRC,NVIRA,
     *            FACSX,F,NVIRB,C(INMY+IPOA(IASYM)),NVIRC,
     *            1.D0,S(INNY+IPOB(ICSYM)),NVIRB)
       CALL DGEMM_('T','N',NVIRC,NVIRA,NVIRB,
     *            FACSX,C(INNY+IPOB(ICSYM)),NVIRB,F,NVIRB,
     *            1.D0,S(INMY+IPOA(IASYM)),NVIRC)
      ENDIF
      GO TO 70
C UPDATED UNTIL HERE
31    IF(ICSYM.GE.IBSYM)GO TO 33
C     CASE 3, ICSYM > OR = IASYM AND IBSYM > ICSYM
      IPF=IPOF(IASYM)+1
      CALL DYAX(IAB,CPL,AIBJ(IPF),1,F,1)
      CALL DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
      IF(MYL.EQ.1) THEN
       IF(IFTA.EQ.0)CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
       IF(IFTA.EQ.1)CALL SQUARN(C(INMY+IPOA(IASYM)),A,NVIRA)
       CALL DGEMM_('N','N',NVIRC,NVIRB,NVIRA,
     *            FACS,A,NVIRC,F,NVIRA,
     *            1.D0,S(INNY+IPOB(IBSYM)),NVIRC)
       CALL DGEMM_('N','T',NVIRA,NVIRC,NVIRB,
     *            FACS,F,NVIRA,C(INNY+IPOB(IBSYM)),NVIRC,
     *            0.D0,A,NVIRA)
       IF(IFTA.EQ.0) CALL SIADD(A,S(INMY+IPOA(IASYM)),NVIRA)
       IF(IFTA.EQ.1) CALL TRADD(A,S(INMY+IPOA(IASYM)),NVIRA)
      ELSE
       FACSX=FACS
       IF(IFTA.EQ.1) FACSX=-FACS
       CALL DGEMM_('T','N',NVIRC,NVIRB,NVIRA,
     *            FACSX,C(INMY+IPOA(ICSYM)),NVIRA,F,NVIRA,
     *            1.D0,S(INNY+IPOB(IBSYM)),NVIRC)
       CALL DGEMM_('N','T',NVIRA,NVIRC,NVIRB,
     *            FACSX,F,NVIRA,C(INNY+IPOB(IBSYM)),NVIRC,
     *            1.D0,S(INMY+IPOA(ICSYM)),NVIRA)
      ENDIF
      GO TO 70
C     CASE 4, ICSYM > OR = IASYM AND ICSYM > OR = IBSYM
33    IPF=IPOF(IBSYM)+1
      CALL DYAX(IAB,CPL,AJBI(IPF),1,F,1)
      CALL DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
      IF(INDA.EQ.INDB)CALL SETZZ(F,NVIRA)
      IF(MYL.EQ.1.AND.NYL.EQ.1) THEN
C
       IF(IFTA.EQ.0) CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
       IF(IFTA.EQ.1) CALL SQUARM(C(INMY+IPOA(IASYM)),A,NVIRA)
       CALL DGEMM_('N','N',NVIRB,NVIRC,NVIRA,
     *            FACS,F,NVIRB,A,NVIRA,
     *            0.D0,B,NVIRB)
       IF(IFTB.EQ.0) CALL SIADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
       IF(IFTB.EQ.1) CALL TRADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
C
      ELSE IF (MYL.EQ.1.AND.NYL.NE.1) THEN
       IF(IFTA.EQ.0) CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIRA)
       IF(IFTA.EQ.1) CALL SQUARM(C(INMY+IPOA(IASYM)),A,NVIRA)
       FACSX=FACS
       IF(IFTB.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','N',NVIRB,NVIRC,NVIRA,
     *            FACSX,F,NVIRB,A,NVIRA,
     *            1.D0,S(INNY+IPOB(ICSYM)),NVIRB)
C
      ELSE IF (MYL.NE.1.AND.NYL.EQ.1) THEN
C
       FACSX=FACS
       IF(IFTA.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','N',NVIRB,NVIRC,NVIRA,
     *            FACSX,F,NVIRB,C(INMY+IPOA(ICSYM)),NVIRA,
     *            0.D0,B,NVIRB)
       IF(IFTB.EQ.0) CALL SIADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
       IF(IFTB.EQ.1) CALL TRADD(B,S(INNY+IPOB(ICSYM)),NVIRB)
      ELSE IF (MYL.NE.1.AND.NYL.NE.1) THEN
       FACSX=FACS
       IF(IFTA+IFTB.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','N',NVIRB,NVIRC,NVIRA,
     *            FACSX,F,NVIRB,C(INMY+IPOA(ICSYM)),NVIRA,
     *            1.D0,S(INNY+IPOB(ICSYM)),NVIRB)
      ENDIF
      IF(INDA.EQ.INDB)GO TO 70
      IPF=IPOF(IASYM)+1
      CALL DYAX(IAB,CPL,AIBJ(IPF),1,F,1)
      CALL DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
C
      IF(NYL.EQ.1.AND.MYL.EQ.1) THEN
C
       IF(IFTB.EQ.0) CALL SQUAR(C(INNY+IPOB(IBSYM)),A,NVIRB)
       IF(IFTB.EQ.1) CALL SQUARM(C(INNY+IPOB(IBSYM)),A,NVIRB)
       CALL DGEMM_('N','N',NVIRA,NVIRC,NVIRB,
     *            FACS,F,NVIRA,A,NVIRB,
     *            0.D0,B,NVIRA)
       IF(IFTA.EQ.0) CALL SIADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
       IF(IFTA.EQ.1) CALL TRADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
C
      ELSE IF (NYL.EQ.1.AND.MYL.NE.1) THEN
       IF(IFTB.EQ.0) CALL SQUAR(C(INNY+IPOB(ICSYM)),A,NVIRB)
       IF(IFTB.EQ.1) CALL SQUARM(C(INNY+IPOB(ICSYM)),A,NVIRB)
       FACSX=FACS
       IF(IFTA.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','N',NVIRA,NVIRC,NVIRB,
     *            FACSX,F,NVIRA,A,NVIRB,
     *            1.D0,S(INMY+IPOA(ICSYM)),NVIRA)
C
      ELSE IF (NYL.NE.1.AND.MYL.EQ.1) THEN
C
       FACSX=FACS
       IF(IFTB.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','N',NVIRA,NVIRC,NVIRB,
     *            FACSX,F,NVIRA,C(INNY+IPOB(ICSYM)),NVIRB,
     *            0.D0,B,NVIRA)
       IF(IFTA.EQ.0) CALL SIADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
       IF(IFTA.EQ.1) CALL TRADD(B,S(INMY+IPOA(ICSYM)),NVIRA)
C
      ELSE IF (NYL.NE.1.AND.MYL.NE.1) THEN
C
       FACSX=FACS
       IF(IFTA+IFTB.EQ.1) FACSX=-FACS
       CALL DGEMM_('N','N',NVIRA,NVIRC,NVIRB,
     *            FACSX,F,NVIRA,C(INNY+IPOB(ICSYM)),NVIRB,
     *            1.D0,S(INMY+IPOA(ICSYM)),NVIRA)
      ENDIF
C
70    CONTINUE
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(INTSYM)
        CALL Unused_integer_array(INDX)
        CALL Unused_integer_array(IBUF)
        CALL Unused_real_array(FSEC)
        CALL Unused_real_array(BUF)
      END IF
      END
      subroutine faibj5(LENBUF,JTURN,IBUF,BUF, AIBJ,ABIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
      DIMENSION ABIJ(NVSQ),AIBJ(NVSQ),
     *          BUF(NBITM3),IBUF(NBITM3+2)
       IF(LENBUF.GT.0) THEN
           IF(JTURN.EQ.1) THEN
             do i=1,LENBUF
               aibj(IBUF(i))=buf(i)
             enddo
           ELSE
             do i=1,LENBUF
               abij(IBUF(i))=buf(i)
             enddo
           END IF
         END IF
       return
       end
