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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE DIIS_CPF(DPI,DPJ,BST,MIT,BIJ,ITP,CN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DPI(*),DPJ(*),BST(MIT,MIT),BIJ(ITP,ITP),
     *CN(*),WHS(50)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
*
      IF(ITPUL.EQ.1)GO TO 26
C
      ITM=ITPUL-1
      DO 5 I=1,ITM
         DO 10 J=1,ITM
            BIJ(J,I)=BST(J,I)
10       CONTINUE
5     CONTINUE
      DO 15 I=1,ITPUL
         BIJ(ITP,I)=-1.0D0
         BIJ(I,ITP)=-1.0D0
15    CONTINUE
      BIJ(ITP,ITP)=0.0D0
C
      DO 20 I=1,ITM
         IAD=IADDP(I+1)
         CALL dDAFILE(Lu_CI,2,DPJ,NCONF,IAD)
         T=DDOT_(NCONF,DPI,1,DPJ,1)
         BIJ(I,ITPUL)=T
         BIJ(ITPUL,I)=T
         BST(I,ITPUL)=T
         BST(ITPUL,I)=T
         IF(I.EQ.1) THEN
            T=DDOT_(NCONF,DPJ,1,DPJ,1)
            BIJ(1,1)=T
            BST(1,1)=T
         END IF
20    CONTINUE
      BIJ(ITPUL,ITPUL)=DDOT_(NCONF,DPI,1,DPI,1)
      BST(ITPUL,ITPUL)=BIJ(ITPUL,ITPUL)
C
      IF(IPRINT.LT.10)GO TO 26
      DO 17 I=1,ITP
      WRITE(6,16)(BIJ(J,I),J=1,ITP)
      CALL XFLUSH(6)
16    FORMAT(6X,'BIJ ',6F12.6)
17    CONTINUE
C
26    IF(IDIIS.EQ.1)GO TO 25
      DO 7 I=1,ITPUL
         IAD=IADDP(I)
         CALL dDAFILE(Lu_CI,2,DPJ,NCONF,IAD)
         DO 6 J=1,NCONF
            DPI(J)=DPI(J)+DPJ(J)
6        CONTINUE
7     CONTINUE
      IF(IPRINT.GE.15)WRITE(6,14)(DPI(I),I=1,NCONF)
14    FORMAT(6X,'C(DIIS)',5F10.6)
      RETURN
C
25    CALL DECOMP(ITP,BIJ,BIJ)
      DO 55 I=1,ITPUL
      WHS(I)=0.0D00
55    CONTINUE
      WHS(ITP)=-1.0D00
      CALL SOLVE(ITP,BIJ,WHS,CN)
C
C     UPDATE P AND DELTA P
C
      CALL NEXT(DPI,DPJ,CN)
C
      ITPUL=0
      RETURN
      END
