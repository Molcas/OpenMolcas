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
      SUBROUTINE DENS_CPF(C,D,ICASE,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),D(*)
      DIMENSION ICASE(*)

#include "SysDef.fh"

#include "cpfmcpf.fh"
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CRL   JO(L)=IAND(ISHFT(QOCC((L+29)/30),-2*((L+29)/30*30-L)),3)
CPAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICASE,L)
C
      ILIM=NORBT*(NORBT+1)/2
      CALL SETZ(D,ILIM)
      C(IREF0)=D0
CRL   CALL DOTPR(C,1,C,1,A,NCONF)
      A=DDOT_(NCONF,C,1,C,1)
      WRITE(6,20)A
      CALL XFLUSH(6)
20    FORMAT(5X,'SUM OF SQUARED CPX(BAR)',F10.4)
      C(IREF0)=D1
      EMA=D1-A
      II1=(IREF0-1)*LN
      DO 5 I=1,LN
      JOJ=JO(II1+I)
      IF(JOJ.GE.2)JOJ=JOJ-1
      II=I*(I+1)/2
      D(II)=JOJ*EMA
5     CONTINUE
      RETURN
      END
