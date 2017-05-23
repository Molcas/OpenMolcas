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
************************************************************************
      SUBROUTINE CHEL(IA,IB,IIM,IEL,ISTOP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IR=IIM-1
C     CHECK FOR A=0 , B=IEL
      IRR=IR-IA
      IF(IRR.LT.0)GO TO 50
      IF(IRR.GE.IB-IEL)GO TO 100
50    IF(IEL.EQ.1)GO TO 90
C     CHECK FOR A=1 , B=IEL-2
      IRR=IR-IA+1
      IF(IRR.LT.0)GO TO 90
      IF(IRR.GE.IB-IEL+2)GO TO 100
90    ISTOP=1
      RETURN
100   ISTOP=0
      RETURN
      END
