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
C
      SUBROUTINE SECORD(A,B,C,FAC,NAL,NBL,NSIJ,IFT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),C(*)
      IAB=0
      NAA=0
      DO 10 NA=1,NAL
      NBB=0
      NA1=NBL
      IF(NSIJ.EQ.1)NA1=NA-1
      IF(NA1.EQ.0)GO TO 15
      DO 20 NB=1,NA1
      IAB=IAB+1
      IF(IFT.EQ.0)C(IAB)=B(NAA+NB)+A(NBB+NA)
      IF(IFT.EQ.1)C(IAB)=B(NAA+NB)-A(NBB+NA)
      NBB=NBB+NAL
20    CONTINUE
15    IF(NSIJ.NE.1)GO TO 35
      IAB=IAB+1
      C(IAB)=FAC*A(NAA+NA)
35    NAA=NAA+NBL
10    CONTINUE
      RETURN
      END
