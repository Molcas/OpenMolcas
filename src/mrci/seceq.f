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
      SUBROUTINE SECEQ(A,B,C,NAL,IFT,FAC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NAL,NAL),B(NAL,NAL),C((NAL*(NAL+1))/2)
      IF(IFT.EQ.0) THEN
        IAB=0
        DO 20 NA=1,NAL
          DO 10 NB=1,NA-1
            IAB=IAB+1
            C(IAB)=B(NB,NA)+A(NA,NB)
10        CONTINUE
          IAB=IAB+1
          C(IAB)=FAC*A(NA,NA)
20      CONTINUE
      ELSE
        IAB=0
        DO 40 NA=1,NAL
          DO 30 NB=1,NA-1
            IAB=IAB+1
            C(IAB)=B(NB,NA)-A(NA,NB)
30        CONTINUE
          IAB=IAB+1
          C(IAB)=0.0D00
40      CONTINUE
      END IF
      RETURN
      END
