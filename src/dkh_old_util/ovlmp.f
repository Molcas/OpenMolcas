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
      REAL*8 FUNCTION OVLMP (NP,ZP,NQ,ZQ)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C         <NP,ZP|NQ,ZQ>
C     Note that both NP and NQ should be either even or odd, but
C     not mixed.
C
C...  auxiliar constant pool:       ready only up to g-valence/g-core
      PARAMETER (lp1=5,lp12=lp1*lp1,lp13=(lp1*lp1+lp1)/2)
      COMMON/CONST/RCA(lp1,lp13),DFAC(lp12),KOSUU(lp13),NYU(lp1,lp13)
      IT11=2*NP
      IT22=2*NQ
      RT1=sqrt(ZP**(IT11+1)*ZQ**(IT22+1))/(DFAC(IT11)*DFAC(IT22))
      RT1=sqrt(RT1)
      IT33=NP+NQ
      ZAV=0.5D0*(ZP+ZQ)
      OVLMP=RT1*DFAC(IT33)/sqrt(ZAV**(IT33+1))
      RETURN
      END
