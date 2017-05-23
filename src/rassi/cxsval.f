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
      SUBROUTINE CXSVAL(ICIS,IXS,NMIDV,NIPWLK,LNOW,LIOW,LNCSF,
     &                 LNOCSF,LIOCSF,NWALK,LICASE,
     &                 MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,
     &           LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='CXSVAL')
#include "Struct.fh"
      Dimension ICIS(nCISize)
      Dimension IXS(nXSize)

C Purpose: Dereference the CI structure and the excitation
C structure arrays and return values and pointers.

      CALL QENTER(ROUTINE)

C CI structure, sizes, addresses...
      nMidV =ICIS(1)
      nIpWlk =ICIS(2)
      lNOW =ICIS(3)
      lIOW =ICIS(4)
      lNCSF =ICIS(5)
      lNOCSF =ICIS(6)
      lIOCSF =ICIS(7)
      nWalk =ICIS(8)
      lICase =ICIS(9)
C Excitation operators, coupling coefficients,...
      MxEO =IXS(1)
      lNOCP =IXS(2)
      lIOCP =IXS(3)
      nICoup =IXS(4)
      lICoup =IXS(5)
      nVTab =IXS(6)
      lVTab =IXS(7)
      lMVL =IXS(8)
      lMVR =IXS(9)
      NT1MX =IXS(10)
      NT2MX =IXS(11)
      NT3MX =IXS(12)
      NT4MX =IXS(13)
      NT5MX =IXS(14)

      CALL QEXIT(ROUTINE)
      RETURN
      END
