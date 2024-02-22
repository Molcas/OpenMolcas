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
      SUBROUTINE SGSVAL(ISGS,SGS,NSYM,NLEV,LISM,NVERT,LDRT,
     &              LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
      use Struct, only: nSGSize, SGStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
      Dimension ISGS(nSGSize)
C Purpose: Dereference the Split Graph structure array
C and return values and pointers.

      NSYM  =ISGS(1)
      NLEV  =ISGS(2)
      LISM  =ISGS(3)
      NVERT =ISGS(4)
      LDRT  =ISGS(5)
      LDOWN =ISGS(6)
      LUP   =ISGS(7)
      MIDLEV=ISGS(8)
      MVSTA =ISGS(9)
      MVEND =ISGS(10)
      LMAW  =ISGS(11)
      LLTV  =ISGS(12)

      NSYM  =SGS%nSym
      NLEV  =SGS%nLev
      LISM  =SGS%lISm
      NVERT =SGS%nVert
      LDRT  =SGS%LDRT
      LDOWN =SGS%lDOWN
      LUP   =SGS%lUP
      MIDLEV=SGS%MidLev
      MVSTA =SGS%MVSta
      MVEND =SGS%MVEnd
      LMAW  =SGS%lMAW
      LLTV  =SGS%lLTV

      END SUBROUTINE SGSVAL
