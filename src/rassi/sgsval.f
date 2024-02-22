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
      SUBROUTINE SGSVAL(SGS,NSYM,NLEV,NVERT,
     &              MIDLEV,MVSTA,MVEND,LMAW,LLTV)
      use Struct, only: SGStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
C Purpose: Dereference the Split Graph structure array
C and return values and pointers.

      NSYM  =SGS%nSym
      NLEV  =SGS%nLev
      NVERT =SGS%nVert
      MIDLEV=SGS%MidLev
      MVSTA =SGS%MVSta
      MVEND =SGS%MVEnd
      LMAW  =SGS%lMAW
      LLTV  =SGS%lLTV

      END SUBROUTINE SGSVAL
