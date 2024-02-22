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
      SUBROUTINE SGSVAL(SGS,NSYM,NLEV)
      use Struct, only: SGStruct
      IMPLICIT None
      Type (SGStruct) SGS
      Integer NSYM,NLEV
C Purpose: Dereference the Split Graph structure user defined type
C and return values.

      NSYM  =SGS%nSym
      NLEV  =SGS%nLev

      END SUBROUTINE SGSVAL
