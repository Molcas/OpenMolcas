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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_GetCutInt(Value)
C
C     Thomas Bondo Pedersen, November 2010.
C
C     Returns the value of CutInt.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 Value
#include "itmax.fh"
#include "info.fh"
      Value=CutInt
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SetCutInt(Value)
C
C     Thomas Bondo Pedersen, November 2010.
C
C     Sets the value of CutInt.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 Value
#include "itmax.fh"
#include "info.fh"
      CutInt=Value
      End
