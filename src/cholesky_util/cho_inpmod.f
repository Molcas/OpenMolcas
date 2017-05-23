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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_InpMod(Mode)
C
C     Thomas Bondo Pedersen, Jan. 2005.
C
C     Purpose: modifiy Cholesky settings for
C              Mode = 'LOW ' : low-accuracy decomposition
C              Mode = 'MEDI' : medium-accuracy decomposition
C              Mode = 'HIGH' : high-accuracy decomposition
C              Mode = '1CCD' : set one-center approximation
C              (All other Mode-values are ignored.)
C
      Implicit None
      Character*4 Mode

      Character*4 Mod2

      Mod2 = Mode
      Call Upcase(Mod2)

      If (Mod2(1:3) .eq. 'LOW') Then
         Call Cho_SetDecompositionThreshold(1.0d-4)
      Else If (Mod2(1:4) .eq. 'MEDI') Then
         Call Cho_SetDecompositionThreshold(1.0d-6)
      Else If (Mod2(1:4) .eq. 'HIGH') Then
         Call Cho_SetDecompositionThreshold(1.0d-8)
      Else If (Mod2(1:4) .eq. '1CCD') Then
         Call Cho_Set1CCD(.True.)
      End If

      End
      Subroutine Cho_SetDecompositionThreshold(Thr)
      Implicit None
      Real*8 Thr
#include "cholesky.fh"

      ThrCom=Thr

      End
      Subroutine Cho_Set1CCD(do1CCD)
      Implicit None
      Logical do1CCD
#include "cholesky.fh"

      Cho_1Center=do1CCD

      End
      Subroutine Cho_SetSpan(val)
      Implicit None
      Real*8 val
#include "cholesky.fh"

      Span=val

      End
