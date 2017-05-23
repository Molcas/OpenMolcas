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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      SubRoutine rKappa_Zeta(rKappa,Zeta,nZeta)
************************************************************************
*                                                                      *
* Object: modify rkappa                                                *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, SWEDEN                               *
*             September '00.                                           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Zeta(nZeta), rKappa(nZeta)
*
      exp32 = -Three/Two
      Do iZeta = 1, nZeta
         rKappa(iZeta) = rKappa(iZeta) * Zeta(iZeta)**exp32
      End Do
*
      Return
      End
