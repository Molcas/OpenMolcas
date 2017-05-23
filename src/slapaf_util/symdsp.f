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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      Logical Function SymDsp(iOper,nIrrep,iBsFnc)
************************************************************************
*                                                                      *
* Object: to establish if a translation or a rotation belongs to the   *
*         total symmetric irreducible representation.                  *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              ICopy                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             April '92                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer   iOper(0:nIrrep-1), jPrmt(0:7)
#include "real.fh"
      Data jPrmt/1,-1,-1,1,-1,1,1,-1/
*
*     Statement function
*
      iPrmt(i,j) = jPrmt(iAnd(i,j))
*
*     Write (*,*) ' iBsFnc=',iBsFnc
      SymDsp = .True.
      mask = 0
      Do 20 i = 0, nIrrep-1
         Do 21 j = 1, 3
            If (iAnd(iOper(i),2**(j-1)).ne.0) mask = iOr(mask,2**(j-1))
 21      Continue
 20   Continue
      jBsFnc = iAnd(mask,iBsFnc)
*     Write (*,*) ' jBsFnc=',jBsFnc
*
*     Loop over operators
*
      iAcc = 0
      Do 10 i = 0, nIrrep-1
         iAcc = iAcc + iPrmt(iOper(i),jBsFnc)
 10   Continue
      If (iAcc.eq.0) SymDsp = .False.
*
      Return
      End
