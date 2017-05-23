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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Integer Function MemSO1(lOper,iCmp,jCmp,iShell,jShell)
************************************************************************
*  Object: to compile the number of SO block which will be generated   *
*          by the current shell doublet.                               *
*          "lOper" is the irreducible representation of which the      *
*          operator of the one-electron matrix element belongs.        *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             February '90                                             *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden.                                         *
*             Modified to general non-symmetric operators January '91  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
*
      MemSO1 = 0
      Do 100 j1 = 0, nIrrep-1
         Do 110 i1 = 1, iCmp
            If (iAnd(IrrCmp(IndS(iShell)+i1),2**j1).eq.0) Go To 110
*old        Do 200 j2 = 0, j1
            Do 200 j2 = 0, nIrrep-1
               j12=iEor(j1,j2)
               If (iAnd(lOper,2**j12).eq.0) Go To 200
               jCmpMx = jCmp
               If (iShell.eq.jShell .and. j1.eq.j2) jCmpMx = i1
               Do 210 i2 = 1, jCmpMx
                  If (iAnd(IrrCmp(IndS(jShell)+i2),2**j2).eq.0)
     &               Go To 210
                  MemSO1 = MemSO1 + 1
 210           Continue
 200        Continue
 110     Continue
 100  Continue
*
      Return
      End
