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
      Subroutine NonSym(nStab,jStab,A,Tx)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 A(3),Tx(3)
      Integer   jStab(0:nStab-1)
*
      Do iStab = 0, nStab-1
         If (A(1).ne.Zero.and.iAnd(jStab(iStab),1).ne.0) Go To 10
         If (A(2).ne.Zero.and.iAnd(jStab(iStab),2).ne.0) Go To 10
         If (A(3).ne.Zero.and.iAnd(jStab(iStab),4).ne.0) Go To 10
         If (iAnd(jStab(iStab),1).ne.0) Tx(1)=Zero
         If (iAnd(jStab(iStab),2).ne.0) Tx(2)=Zero
         If (iAnd(jStab(iStab),4).ne.0) Tx(3)=Zero
10       Continue
      End Do
*
      Return
      End
