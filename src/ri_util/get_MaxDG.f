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
      Subroutine Get_maxDG(SDG,nnSkal,MxBasSh)
************************************************************************
*     Compute Sqrt(Abs( (mu,nu|mu,nu) ) )                              *
*     Make a list of the largest such element for each shell-pair      *
*     Store in SDG.                                                    *
************************************************************************
      use ChoArr, only: iSOShl, iRS2F
      Implicit Real*8 (a-h,o-z)
      Integer nnSkal, MxBasSh
      Real*8 SDG(nnSkal)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
*                                                                      *
************************************************************************
*                                                                      *
      Call FZero(SDG(1),nnSkal)
*
      iLoc=1 ! point to 1st reduced set in index arrays
      CALL GETMEM('diagI','Allo','Real',ipDIAG,NNBSTRT(iLoc))
*
*     Read the diagonal of the integrals, (mu,nu|mu,nu)
*
      CALL CHO_IODIAG(Work(ipDIAG),2)
*
      Do jSym=1,nSym
*
         Do jRab=1,nnBstR(jSym,iLoc)
*
            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set
*
            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)
*
            iaSh = iSOShl(iag) ! shell to which it belongs
            ibSh = iSOShl(ibg)
*
            iabSh= iTri(iaSh,ibSh)
*
            SDG(iabSh)= Max(SDG(iabSh),sqrt(abs(Work(ipDIAG+kRab-1))))
*
         End Do  ! jRab loop
      End Do
*
      CALL GETMEM('diagI','Free','Real',ipDIAG,NNBSTRT(1))
*
      MxBasSh = MxOrSh
*
      Return
      End
