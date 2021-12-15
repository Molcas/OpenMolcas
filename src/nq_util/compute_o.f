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
      Subroutine Compute_O(ZA,RA,nAtoms,Z_Tot,T,O,Lambda)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 ZA(nAtoms), RA(3,nAtoms), T(3), O(3,3), Lambda(3)
*     Local arrays
      Real*8 EVal(6), M(3,3)
*                                                                      *
************************************************************************
*                                                                      *
*---- Form the nuclear charge moment tensor
*
      Call Compute_M(ZA,nAtoms,RA,Z_Tot,T,M)
*                                                                      *
************************************************************************
*                                                                      *
*---- Diagonalize the nuclear charge momentum tensor to get
*     the principle axis system.
*
      Call FZero(O,9)
      call dcopy_(3,[One],0,O,4)
      EVal(1)=M(1,1)
      EVal(2)=M(2,1)
      EVal(3)=M(2,2)
      EVal(4)=M(3,1)
      EVal(5)=M(3,2)
      EVal(6)=M(3,3)
      Call Jacob(EVal,O,3,3)
C     Call JacOrd(EVal,O,3,3)
#ifdef _DEBUGPRINT_
      Call TriPrt('RotGrd: EVal',' ',EVal,3)
      Call RecPrt('RotGrd: O',' ',O,3,3)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Lambda(1)=EVal(1)
      Lambda(2)=EVal(3)
      Lambda(3)=EVal(6)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
