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
      Function D_Bond(Ind,iOp_,nSym)
      use Slapaf_Info, only: jStab, nStab
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer Ind(2), iOp_(2)
      Real*8 D_Bond
*                                                                      *
************************************************************************
*                                                                      *
*---- B O N D  (iAtom,jAtom)
*
*---- Compute the stabilizer of P(A) & R(B), this is done in two ways.
*
*     P(A)-R(B) = A-PR(B)
*
*     A=/=B, the stabilizer is formed as the intersection of
*            the stabilizers of A and B.
*
*     A=B, the stabilizer is formed as union of U and R(U)
*
      D_Bond=Zero
*
      iAtom=Ind(1)
      jAtom=Ind(2)
      iOp_E=iOp_(1)
      iOp_R=iOp_(2)
*
      nU_A=nStab(iAtom)
      iU_A=iU(jStab(0,iAtom),nU_A)
      nU_B=nStab(jAtom)
      iU_B=iU(jStab(0,jAtom),nU_B)
*
      If (iAtom.eq.jAtom) Then
         iOp_ER=iEor(iOp_E,iOp_R)
         iU_AB=iOr(iU_A,iUR(iOp_ER,iU_A))
      Else
         iU_AB=iAnd(iU_A,iU_B)
      End If
      nU_AB=nU(iU_AB)
*
*-----Now evaluate the degeneracy of the bond.
*
      iDeg=nSym/nU_AB
      D_Bond=DBLE(iDeg)
*
*     Write (*,*) ' D_Bond=',D_Bond
*
      Return
      End
