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
      Function D_Bend(Ind,iOp_,nSym)
      use Slapaf_Info, only: jStab, nStab
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer Ind(3), iOp_(3)
      Real*8 D_Bend
*                                                                      *
************************************************************************
*                                                                      *
*---- B E N D  (iAtom,mAtom,jAtom)
*
      D_Bend=Zero
*
      iAtom=Ind(1)
      mAtom=Ind(2)
      jAtom=Ind(3)
      iOp_E=iOp_(1)
      iOp_R=iOp_(2)
      iOp_T=iOp_(3)
*
*     Pick up the dimension of the stabilizer and store the operator
*     indices in iU_*
*
      nU_A=nStab(iAtom)
      iU_A=iU(jStab(0,iAtom),nU_A)
*
      nU_B=nStab(mAtom)
      iU_B=iU(jStab(0,mAtom),nU_B)
*
      nU_C=nStab(jAtom)
      iU_C=iU(jStab(0,jAtom),nU_C)
*
*     Write (*,*) ' U_A',iU_A,nU_A
*     Write (*,*) ' U_B',iU_B,nU_B
*     Write (*,*) ' U_C',iU_C,nU_C
*
*---- Form the stabilizer for ((iAtom,mAtom),jAtom)
*
      If (iAtom.eq.mAtom.and.iAtom.eq.jAtom) Then
*
*        A-R(A)-T(A)
*
         iOp_ER=iEor(iOp_E,iOp_R)
         iU_AB=iOr(iU_A,iUR(iOp_ER,iU_A))
         iOp_ET=iEor(iOp_E,iOp_T)
         iU_ABC=iOr(iU_AB,iUR(iOp_ET,iU_C))
      Else If (iAtom.eq.jAtom) Then
*
*        A-R(B)-T(A)
*
         iOp_ET=iEor(iOp_E,iOp_T)
         iU_AC=iOr(iU_A,iUR(iOp_ET,iU_C))
         iU_ABC=iAnd(iU_AC,iU_B)
      Else If (mAtom.eq.jAtom) Then
*
*        A-R(B)-T(B)
*
         iU_ABC=iAnd(iU_A,iU_C)
      Else If (iAtom.eq.mAtom) Then
*
*        A-R(A)-T(C)
*
         iU_ABC=iAnd(iU_A,iU_C)
      Else
*
*        A-R(B)-T(C)
*
         iU_AB=iAnd(iU_A,iU_B)
         iU_ABC=iAnd(iU_AB,iU_C)
      End If
      nU_ABC=nU(iU_ABC)
*
*---- Compute the degeneracy of the angle
*
      ideg=nSym/nU_ABC
      D_Bend=DBLE(iDeg)
*
C     Write (*,*) ' D_Bend=',D_Bend
*
      Return
      End
