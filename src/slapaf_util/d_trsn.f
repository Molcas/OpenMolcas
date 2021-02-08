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
      Function D_Trsn(Ind,iOp_,nSym)
      use Slapaf_Info, only: jStab, nStab
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer Ind(4), iOp_(4)
      Real*8 D_Trsn
*                                                                      *
************************************************************************
*                                                                      *
*     T O R S I O N S ((iAtom,R(jAtom)),T(kAtom,S(lAtom)))
*
      D_Trsn=Zero
*
      iAtom=Ind(1)
      jAtom=Ind(2)
      kAtom=Ind(3)
      lAtom=Ind(4)
      iOp_E=iOp_(1)
      iOp_R=iOp_(2)
      iOp_T=iOp_(3)
      iOp_TS=iOp_(4)
      iOp_S=iEor(iOp_T,iOp_TS)
*     Write (*,*) iAtom,jAtom,kAtom,lAtom
*     Write (*,*) iOp_E, iOp_R, iOp_T, iOp_TS, iOp_S
*
      nU_A=nStab(iAtom)
      iU_A=iU(jStab(0,iAtom),nU_A)
      nU_B=nStab(jAtom)
      iU_B=iU(jStab(0,jAtom),nU_B)
      nU_C=nStab(kAtom)
      iU_C=iU(jStab(0,kAtom),nU_C)
      nU_D=nStab(lAtom)
      iU_D=iU(jStab(0,lAtom),nU_D)
*
*     Write (*,*) ' U_A:',iU_A,nU_A
*     Write (*,*) ' U_B:',iU_B,nU_B
*     Write (*,*) ' U_C:',iU_C,nU_C
*     Write (*,*) ' U_D:',iU_D,nU_D
*
*---- Form stabilizer for (iAtom,jAtom)
*
      iOp_ER=iEor(iOp_E,iOp_R)
      If (iAtom.eq.jAtom) Then
         iU_AB=iOr(iU_A,iUR(iOp_ER,iU_B))
      Else
         iU_AB=iAnd(iU_A,iU_B)
      End If
*     Write (*,*) iAtom.eq.jAtom
*     Write (*,*) ' U_AB:',iU_AB,nU(iU_AB)
*
*-----Form stabilizer for (kAtom,lAtom)
*
      iOp_ES=iEor(iOp_E,iOp_S)
      If (kAtom.eq.lAtom) Then
         iU_CD=iOr(iU_C,iUR(iOp_ES,iU_D))
      Else
         iU_CD=iAnd(iU_C,iU_D)
      End If
*     Write (*,*) kAtom.eq.lAtom
*     Write (*,*) ' U_CD:',iU_CD,nU(iU_CD)
*
*---- Form the stabilizer for the torsion
*
      If (iAtom.ne.lAtom .or.  jAtom.ne.kAtom .or.
     &   (iAtom.eq.lAtom .and. jAtom.eq.kAtom .and.
     &     iOp_ER.ne.iOp_ES) ) Then
         iU_ABCD=iAnd(iU_AB,iU_CD)
*        Write (*,*) ' Ops!'
      Else
         iOp_ET=iEor(iOp_E,iOp_T)
         iU_ABCD=iOr(iU_AB,iUR(iOp_ET,iU_CD))
      End If
*     Write (*,*) iAtom.ne.lAtom
*     Write (*,*) jAtom.ne.kAtom
*     Write (*,*) (iAtom.eq.lAtom .and. jAtom.eq.kAtom .and.
*    &     iOp_ER.ne.iOp_ES)
      nU_ABCD=nU(iU_ABCD)
*     Write (*,*) ' U_ABCD:',iU_ABCD,nU(iU_ABCD)
*
*-----Compute the degeneracy of the torsion
*
      iDeg=nSym/nU_ABCD
      D_Trsn=DBLE(iDeg)
*
*     Write (*,*) ' D_Trsn=',D_Trsn
*
      Return
      End
