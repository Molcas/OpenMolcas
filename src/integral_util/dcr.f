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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine DCR(Lambda,iStab1,nStab1,iStab2,nStab2,iDCR,mDCR)
      Use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (A-H,O-Z)
      Integer   iStab1(0:nStab1-1),iStab2(0:nStab2-1), iDCR(0:7)
      Integer   Index(50), Lambda_all(1275), mDCR_all(1275),
     &          iDCR_all(0:7,1275)
      Save   Index, Lambda_all, mDCR_all, iDCR_all
      Logical   Done(1275)
      Common /DCR_Done/ Done
      Common /DCR_nIndex/ nIndex
*
      Ind1=0
      Do i = 1, nStab1-1
         Do iIrrep = 1, nIrrep-1
            If (iStab1(i).eq.iOper(iIrrep)) Then
               Ind1=Ind1+2**(iIrrep-1)
               Go To 11
            End If
         End Do
 11      Continue
      End Do
      Do k = 1, nIndex
         If (Ind1.eq.Index(k)) Then
            Ind1=k
            Go To 10
         End If
      End Do
      nIndex=nIndex+1
      Index(nIndex)=Ind1
      Ind1=nIndex
10    Continue
*
      Ind2=0
      Do i = 1, nStab2-1
         Do iIrrep = 1, nIrrep-1
            If (iStab2(i).eq.iOper(iIrrep)) Then
               Ind2=Ind2+2**(iIrrep-1)
               Go To 21
            End If
         End Do
 21      Continue
      End Do
      Do k = 1, nIndex
         If (Ind2.eq.Index(k)) Then
            Ind2=k
            Go To 20
         End If
      End Do
      nIndex=nIndex+1
      Index(nIndex)=Ind2
      Ind2=nIndex
20    Continue
*
      ij = Max(Ind1,Ind2)*(Max(Ind1,Ind2)-1)/2 + Min(Ind1,Ind2)
*
      If (.Not.Done(ij)) Then
         Call DCR_(Lambda_all(ij),iStab1,nStab1,iStab2,nStab2,
     &             iDCR_all(0,ij),mDCR_all(ij))
         Done(ij)=.True.
      End If
      Lambda=Lambda_all(ij)
      mDCR  =mDCR_all(ij)
      Call ICopy(mDCR,iDCR_all(0,ij),1,iDCR,1)
*
      Return
      End
*
      Subroutine DCR_Init
      Logical   Done(1275)
      Common /DCR_Done/ Done
      Common /DCR_nIndex/ nIndex
      nindex=0
      Do I=1,1275
         Done(I)=.False.
      End Do
      Return
      End
      SubRoutine DCR_(Lambda,iStab1,nStab1,iStab2,nStab2,iDCR,mDCR)
************************************************************************
* Oject: to compute the double coset representatives (DCR) and Lambda. *
*                                                                      *
* Called from: OneEl                                                   *
*              NucAtt                                                  *
*              TwoEl                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
************************************************************************
      Use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (A-H,O-Z)
      Integer   iStab1(0:nStab1-1),iStab2(0:nStab2-1),iDCR(0:7)
      Integer   iScrt(0:7,0:7)
*
      iScrt(:,:)=0
*
*     Construct all UGV subgroups. We accumulate the number of times an
*     operator appears in each UGV subgroup.
*
      Do 20 i = 0, nIrrep-1
         Do 30 j = 0, nStab1-1
            Do 40 k = 0, nStab2-1
               jik=iEor(iStab1(j),iEor(iOper(i),iStab2(k)))
               iScrt(i,jik) = iScrt(i,jik) + 1
 40         Continue
 30      Continue
 20   Continue
*
*     Find Lambda. Lambda is the number of times an operator is in the
*     subgroup UGV. Look only in the first subgroup.
*
      Do 100 i = 0, 7
         If(iScrt(0,i).ne.0) Lambda = iScrt(0,i)
 100  Continue
*
*     Find the unique double cosets (DC) and construct the double coset
*     representatives (DCR)
*
*     Move the first double coset representative, i.e. take the first
*     operator which appears in the subgroup UGV.
*
      mDCR = 0
      Do 200 i = 0, 7
         If(iScrt(0,iOper(i)).ne.0) Then
            iDCR(mDCR)=iOper(i)
            mDCR = mDCR + 1
            Go To 201
         End If
 200  Continue
 201  Continue
*
*     Now look through the remaining subgroups UGV, if any. If a new
*     unique subgroup is found take a representative from this set.
*     Observe that the subgroups UGV are either totally identical or
*     completely disjoint.
*
      Do 210 i=1, nIrrep-1
         Do 211 k = 0, nIrrep-1
*           Check if operator exists in UGV.
            If (iScrt(i,iOper(k)).eq.0) Go To 211
            Do 212 j = 0, mDCR-1
*              See that no element of UGV is already in DCR.
               If (iDCR(j).eq.iOper(k)) Go To 210
 212        Continue
 211     Continue
*        Here if new unique subgroup UGV was found.
         Do 220 k = 0, nIrrep-1
*           Move a representative to the DCR set.
            If(iScrt(i,iOper(k)).ne.0) Then
               iDCR(mDCR)=iOper(k)
               mDCR = mDCR + 1
               Go To 221
            End If
 220     Continue
 221     Continue
*
 210  Continue
*
      Return
      End
