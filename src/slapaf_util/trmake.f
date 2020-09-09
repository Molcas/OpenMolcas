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
      Subroutine TRMake(TRVec,Coor,nAtoms,nTR,uMtrx,Smmtrc,nDim,dMass,
     &                  CofM)
      Implicit Real*8 (a-h,o-z)
#include "sbs.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 TRVec(6,3*nAtoms), Coor(3,nAtoms),uMtrx(3*nAtoms),
     &       CM(3), dMass(nAtoms)
      Logical SymDsp, Smmtrc(3,nAtoms), TransVar, RotVar, CofM
*
      iRout = 131
      iPrint=nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In TRMake: Coor',' ',Coor,3,nAtoms)
         Write (6,*) ' nDim=',nDim
      End If
*
      call dcopy_(6*3*nAtoms,[Zero],0,TRVec,1)
      TransVar=iAnd(iSBS,2**7).eq. 2**7
      RotVar  =iAnd(iSBS,2**8).eq. 2**8
      nTR = 0
*                                                                      *
************************************************************************
*                                                                      *
*     Translation
*
      If (TransVar) Go To 100
      Do i = 1, 3
         iCmp=2**(i-1)
         If (SymDsp(iCmp)) Then
            nTR = nTR + 1
            call dcopy_(nAtoms,[One],0,TRVec(nTR,i),18)
         End If
      End Do
 100  Continue
*                                                                      *
************************************************************************
*                                                                      *
*-----Rotation around some center
*
*     Loop over axis
*
      If (RotVar) Go To 200
      Do i = 1, 3
         CM(i)=Zero
         rNorm=Zero
         Do iAtom = 1, nAtoms
            j = (iAtom-1)*3+i
            If (CofM) Then
               rNorm=rNorm+uMtrx(j)*dMass(iAtom)
               If (Smmtrc(i,iAtom)) Then
                  CM(i)=CM(i)+uMtrx(j)*Coor(i,iAtom)*dMass(iAtom)
               End If
            Else
               rNorm=rNorm+uMtrx(j)
               If (Smmtrc(i,iAtom)) Then
                  CM(i)=CM(i)+uMtrx(j)*Coor(i,iAtom)
               End If
            End If
         End Do
         CM(i)=CM(i)/rNorm
      End Do
C     Write (6,*) 'TrMake CM=',CM
*
      Do i = 1, 3
         j = i + 1
         If (j.gt.3) j = j - 3
         k = i + 2
         If (k.gt.3) k = k - 3
*--------j and k are the index of the plane perpendicular to the axis
*
*------- Check the rotation has any mirror plane parallel to the axis
*        of the rotation. If not then the rotation will not break the
*        symmetry.
*
         iCmp = 2**(j-1) + 2**(k-1)
         If (SymDsp(iCmp)) Then
            nTR = nTR + 1
            Call DYaX (nAtoms, One,Coor(j,1),3,TRVec(nTR,k),18)
            Call DaXpY_(nAtoms,-One,CM(j)    ,0,TRVec(nTR,k),18)
            Call DYaX (nAtoms,-One,Coor(k,1),3,TRVec(nTR,j),18)
            Call DaXpY_(nAtoms, One,CM(k)    ,0,TRVec(nTR,j),18)
         End If
      End Do
 200  Continue
*                                                                      *
************************************************************************
*                                                                      *
*-----Normalize vectors
*
      Do i = 1, nTR
         rii = Zero
         Do iAtom = 1, 3*nAtoms
            rii = rii + uMtrx(iAtom)*TRVec(i,iatom)**2
         End Do
         If(rii.gt.1.d-15) Then
           Call DScal_(3*nAtoms,One/Sqrt(rii),TRVec(i,1),6)
         Else
           call dcopy_(3*nAtoms,[Zero],0,TRVec(i,1),6)
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.99) Call RecPrt(' In TRMake: TRVec',' ',
     &                           TRVec,6,3*nAtoms)
      Call TROrder(TRVec,nTR,3*nAtoms)
      If (iPrint.ge.99) Call RecPrt(' In TRMake: TRVec',' ',
     &                           TRVec,nTR,3*nAtoms)
      Call TRComp(TRVec,nTR,3*nAtoms,SmmTrc)
*
      If (iPrint.ge.99) Call RecPrt(' In TRMake: TRVec',' ',
     &                           TRVec,nTR,nDim)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine TROrder(TRVec,nTR,nX)
      Implicit Real*8 (a-h,o-z)
      Real*8 TRVec(6*nX)
      If (nTR.eq.6) Return
      Do i = 2, nX
         Do iTR = 1, nTR
            iFrom=(i-1)*6 + iTR
            iTo  =(i-1)*nTR + iTR
            TRVec(iTo)=TRVec(iFrom)
         End Do
      End Do
*
      Return
      End
      Subroutine TRComp(TRVec,nTR,nX,Smmtrc)
      Implicit Real*8 (a-h,o-z)
      Real*8 TRVec(nTR,nX)
      Logical Smmtrc(nX)
*
      If (nTR.eq.0) Return
      iDim = 0
      Do iX = 1, nX
         If (Smmtrc(iX)) Then
            iDim = iDim + 1
            If (iDim.ne.iX) call dcopy_(nTR,TRVec(1,iX),1,
     &                                 TRVec(1,iDim),1)
         End If
      End Do
*
      Return
      End
