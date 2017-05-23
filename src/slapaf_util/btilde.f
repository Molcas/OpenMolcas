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
      Subroutine BTilde(Coor,nAtom,nInt,B,BTILDA,A,nDim,AB,
     &                  Index,iOper,nSym,Smmtrc)
************************************************************************
*                                                                      *
*     Object: To construct the truncated B inverse matrix.             *
*                                                                      *
*             Definition |dq> = B |dx>, i.e. B(nInt,3*nAtom)           *
*                                                                      *
*             Observe that B is stored in transpose form.              *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "sbs.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 B(3,nAtom,nInt), BTILDA(3,nAtom,nInt), Coor(3,nAtom),
     &       A(3*nAtom,nDim), AB(nDim,nDim), Det(2)
      Integer   Index(3,nAtom), iOper(0:7)
      Logical SymDsp, Smmtrc(3,nAtom), TransVar, RotVar
*
      iRout = 31
      iPrint=nPrint(iRout)
      Call QEnter('BTilde')
      If (iPrint.ge.99) Call RecPrt(' In BTilde: B+',' ',B,3*nAtom,nInt)
*
*     Complement the B matrix with the vectors of the symmetry
*     preserving translations and rotations
*
      call dcopy_(3*nAtom*nDim,Zero,0,A,1)
      call dcopy_(3*nAtom*nInt,B,1,A,1)
*
      TransVar=iAnd(iSBS,2**7).eq. 2**7
      RotVar  =iAnd(iSBS,2**8).eq. 2**8
*     Write (*,*) ' Putting in symmetric translations'
*
*     Translation
*
      iAdd = nInt
      If (TransVar) Go To 101
      Do 100 i = 1, 3
*        Write (*,*) ' Translation=',i
         iCmp=2**(i-1)
         If (SymDsp(iOper,nSym,iCmp)) Then
*           Write (*,*) ' Adding translation, iCmp=',iCmp
            iAdd = iAdd + 1
*           call dcopy_(3*nAtom,Zero,0,A(1,iAdd),1)
            call dcopy_(nAtom,One,0,A(i,iAdd),3)
         End If
 100  Continue
 101  Continue
*
*-----Rotation
*     Loop over axis (yz, zx, and, xy)
*
      If (RotVar) Go To 201
      Do 200 i = 1, 3
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
         If (SymDsp(iOper,nSym,iCmp)) Then
*           Write (*,*) ' Adding rotation, iCmp=',iCmp
            iAdd = iAdd + 1
*-----------Clear the column  (iAdd)
*           call dcopy_(3*nAtom,Zero,0,A(1,iAdd),1)
            Call DYaX(nAtom,One,Coor(j,1),3,A(k,iAdd),3)
            Call DYaX(nAtom,-One,Coor(k,1),3,A(j,iAdd),3)
         End If
 200  Continue
 201  Continue
*
*---- Normalize the added vectors
*
*.... more to come
      If (iPrint.ge.99) Call RecPrt(
     &' In BTilde: B+ complemented with translation and rotation',
     &                           ' ',A,3*nAtom,nDim)
*
*     Compress the B matrix to symmetry nonredundant coordinates
*
*     Start loop over symmetry distinct centers
      ibc = 0
      Do 623 isAtom = 1, nAtom
*        Start loop over components
         Do 624 k = 1, 3
*-----------Process if not redundant, .i.e. this coordinate will
*           form a symmetrical displacement
*
            Index(k,isAtom) = -ibc
            If (Smmtrc(k,isAtom)) Then
               jsAtom = k + (isAtom-1)*3
               ibc = ibc + 1
               Index(k,isAtom) = ibc
               call dcopy_(nDim,A(jsAtom,1),3*nAtom,AB(ibc,1),nDim)
            End If
624      Continue
623   Continue
      If (iPrint.ge.99) Call RecPrt(' In BTilde: B+(compressed)',' ',
     &                              AB,nDim,nDim)
*
*     Compute the inverse
*
      nAux = 33*nDim
      iOpt=0
      Call GetMem(' Aux','Allo','Real',ipAux,nAux)
      Call DGeICD(AB,nDim,nDim,iOpt,rcond,Det,Work(ipAux),nAux)
      Call GetMem(' Aux','Free','Real',ipAux,nAux)
*     If (iPrint.ge.99) Call RecPrt(' In BTilde: (B-1)+',' ',
*    &                                AB,nDim,nDim)
      If (iPrint.ge.99) Call RecPrt(' In BTilde: (B-1)+',
     &        ' ',           AB,nDim,nDim)
*
*     Expand the B tilde matrix to symmetry district coordinates,
*     dimension is (nAtom*3 x nInt). Rows with zeros will correspond
*     to redundant coordinates
*
      Do isAtom = 1, nAtom
*        Start loop over components
         Do K = 1, 3
*           Process if not redundant
            ibc = Index(k,isAtom)
            If (ibc.gt.0) Then
               call dcopy_(nInt,AB(1,ibc),     1,BTilda(k,isAtom,1),
     &                    3*nAtom)
            Else
               call dcopy_(nInt,Zero,0,BTilda(k,isAtom,1),3*nAtom)
            End If
         End Do
      End Do
*
      If (iPrint.ge.99) Call RecPrt(' BTilda',' ',BTilda,3*nAtom,nInt)
      Call QExit('BTilde')
      Return
      End
