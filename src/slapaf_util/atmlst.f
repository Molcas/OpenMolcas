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
      Subroutine AtmLst(Cart,nAtom,Coor,mAtom)
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Cart(3,nAtom), Coor(3,mAtom), r(3)
      Logical New
*
*     Call RecPrt(' In AtmLst:Cart',' ',Cart,3,nAtom)
*
*-----Loop over list of symmetry unique centers
*
      iSt=1
      Do iAtom = 1, nAtom
         iEnd=iSt
         call dcopy_(3,Cart(1,iAtom),1,Coor(1,iSt),1)
*
*-----Loop over the operators of the point group
*
         Do ig = 1, nIrrep-1
            r(1)=One
            If (iAnd(iOper(ig),1).ne.0) r(1)=-One
            r(2)=One
            If (iAnd(iOper(ig),2).ne.0) r(2)=-One
            r(3)=One
            If (iAnd(iOper(ig),4).ne.0) r(3)=-One
            x=r(1)*Cart(1,iAtom)
            y=r(2)*Cart(2,iAtom)
            z=r(3)*Cart(3,iAtom)
*
            New=.True.
            Do iGo = iSt, iEnd
               If (New .and. x.eq.Coor(1,iGo)
     &                 .and. y.eq.Coor(2,iGo)
     &                 .and. z.eq.Coor(3,iGo)) New=.False.
            End Do
            If (New) Then
               iEnd = iEnd + 1
               Coor(1,iEnd)=x
               Coor(2,iEnd)=y
               Coor(3,iEnd)=z
            End If
         End Do      ! End loop over operators
         iSt = iEnd + 1
      End Do         ! End loop over centers
*
*     Call RecPrt(' In AtmLst: Coor',' ',Coor,3,mAtom)
      Return
      End
