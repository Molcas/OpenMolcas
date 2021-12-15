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
      Subroutine GenCoo(Cart,nAtom,Coor,mAtom,Vctrs,Smmtrc,
     &                  nDim,iAnr,jAnr,iTabAI,Degen)
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Cart(3,nAtom), Coor(3,mAtom), Vctrs(3*mAtom,nDim), r(3),
     &       Degen(3*nAtom)
      Integer iAnr(nAtom), jAnr(mAtom), iTabAI(2,mAtom)
      Logical New, SmmTrc(3,nAtom)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('GenCoo: Cart',' ',Cart,3,nAtom)
      Call RecPrt('GenCoo: Degen',' ',Degen,3,nAtom)
#endif
*
*-----Loop over list of symmetry unique centers
*
      iSt=1
      iDim = 0
      Do iAtom = 1, nAtom
         Fact = One/Sqrt(Degen((iAtom-1)*3 + 1))
         iEnd=iSt
         jDim=iDim
         call dcopy_(3,Cart(1,iAtom),1,Coor(1,iSt),1)
         iTabAI(1,iEnd) = iAtom
         iTabAI(2,iEnd) = iOper(0)
         jAnr(iEnd)=iAnr(iAtom)
         Do ix = 1, 3
            If (Smmtrc(ix,iAtom)) Then
               jDim=jDim+1
               call dcopy_(3*mAtom,[Zero],0,Vctrs(1,jDim),1)
               Vctrs((iEnd-1)*3+ix,jDim)=Fact
            End If
         End Do
*
*-----Loop over the operators of the point group
*
         iElem=1
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
               iElem=iElem+1
               iEnd = iEnd + 1
               Coor(1,iEnd)=x
               Coor(2,iEnd)=y
               Coor(3,iEnd)=z
               iTabAI(1,iEnd) = iAtom
               iTabAI(2,iEnd) = iOper(ig)
               jAnr(iEnd)=iAnr(iAtom)
               jDim=iDim
               Do ix = 1, 3
                  If (Smmtrc(ix,iAtom)) Then
                     jDim=jDim+1
                     Vctrs((iEnd-1)*3+ix,jDim) = r(ix)*Fact
                  End If
               End Do
            End If
         End Do      ! End loop over operators
*
         Do ix = 1, 3
            If (Smmtrc(ix,iAtom)) Then
               iDim=iDim+1
            End If
         End Do
         iSt = iEnd + 1
      End Do         ! End loop over centers
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In GenCoo: Coor',' ',Coor,3,mAtom)
      Call RecPrt(' In GenCoo: Vctrs',' ',Vctrs,3*mAtom,nDim)
      Write (6,*)
      Write (6,*) ' iTabAI'
      Write (6,*)
      Do iAtom = 1, mAtom
         Write (6,*) iTabAI(1,iAtom),iTabAI(2,iAtom)
      End Do
#endif
      Return
      End
