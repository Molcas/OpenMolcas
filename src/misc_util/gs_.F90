!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine GS_(T,nInter,nVec,Thr)
      Implicit Real*8 (a-h,o-z)
      Real*8 T(nInter,nVec)
#include "real.fh"
!
!
      Do i = 1, nVec
!
!        Order the vectors according to the diagonal value.
!
         Call GS_Order(T(1,i),nInter,nVec-i+1)
!
!        Normalize the vector
!
         XX=Sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
!        Write (6,*) 'GS_: i,XX=',i,XX
         If (XX.gt.Thr) Then
            Call DScal_(nInter,One/XX,T(1,i),1)
         Else
            Call FZero(T(1,i),nInter)
            Go To 100
         End If
!
!        Orthogonalize against the previous vectors
!
!        |Y(new)>=|Y> - <X|Y> * |X>
!
         Do j = 1, i-1
            XY=DDot_(nInter,T(1,i),1,T(1,j),1)
!           If (Abs(XY).gt.Thr) Then
!              Write (6,*) 'GS_: j,XY=',j,XY
               Call DaXpY_(nInter,-XY,T(1,j),1,T(1,i),1)
!              XY=DDot_(nInter,T(1,i),1,T(1,j),1)
!              Write (6,*) 'GS_: j,XY=',j,XY
!           End If
         End Do
!
!        Renormalize
!
         XX=Sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
!        Write (6,*) 'GS_: i,XX=',i,XX
         If (XX.gt.Thr) Then
            Call DScal_(nInter,One/XX,T(1,i),1)
         Else
            Call FZero(T(1,i),nInter)
         End If
 100     Continue
#ifdef _DEBUGPRINT_
         Call RecPrt('GS_: T',' ',T,nInter,nInter)
#endif
      End Do
!
      Return
      End
