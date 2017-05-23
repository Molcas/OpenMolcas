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
      Subroutine GS(drdq,nLambda,T,nInter,Swap,RD)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 drdq(nInter,nLambda),T(nInter,nInter)
      Logical Swap,RD
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*
*     Be careful here so that noise is not converted to a basis!
*
      Thr=1.0D-12
#ifdef _DEBUG_
      Call RecPrt('GS: dRdQ',' ',drdq,nInter,nLambda)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Initial check to see if the dRdQ vectors are independent.
*
      Call Allocate_Work(ipTemp,nInter*nLambda)
      call dcopy_(nInter*nLambda,drdq,1,Work(ipTemp),1)
      Call GS_(drdq,nInter,nLambda,Thr)
      jLambda=0
      Do i = 1, nLambda
         XX=Sqrt(DDot_(nInter,drdq(1,i),1,drdq(1,i),1))
C        Write (6,*) 'i,XX=',i,XX
         If (XX.gt.Thr) Then
            jLambda=jLambda+1
*           RD = remove degeneracies
            If (RD) Then
               If (jLambda.ne.i) Then
                  Call dCopy_(nInter,drdq(1,i),1,drdq(1,jLambda),1)
               End If
            End If
         End If
      End Do
#ifdef _DEBUG_
      Call RecPrt('GS: dRdQ(orth)',' ',drdq,nInter,nLambda)
#endif
      If ((.not.RD).and.(jLambda.ne.nLambda)) Then
         Write (6,*) ' Constraints are linear dependent!'
         Call abend()
      End If
      nLambda=jLambda
*                                                                      *
************************************************************************
*                                                                      *
*     Project away the space which is spanned by the constraints.
*
      Call FZero(T,nInter**2)
*
      call dcopy_(nInter,One,0,T,1+nInter)
#ifdef _DEBUG_
      Call RecPrt('T(orig)',' ',T,nInter,nInter)
#endif
*
*     Form 1 - P
*
      Do iLambda = 1, nLambda
         Do iInter = 1, nInter
            Do jInter = 1, nInter
               T(iInter,jInter) = T(iInter,jInter)
     &                          - drdq(iInter,iLambda)
     &                          * drdq(jInter,iLambda)
            End Do
         End Do
      End Do
#ifdef _DEBUG_
      Call RecPrt('1-P',' ',T,nInter,nInter)
#endif
*
*     Orthonormalize
*
      Call GS_(T,nInter,nInter,Thr)
*
*     Set the trailing vectors to null vectors.
*     We might have noise here.
*
      If (nLambda.ne.0) Then
         iStart = nInter-nLambda + 1
         Call FZero(T(1,iStart),nInter*nLambda)
      End If
#ifdef _DEBUG_
      Call RecPrt('1-P(GS)',' ',T,nInter,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Restore dRdQ
*
      If (.not.RD) Then
         call dcopy_(nInter*nLambda,Work(ipTemp),1,drdq,1)
      End If
      Call Free_Work(ipTemp)
*                                                                      *
************************************************************************
*                                                                      *
*     Reorder, place the null vectors at the start.
*
      j=nInter
      Do i = nInter, 1, -1
         XX=DDot_(nInter,T(1,i),1,T(1,i),1)
         If (XX.gt.Zero.and.i.ne.j) Then
            call dcopy_(nInter,T(1,i),1,T(1,j),1)
            j = j - 1
         Else If (XX.gt.Zero) Then
            j = j - 1
         End If
      End Do
*
*     Put drdq at the start
*
      call dcopy_(nInter*nLambda,drdq,1,T,1)
#ifdef _DEBUG_
      Call RecPrt('T(ReOrdered)',' ',T,nInter,nInter)
#endif
*
*     Call GS_(T,nInter,nInter,Thr)
      If (Swap) call dswap_(nInter,T(1,1),1,T(1,3),1)
*
      Return
      End
      Subroutine GS_(T,nInter,nVec,Thr)
      Implicit Real*8 (a-h,o-z)
      Real*8 T(nInter,nVec)
#include "real.fh"
*
*
      Do i = 1, nVec
*
*        Order the vectors according to the diagonal value.
*
         Call GS_Order(T(1,i),nInter,nVec-i+1)
*
*        Normalize the vector
*
         XX=Sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
C        Write (6,*) 'GS_: i,XX=',i,XX
         If (XX.gt.Thr) Then
            Call DScal_(nInter,One/XX,T(1,i),1)
         Else
            Call FZero(T(1,i),nInter)
            Go To 100
         End If
*
*        Orthogonalize against the previous vectors
*
*        |Y(new)>=|Y> - <X|Y> * |X>
*
         Do j = 1, i-1
            XY=DDot_(nInter,T(1,i),1,T(1,j),1)
C           If (Abs(XY).gt.Thr) Then
C              Write (6,*) 'GS_: j,XY=',j,XY
               Call DaXpY_(nInter,-XY,T(1,j),1,T(1,i),1)
C              XY=DDot_(nInter,T(1,i),1,T(1,j),1)
C              Write (6,*) 'GS_: j,XY=',j,XY
C           End If
         End Do
*
*        Renormalize
*
         XX=Sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
C        Write (6,*) 'GS_: i,XX=',i,XX
         If (XX.gt.Thr) Then
            Call DScal_(nInter,One/XX,T(1,i),1)
         Else
            Call FZero(T(1,i),nInter)
         End If
 100     Continue
C        Call RecPrt('GS: T',' ',T,nInter,nInter)
      End Do
*
      Return
      End
      Subroutine GS_Order(T,nInter,nVec)
      Implicit Real*8 (a-h,o-z)
      Real*8 T(nInter,nVec)
*
#ifdef _DEBUG_
      Call RecPrt('GS_Order: T(orig)','(12F6.2)',T,nInter,nVec)
#endif
      Do j = 1, nVec-1
         DiagMax=DDot_(nInter,T(1,j),1,T(1,j),1)
         iDiag=j
         Do i = j+1, nVec
            Diag=DDot_(nInter,T(1,i),1,T(1,i),1)
            If (T(i,i).gt.DiagMax) Then
               DiagMax=Diag
               iDiag=i
            End If
         End Do
         If (iDiag.ne.j) call dswap_(nInter,T(1,iDiag),1,T(1,j),1)
      End Do
#ifdef _DEBUG_
      Call RecPrt('GS_Order: T(ordered)','(12F6.2)',T,nInter,nVec)
#endif
      Return
      End
