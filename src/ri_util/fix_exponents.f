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
       Subroutine Fix_Exponents(nP,mP,nC,Exp,CoeffC,CoeffP)
       Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
       Real*8, Allocatable:: Exp(:), CoeffC(:,:,:), CoeffP(:,:,:)
       Real*8, Allocatable:: Scr(:,:,:)
*
*#define _DEBUG_
#ifdef _DEBUG_
       Call RecPrt('Fix_Exponents: Exp',' ',Exp,1,nP)
       Call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(1,1,1),nP,nC)
       Call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(1,1,2),nP,nC)
       Call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(1,1,1),nP,nP)
       Call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(1,1,2),nP,nP)
#endif
*
       mP = nP
       Call Fix_Exp()
*
#ifdef _DEBUG_
       Write (6,*) 'After Fix_Exp'
       Call RecPrt('Fix_Exponents: Exp',' ',Exp,1,nP)
       Call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(1,1,1),nP,nC)
       Call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(1,1,2),nP,nC)
       Call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(1,1,1),nP,nP)
       Call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(1,1,2),nP,nP)
#endif
*
*      Reallocate arrays if the number of primitives is reduced.
*
       If (mP.ne.nP) Then
          Call mma_allocate(Scr,mP,1,1,Label='Scr')
          Scr(1:mP,1,1)=Exp(1:mP)
          Call mma_deallocate(Exp)
          Call mma_allocate(Exp,nP,Label='Exp')
          Exp(:)=Scr(:,1,1)
          Call mma_deallocate(Scr)
*
          Call mma_allocate(Scr,mP,nC,2,Label='Scr')
          Scr(1:mP,1:nC,:) = CoeffC(1:mP,1:nC,:)
          Call mma_deallocate(CoeffC)
          Call mma_allocate(CoeffC,mP,nC,2,Label='CoeffC')
          CoeffC(:,:,:)=Scr(:,:,:)
          Call mma_deallocate(Scr)
*
          Call mma_allocate(Scr,mP,mP,2,Label='Scr')
          Scr(1:mP,1:mP,:) = CoeffP(1:mP,1:mP,:)
          Call mma_deallocate(CoeffP)
          Call mma_allocate(CoeffP,mP,mP,2,Label='CoeffP')
          CoeffP(:,:,:)=Scr(:,:,:)
          Call mma_deallocate(Scr)
       End If
*
#ifdef _DEBUG_
       Write (6,*) 'After Reallocation'
       Call RecPrt('Fix_Exponents: Exp',' ',Exp,1,mP)
       Call RecPrt('Fix_Exponents: CoeffC(1)',' ',CoeffC(1,1,1),mP,nC)
       Call RecPrt('Fix_Exponents: CoeffC(2)',' ',CoeffC(1,1,2),mP,nC)
       Call RecPrt('Fix_Exponents: CoeffP(1)',' ',CoeffP(1,1,1),mP,mP)
       Call RecPrt('Fix_Exponents: CoeffP(2)',' ',CoeffP(1,1,2),mP,mP)
#endif
       Return
*                                                                      *
************************************************************************
*                                                                      *
       Contains
*                                                                      *
************************************************************************
*                                                                      *
       Subroutine Fix_Exp
*
*      First, put the exponents with all coefficients less than the
*      threshold, Thr_Skip, at the end.
*
       Thr_Skip = 1.0D-13
       Do iP = nP, 1, -1
*
          iSkip=1
          Do iC = 1, nC
             If (Abs(CoeffC(iP,iC,1)).ge.Thr_Skip) iSkip=0
          End Do
*
          If (iSkip.eq.1) Then
             If (iP.lt.mP) Then
                Temp   =Exp(iP)
                Exp(iP)=Exp(mP)
                Exp(mP)=Temp
                Do i = 1, 2
                   Temp           =CoeffP(iP,iP,i)
                   CoeffP(iP,iP,i)=CoeffP(mP,mP,i)
                   CoeffP(mP,mP,i)=Temp
                   Do iC = 1, nC
                      Temp            = CoeffC(iP,iC,i)
                      CoeffC(iP,iC,i) = CoeffC(mP,iC,i)
                      CoeffC(mP,iC,i) = Temp
                   End Do
                End Do
             End If
             mP = mP -1
          End If
*
       End Do
*
*      Second, order from largest to smallest exponent
*
       Do iP = 1, mP-1
          Do jP = iP+1, mP
             If (Exp(jP).gt.Exp(ip)) Then
                Temp   =Exp(iP)
                Exp(iP)=Exp(jP)
                Exp(jP)=Temp
                Do i = 1, 2
                   Temp           =CoeffP(iP,iP,i)
                   CoeffP(iP,iP,i)=CoeffP(jP,jP,i)
                   CoeffP(jP,jP,i)=Temp
                   Do iC = 1, nC
                      Temp            = CoeffC(iP,iC,i)
                      CoeffC(iP,iC,i) = CoeffC(jP,iC,i)
                      CoeffC(jP,iC,i) = Temp
                   End Do
                End Do
             End If
          End Do
       End Do
*
       Return
       End Subroutine Fix_Exp
*                                                                      *
************************************************************************
*                                                                      *
       End Subroutine Fix_Exponents
