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
      Subroutine FixHess(H,nH,iOptC,MF,GNrm,nsAtom,AnalHess,AllowFindTS)
      use Slapaf_Parameters, only: iNeg, GNrm_Threshold, Mode
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "stdalloc.fh"
      Real*8 H(nH,nH), MF(3,nsAtom)
      Logical AnalHess, AllowFindTS, Corrected, Found
      Real*8, Allocatable:: EVal(:), LowVal(:), LowVec(:,:), Tmp(:,:),
     &                      FixVal(:), Rx(:,:), Vect(:)
#ifdef _DEBUGPRINT_
      Logical Too_Small
#endif
*
      iRout=211
      iPrint=nPrint(iRout)
*
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'AnalHess=',AnalHess
      Call RecPrt('FixHess: H(Start)',' ',H,nH,nH)
      Lu=6
      Too_Small=.False.
#endif
*
      Corrected=.False.
      HTh=1.0d-3
      ZTh=1.0D-12
      HHigh=1.0D0
*
      Call mma_allocate(EVal,nH*(nH+1)/2,Label='EVal')
*
*---- Copy elements for H
*
      SumHii=Zero
      Do i = 1, nH
         Do j = 1, i
            ij=i*(i-1)/2 + j
            EVal(ij)=H(i,j)
         End Do
         SumHii=SumHii+H(i,i)
      End Do
*
#ifdef _DEBUGPRINT_
      Write (Lu,*) 'FixHess: SumHii=',SumHii
      Call RecPrt('FixHess: Hessian',' ',H,nH,nH)
      Call TriPrt('FixHess: H',' ',EVal,nH)
#endif

*
*---- Compute eigenvalues and eigenvectors
*
*---- Davidson procedure to compute only the lowest eigenvalues
*     For small matrices we can afford to solve it directly
*      (tested with NumVal=2 in all cases)
      If (nH .le. 30) Then
        NumVal=nH
      Else
        NumVal=2
      End If
      nVStep=2
      Found=.False.
      Call mma_allocate(LowVal,NumVal,Label='LowVal')
      Call mma_allocate(LowVec,nH,NumVal,Label='LowVec')
      LowVec(:,:) = Zero
*---- Stop when the highest eigenvalue found is larger than HTh * 10
      Do While (.Not.Found)
        Call Davidson(EVal,nH,NumVal,LowVal,LowVec,iStatus)
#ifdef _DEBUGPRINT_
        Call RecPrt(' Eigenvalues',' ',LowVal,1,NumVal)
        Call RecPrt(' Eigenvectors',' ',LowVec,nH,NumVal)
#endif
        If (iStatus.gt.0) Then
          Call SysWarnMsg('FixHess',
     &      'Davidson procedure did not converge','')
        End If
        If ((LowVal(NumVal).gt.Ten*HTh).or.(NumVal.ge.nH)) Then
          Found=.True.
        Else
*----     Increase the number of eigenpairs to compute
          Call mma_allocate(Tmp,nH,NumVal,Label='Tmp')
          Tmp(:,:)=LowVec(:,:)
          Call mma_deallocate(LowVal)
          Call mma_deallocate(LowVec)
*----     At some point, start doubling the number
          If (NumVal.ge.16) NVStep=NumVal
          nVStep=Min(nVStep,nH-NumVal)
          NumVal=NumVal+nVStep
          Call mma_allocate(LowVal,NumVal,Label='LowVal')
          LowVal(:)=Zero
          Call mma_allocate(LowVec,nH,NumVal,Label='LowVec')
          LowVec(:,:)=Zero
          LowVec(:,1:NumVal-nVStep) = Tmp(:,:)
          Call mma_deallocate(Tmp)
        End If
      End Do
      Call mma_deallocate(EVal)
*
*---- Apply corrections if any ...
*
      Call mma_allocate(FixVal,NumVal,Label='FixVal')
#ifdef _DEBUGPRINT_
      Call RecPrt(' Eigenvalues',' ',LowVal,1,NumVal)
      Call RecPrt(' Eigenvectors',' ',LowVec,nH,NumVal)
#endif
      iNeg(1)=0
      jNeg=0
      rLow=Ten
      iLow=0
*     with sorted eigenvalues, jNeg=iNeg, iLow=1
      Do i = 1, NumVal
         temp=LowVal(i)
         FixVal(i)=temp
         If (temp.lt.rlow) Then
            rlow=temp
            iLow=i
         End If
*        No fixes if the Hessian is analytical
         If (.Not.AnalHess.and.(Abs(temp).lt.HTh)) Then
#ifdef _DEBUGPRINT_
            Too_Small=.True.
#endif
            Corrected=.True.
*
*           For redundant coordinates we will have some
*           eigenvalues which are zero due to the
*           redundancy.
*
            If (Abs(temp).lt.ZTh) Then
               temp=Zero
               FixVal(i)=Zero
            Else
               FixVal(i)=Sign(HTh,temp)
            End If
         End If
         If (temp.lt.Zero) Then
            iNeg(1)=iNeg(1)+1
            jNeg=i
            If ((.Not.AnalHess.or.iAnd(iOptC,256).eq.256) .and.
     &          iAnd(iOptC,128).eq.128 .and.
     &          iAnd(iOptC,4096).ne.4096) Then
*
*             Change the sign and if just too large reduce the value
*             to the default of HHigh.
*
               FixVal(i)=Min(HHigh,Abs(FixVal(i)))
               Corrected=.True.
            End If
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     If FindTS and we have got a negative eigenvalue signal that
*     we are in the TS regime.
*
      If (AllowFindTS) Then
         Call qpg_darray('TanVec',Found,nRP)
         If (iAnd(iOptC,4096).eq.4096 .and.
     &       (iNeg(1).ge.1.or.Mode.ge.0) .and.
     &       ((GNrm.le.GNrm_Threshold).or.Found)) Then
            If (iAnd(iOptC,8192).ne.8192) Then
               Write (6,*) '**************************'
               Write (6,*) '* Enable TS optimization *'
               Write (6,*) '**************************'
            End If
            iOptC=iOr(iOptC,8192)
         End If
         If (iAnd(iOptC,8192).eq.8192) Then
            Mask=2**30-1  - 2**7
            iOptC=iAnd(Mask,iOptC)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
      If (Too_Small) Then
         Write (Lu,*)
         Write (Lu,*) ' Some too small eigenvalues has been corrected'
         Write (Lu,*)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*----- M I N I M A
*                                                                      *
************************************************************************
*                                                                      *
      If (iAnd(iOptC,128).eq.128) Then
*
         If (iNeg(1).ne.0.and.iAnd(iOptC,256).ne.256) Then
            Corrected=.True.
#ifdef _DEBUGPRINT_
            Write (Lu,*) ' Some negative eigenvalues has been corrected'
            Write (Lu,*) 'iNeg=',iNeg(1)
            Write (Lu,*)
#endif
         End If
*                                                                      *
************************************************************************
*                                                                      *
*---- T R A N S I T I O N  S T A T E  S E A R C H
*                                                                      *
************************************************************************
*                                                                      *
      Else If (iAnd(iOptC,8).eq.8 .or.
     &         iAnd(iOptC,8192).eq.8192 ) Then
*                                                                      *
************************************************************************
*                                                                      *
*        1 negative eigenvalue
*
         If (iNeg(1).eq.1) Then
*
            If (Mode.le.0) Then
*
*------------- Store the eigenvector which we are following
*
               Mode=jNeg
               Call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#ifdef _DEBUGPRINT_
               Write (Lu,'(A,I3)') ' Store Original mode:',Mode
               Call RecPrt(' Reaction mode',' ',MF,3,nsAtom)
#endif
*
            Else
*
*------------- Check that it is the correct eigenvector!
*
#ifdef _DEBUGPRINT_
               Call RecPrt(' Old Reaction mode',' ',MF,3,nsAtom)
#endif
               iTest=0
               Test=Zero
               Call mma_allocate(Rx,3,nsAtom,Label='Rx')
               Do i = 1, NumVal
                  Call ReacX(LowVec(:,i),nH,Rx,3*nsAtom)
                  dRx=Sqrt(DDot_(3*nsAtom,Rx,1,Rx,1))
                  rq=Abs(DDot_(3*nsAtom,MF,1,Rx,1))/dRx
                  If (rq.gt.Test) Then
                     iTest=i
                     Test=rq
                  End If
                  Temp=FixVal(i)
#ifdef _DEBUGPRINT_
                  Write (6,*) '<old|new>,H_new=',rq,Temp
#endif
               End Do
               Call mma_deallocate(Rx)
*
*              Only iTest and jNeg may be touched
               If (iTest.eq.jNeg) Then
                  Mode = jNeg
               Else
#ifdef _DEBUGPRINT_
                  Write (Lu,*) ' Warning: wrong eigenvector has'
     &                       //' negative eigenvalue.'
#endif
*                 Keep the old vector if there is significant overlap
*                 Note: there could be a better vector not in the computed set
                  If (.Not.AnalHess.and.(Test.gt.0.50d0)) Then
                     Mode = iTest
#ifdef _DEBUGPRINT_
                     Write (Lu,*) 'Keep old eigenvector!',Mode
#endif
                     FixVal(jNeg) = Abs(FixVal(jNeg))
                     Corrected=.True.
*                 Prefer the new eigenvector if the Hessian is analytical
*                 or if the best overlap is poor
                  Else
                     Mode = jNeg
#ifdef _DEBUGPRINT_
                     Write (Lu,*) 'Take new eigenvector!',Mode
#endif
                  End if
               End If
*
               FixVal(Mode) = -Abs(FixVal(Mode))
               Call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#ifdef _DEBUGPRINT_
               Write (Lu,'(A,1X,I3)') ' Store mode:',Mode
               Call RecPrt(' New Reaction mode',' ',MF,3,nsAtom)
#endif
*
            End If
*                                                                      *
************************************************************************
*                                                                      *
*        0 negative eigenvalues
*
         Else If (iNeg(1).eq.0) Then
*
            If (Mode.lt.1) Then
*
*------------- Store the eigenvector which we are following
*
               Mode=iLow
               Call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#ifdef _DEBUGPRINT_
               Write (Lu,'(A,I3)') ' Store Original mode:',Mode
               Call RecPrt(' Reaction mode',' ',MF,3,nsAtom)
#endif
*
            Else
*
*------------- Find the eigenvector with the best overlap
*
#ifdef _DEBUGPRINT_
               Call RecPrt(' Old Reaction mode',' ',MF,3,nsAtom)
#endif
               iTest=0
               Test=Zero
               Call mma_allocate(Rx,3,nsAtom,Label='Rx')
               Do i = 1, NumVal
                  Call ReacX(LowVec(:,i),nH,Rx,3*nsAtom)
                  dRx=Sqrt(DDot_(3*nsAtom,Rx,1,Rx,1))
                  rq=Abs(DDot_(3*nsAtom,MF,1,Rx,1))/dRx
                  If (rq.gt.Test) Then
                     iTest=i
                     Test=rq
                  End If
                  Temp=FixVal(i)
#ifdef _DEBUGPRINT_
                  Write (6,*) '<old|new>,H_new=',rq,Temp
#endif
               End Do
               Call mma_deallocate(Rx)
*
*              Keep the old vector if there is significant overlap
*              Note: there could be a better vector not in the computed set
               If (Test.gt.0.50d0) Then
                  Mode = iTest
*              Prefer the lowest eigenvector if the best overlap is poor
               Else
                  Mode = iLow
#ifdef _DEBUGPRINT_
                  Write (Lu,*) ' Warning: no good overlap among'
     &                       //' the computed set of eigenvectors.'
#endif
               End if
*
               Call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#ifdef _DEBUGPRINT_
               Write (Lu,'(A,1X,I3)') ' Store mode:',Mode
               Call RecPrt(' New Reaction mode',' ',MF,3,nsAtom)
#endif
*
            End If
*
            FixVal(Mode) = - Half * Abs(FixVal(Mode))
            Corrected=.True.
#ifdef _DEBUGPRINT_
            Write (Lu,'(A,I2,A)')
     &                ' No negative eigenvalue, correction: mode ',
     &                Mode,' was changed to negative'
            Write (Lu,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        2 or more negative eigenvalues
*
         Else If (iNeg(1).ge.2) Then
*
            If (Mode.lt.1) Then
*
*------------- Store the eigenvector which we are following
*
               Mode=iLow
               Call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#ifdef _DEBUGPRINT_
               Write (Lu,'(A,I3)') ' Store Original mode:',Mode
               Call RecPrt(' Reaction mode',' ',MF,3,nsAtom)
#endif
*
            Else
*
*------------- Find the eigenvector with the best overlap
*
#ifdef _DEBUGPRINT_
              Call RecPrt(' Old Reaction mode',' ',MF,3,nsAtom)
#endif
               iTest=0
               Test=Zero
               Call mma_allocate(Rx,3,nsAtom,Label='Rx')
               Do i = 1, NumVal
                  Call ReacX(LowVec(:,i),nH,Rx,3*nsAtom)
                  dRx=Sqrt(DDot_(3*nsAtom,Rx,1,Rx,1))
                  rq=Abs(DDot_(3*nsAtom,MF,1,Rx,1))/dRx
                  If (rq.gt.Test) Then
                     iTest=i
                     Test=rq
                  End If
                  Temp=FixVal(i)
#ifdef _DEBUGPRINT_
                  Write (6,*) '<old|new>,H_new=',rq,Temp
#endif
               End Do
               Call mma_deallocate(Rx)
*
*              Keep the old vector if there is significant overlap
*              Note: there could be a better vector not in the computed set
               If (Test.gt.0.50d0) Then
                  Mode = iTest
*              Prefer the lowest eigenvector if the best overlap is poor
               Else
                  Mode = iLow
#ifdef _DEBUGPRINT_
                  Write (Lu,*) ' Warning: no good overlap among'
     &                       //' the computed set of eigenvectors.'
#endif
               End if
*
               Call ReacX(LowVec(:,Mode),nH,MF,3*nsAtom)
#ifdef _DEBUGPRINT_
               Write (Lu,'(A,1X,I3)') ' Store mode:',Mode
               Call RecPrt(' New Reaction mode',' ',MF,3,nsAtom)
#endif
*
            End If
*
*           Caution! Negative eigenvalues which are not assigned
*           to the reaction mode will be increased by two orders of magnitude
*           magnitude!
*
            Fact=Ten**2
            Do i = 1, NumVal
               Temp = FixVal(i)
               If (i.eq.Mode) Then
                  FixVal(i) = -Abs(Temp)
               Else If (Temp.lt.0) Then
                  FixVal(i) = Abs(Temp) * Fact
               Else
                  FixVal(i) = Abs(Temp)
               End If
            End Do
            Corrected=.True.
#ifdef _DEBUGPRINT_
            Write (Lu,'(A,I2,A)')
     &            ' Too many negative eigenvalue, correction: mode ',
     &            Mode,' was kept'
            Write (Lu,*)
#endif
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Else
         Write (6,*) 'No Hessian massage!'
#endif
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
      Write (Lu,*)
      Write (Lu,*)' Analysis of the Hessian'
      Write (Lu,*)
      Call RecPrt(' Eigenvalues',' ',FixVal,1,NumVal)
      Call RecPrt(' Eigenvectors',' ',LowVec,nH,NumVal)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*.... Recompute the Hessian if needed
*
      iNeg(2)=iNeg(1)
      If (Corrected) Then
*
         If (iPrint.ge.99) Then
            Call RecPrt(' Corrected eigenvalues',' ',FixVal,1,NumVal)
            Call RecPrt(' Hessian',' ',H,nH,nH)
         End If
*
         Call mma_allocate(Vect,nH,Label='Vect')
         iNeg(1)=0
         Do i=1,NumVal
           If (FixVal(i).lt.Zero) iNeg(1)=iNeg(1)+1
           Fix_Val=FixVal(i)-LowVal(i)
           If (Abs(Fix_Val).gt.1.0D-12) Then
             Fix_Val=FixVal(i)+LowVal(i)
*
*            H |i>
*
             Call dGeMV_('N',nH,nH,One,H,nH,LowVec(:,i),1,Zero,Vect,1)
*
*            H' = (I-|i><i|) H (I-|i><i|) + val_new |i><i|
*               = H - H|i> <i| - |i> <i|H + (val_old + val_new) |i> <i|
*
*  (since |i> is an eigenvector this could be used instead):
*  H' = H + (val_new - val_old) |i> <i|
*
             Do j=1,nH
               Do k=1,nH
                 H(j,k)=H(j,k)-Vect(j)*LowVec(k,i)
     &                        -Vect(k)*LowVec(j,i)
     &                        +Fix_Val*LowVec(j,i)*LowVec(k,i)
               End Do
             End Do
           End If
         End Do
         Call mma_deallocate(Vect)
*
      End If
*
      If (iPrint.ge.99) Then
         Call RecPrt('FixHess: Hessian',' ',H,nH,nH)
      End If
*
      Call mma_deallocate(FixVal)
      Call mma_deallocate(LowVal)
      Call mma_deallocate(LowVec)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
