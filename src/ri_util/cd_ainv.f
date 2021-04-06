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
      Subroutine CD_AInv(A,n,AInV,Thr_CD)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 A(n,n), AInv(n,n)
      Real*8, Allocatable :: ADiag(:), QVec(:,:)
      Integer, Allocatable :: iADiag(:)
#ifdef _ACCURACY_
      Real*8, Allocatable :: Tmp(:,:), Tmp2(:,:)
#endif
*
      Call mma_allocate(ADiag,n,Label='ADiag')
      Call mma_allocate(iADiag,n,Label='iADiag')
*
      iSeed=77
      Lu_A=IsFreeUnit(iSeed)
      Call DaName_MF_WA(Lu_A,'AMat09')
*
      iDisk=0
      Call dDaFile(Lu_A,1,A,n**2,iDisk)
*
C     Call RecPrt('A',' ',A,n,n)
*
      iSeed=iSeed+1
      Lu_Q=IsFreeUnit(iSeed)
      Call DaName_MF_WA(Lu_Q,'QMat09')
*
      call dcopy_(n,A,n+1,ADiag,1)
*
      Call CD_AInv_(n,m,ADiag,iADiag,Lu_A,Lu_Q,Thr_CD)
*
      Call mma_deallocate(ADiag)
      Call mma_deallocate(iADiag)
*
      Call mma_allocate(QVec,n,m,Label='QVec')
*
      iDisk=0
      Call dDaFile(Lu_Q,2,QVec,n*m,iDisk)
*
C     Call RecPrt('QVec','(6G20.10)',QVec,n,m)
      Call DGEMM_('N','T',n,n,m,
     &            One,QVec,n,
     &                  QVec,n,
     &            Zerp,AInv,n)
C     Call RecPrt('AInv',' ',AInv,n,n)
      Call DaEras(Lu_Q)
      Call mma_deallocate(QVec)
*                                                                      *
************************************************************************
*                                                                      *
*     Check the accuracy I-AA^1
*
#ifdef _ACCURACY_
      Call mma_allocate(Tmp,n,n,Label='Tmp')
*---
      Tmp(:,:)=Zero
*     I
      call dcopy_(n,One,0,Tmp,n+1)
*     I-AA^-1
      Call DGEMM_('N','N',n,n,n,
     &           -One,A,n,
     &                  AInv,n,
     &            One,Tmp,n)
      Call RecPrt('I-AA^-1','(6G20.12)',Tmp,n,n)
*
      Call DGEMM_('N','N',n,n,n,
     &            One,A,n,
     &                  AInv,n,
     &            Zero,Tmp,n)

      Call mma_allocate(Tmp2,n,n,Label='Tmp2')
      Tmp2(:,:)=Zero
      call dcopy_(n,One,0,Tmp2,n+1)
      Call DGEMM_('N','N',n,n,n,
     &           -One,Tmp,n,
     &                Tmp,n,
     &            One,Tmp2,n)
      Call RecPrt('I-AA^-1AA^-1','(6G20.12)',Tmp2,n,n)
*---
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine CD_AInv_(n,m,ADiag,iADiag,Lu_A,Lu_Q,Thr_CD)
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "real.fh"
      Real*8 ADiag(n)
      Integer iADiag(n)
      Logical Out_of_Core

      Real*8, Allocatable :: Scr(:), Z(:), X(:)
      Real*8, Allocatable, Target :: Qm(:), Am(:), Q_k(:), A_k(:)
      Real*8, Pointer :: Q_l(:)=>Null(), A_l(:)=>Null()
*
      nScr=3*n
      Call mma_maxDBLE(MaxMem)
      lScr=Min(MaxMem,nScr)
      Call mma_allocate(Scr,lScr,Label='Scr')
*
      nDim=n
*
      Thr=Thr_CD*1.0D-1
      Lu_Z=7
      Call DaName_MF_WA(Lu_Z,'ZMAT09')
      Call Get_Pivot_idx(ADiag,nDim,nVec,Lu_A,Lu_Z,iADiag,Scr,lScr,Thr)
      m=nVec
      If (nDim.ne.nVec) Then
         Write (6,*)
         Write (6,*) 'Detected lin. dep. in the auxiliary basis'
         Write (6,'(A,I6)')' # of aux. bfns before lin. dep. removal: ',
     &                     nDim
         Write (6,'(A,I6)')' # of aux. bfns after  lin. dep. removal: ',
     &                     nVec
      End If
*
      Call Pivot_Mat(nDim,nVec,Lu_A,Lu_Z,iADiag,Scr,lScr)
*
      Call mma_deallocate(Scr)
*
************************************************************************
*     A-vectors are now on disk. Go ahead and compute the Q-vectors!
************************************************************************
*
      ThrQ=Thr_CD*1.0D-1 ! Threshold for Inv_Cho_Factor
*
      nB=nVec
      If (nB.eq.0) Go To 777
      nQm=nB*(nB+1)/2
*
      nXZ=nB
      nQm_full= nB*(nB+1)/2
*
      Out_of_Core=2*nQm_full+5*nXZ.gt.MaxMem
*
      If (Out_Of_Core) Then
         mQm=(nQm*MaxMem-5*nXZ)/(2*nQm_full)
         a=One
         b=-Two*DBLE(mQm)
         mB=INT(-a/Two + Sqrt( (a/Two)**2 - b ))
         kQm=mB*(mB+1)/2
         If (kQm.gt.mQm) Then
            Call WarningMessage(2,'Error in CD_AInv')
            Write (6,*) 'kQm.gt.mQm!'
            Write (6,*) 'MaxMem=',MaxMem
            Write (6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
            Write (6,*) 'nB,mB=',nB,mB
            Call Abend()
         End If
      Else
         mB = nB
         kQm = nQm
      End If
*
      lQm=kQm
      lAm=lQm
*
      If (lQm.lt.1) Then
         Call WarningMessage(2,'Error in CD_AInv')
         Write (6,*) 'lQm.lt.1'
         Call Abend()
      End If
*
*     Some of memory for scratch arrays for Inv_Cho_Factor
*     Allocate memory for the A- and Q-vectors and initialize.
*
      lScr=nXZ
      Call mma_allocate(Scr,lScr,Label='Scr')
      Call mma_allocate(Qm,lQm,Label='Qm')
      Call mma_allocate(Am,lAm,Label='Am')
      Call mma_allocate(A_k,nXZ,Label='A_k')
      Call mma_allocate(Q_k,nXZ,Label='Q_k')
      Call mma_allocate(X,nXZ,Label='X')
      Call mma_allocate(Z,nXZ,Label='Z')
*
      Am(:)=Zero
      Qm(:)=Zero
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Process the A_ks to generate Q_ks.
*
      iAddr=0
      nMem=mB
      Do kCol = 1, nB
*
         iAddr_=iAddr
         If (kCol.le.nMem) Then
            iOff = (kCol-1)*kCol/2
*           Point to A_k in Am
            A_l(1:kCol) => Am(iOff+1:iOff+kCol)
            If (kCol.eq.1) Then
               nAm=nMem*(nMem+1)/2
               Call dDaFile(Lu_Z,2,Am,nAm,iAddr_)
            End If
*           Point to Q_k in Qm
            Q_l(1:kCol) => Qm(iOff+1:iOff+kCol)
         Else If (kCol.gt.nMem) Then
*           Use special scratch for A_k
            A_l(1:kCol) => A_k(1:kCol)
            Call dDaFile(Lu_Z,2,A_l,kCol,iAddr_)
*           Use special scratch for Q_k
            Q_l(1:kCol) => Q_k(1:kCol)
         End If
*
         LinDep=2
         Call Inv_Cho_Factor(A_l,kCol,
     &                       Am,Qm,nMem,
     &                       Lu_Z,Lu_Q,
     &                       Scr,lScr,
     &                       Z,X,ThrQ,
     &                       Q_l,LinDep)
*
         If (LinDep.ne.0) Then
            Call WarningMessage(2,'Error in CD_AInv')
            Write(6,*) 'Inv_Cho_Factor found linear dependence!'
            Call Abend()
         End If
*
*        Write the new A/Q-vector to file
*
         iAddr_=iAddr
         If (kCol.eq.nMem) Then
            nQm=kCol*(kCol+1)/2
            Call dDaFile(Lu_Q,1,Qm,nQm,iAddr )
            Call dDaFile(Lu_Z,1,Am,nQm,iAddr_)
         Else If (kCol.gt.nMem) Then
            Call dDaFile(Lu_Q,1,Q_l,kCol,iAddr )
            Call dDaFile(Lu_Z,1,A_l,kCol,iAddr_)
         End If
*
      End Do
*
      Q_l=>Null()
      A_l=>Null()
      Call mma_deallocate(X)
      Call mma_deallocate(Z)
      Call mma_deallocate(Q_k)
      Call mma_deallocate(A_k)
      Call mma_deallocate(Am)
      Call mma_deallocate(Qm)
      Call mma_deallocate(Scr)
 777  Continue
      Call DaEras(Lu_Z)
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Sort the Q-matrix back to the original order.
*
      Call mma_maxDBLE(MaxMem2)
*
      nBfnTot=n
      nBfn2=n**2
      lScr=Min(MaxMem2,Max(nBfn2,2*nBfnTot))
      Call mma_allocate(Scr,lScr,Label='Scr')
*
      Call Restore_Mat(nDim,nVec,Lu_Q,Lu_A,iADiag,Scr,lScr,.true.)
      Call DaEras(Lu_Q)
      Lu_Q=Lu_A
*
*     Note: after the 'Restore' call to Sort_mat, the Q-matrix is
*           no longer stored as upper-triangular but as squared
*           (zeros have been added because the corresponding argument
*            is set to .true.). The column index is still pivoted.
*
      Call mma_deallocate(Scr)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Return
      End
