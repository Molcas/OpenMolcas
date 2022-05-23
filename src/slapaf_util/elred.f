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
      Subroutine ElRed(Bmtrx,nq,nx,Gmtrx,EVal,EVec,nK,uMtrx,Scrt,g12K,
     &                 Thr)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Bmtrx(nq,nx), Gmtrx(nq,nq), EVec(nq,nq),
     &       EVal(nq*(nq+1)/2), uMtrx(nX), Scrt(nq,nX)
      Logical g12K, Diagonal
      Real*8, Parameter:: Zero_Approx=0.1D-9
      Real*8, Allocatable:: Work(:), W(:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
      Do i = 1, nq
         Do j = 1, nx
            If (Abs(Bmtrx(i,j)).lt.1.0D-10) Bmtrx(i,j)=Zero
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed: The B matrix','(5e21.12)',Bmtrx,nq,nx)
      Call RecPrt('ElRed: The u matrix','(5e21.12)',umtrx,nx,1)
#endif
      If (nq.eq.0) Then
         nK=0
         Go To 99
      End If
*                             T
*---- Form the G matrix, G=BuB
*
      Do j = 1, nX
         Do i = 1, nq
            Scrt(i,j) = BMtrx(i,j) * umtrx(j)
         End Do
      End Do
      Call DGEMM_('N','T',
     &            nq,nq,nX,
     &            1.0d0,Scrt,nq,
     &            Bmtrx,nq,
     &            0.0d0,Gmtrx,nq)
*
      Diagonal = .True.
      Do i = 1, nq
         Sum = 0.0D0
         Do j = 1, nq
            If (Abs(Gmtrx(i,j)).lt.1.0D-10) Gmtrx(i,j)=Zero
            If (j.ne.i) Sum = Sum + GMtrx(i,j)
         End Do
         Diagonal = Diagonal .and. Sum.eq.0.0D0
      End Do
*
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed: The G Matrix (nq x nq)',
     &            '(5e21.12)',Gmtrx,nq,nq)
      Write (6,*) 'Diagonal=',Diagonal
#endif
*
*---- Set up a unit matrix
*
      call dcopy_(nq*nq,[Zero],0,EVec,1)
      call dcopy_(nq,[One],0,EVec,nq+1)
*
*---- Set up the Hessian in lower triangular form, the elements
*     are symmetrized.
*
      Do i = 1, nQ
         Do j = 1, i
            ijTri = i*(i-1)/2 + j
            EVal(ijTri) = Half*(Gmtrx(i,j)+Gmtrx(j,i))
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call TriPrt('Eval prediagonalization',' ',EVal,nQ)
#endif
*
*                                                        |g 0|
*---- Compute eigenvalues and eigenvectors G(K L) = (K L)|0 0|
*     K: nonredundant vectors with eigenvalues g
*     L: redundant vectors
*
      If (.NOT.Diagonal) Then
         N=nQ
         LDZ=Max(1,N)
         Call mma_allocate(Work,3*N,Label='Work')
         Work(:)=Zero
         Call mma_allocate(W,N,Label='W')
         W(:)=Zero
         Info=0
        call dspev_('V','U',N,Eval,W,EVec,LDZ,Work,Info)
         If (Info.ne.0) Then
            Write (6,*) 'Info.ne.0'
            Write (6,*) 'Info=',Info
            Call Abend()
         End If
         Call FZero(EVal,N*(N+1)/2)
         Do i = 1, N
            ii = i*(i+1)/2
            EVal(ii)=W(i)
         End Do
         Call mma_deallocate(W)
         Call mma_deallocate(Work)
      End If
      Call DScal_(nQ*(nQ+1)/2,-1.0D0,EVal,1)
      Call JacOrd(EVal,EVec,nQ,nQ)
*     Fix standard direction.
      Do iQ = 1, nQ
         tmp=OrbPhase(EVec(1,iQ),nQ)
      End Do
      Call DScal_(nQ*(nQ+1)/2,-1.0D0,EVal,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed: Eigenvectors',' ',EVec,nQ,nQ)
      Call TriPrt('ElRed: Eigenvalues',' ',EVal,nQ)
#endif
*
*                                        -1/2
*---- Remove redundant vectors and form g     K
*
      nK = 0
      Do i = 1, nQ
         ii=i*(i+1)/2
         If (EVal(ii).gt.Thr) Then
            nK = nK + 1
         End If
         EVal(i) = EVal(ii)
c        If (g12K .and. Abs(EVal(i)).gt.Zero)
         If (g12K .and. Abs(EVal(i)).gt.Zero_Approx)
     &      Call DScal_(nQ,One/Sqrt(EVal(i)),EVec(1,i),1)
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed: The NonRedundant eigenvectors',
     &            '(5e21.12)',EVec,nQ,nK)
      Call RecPrt('ElRed: eigenvalues ','(8E12.4)',
     &            EVal,1,nK)
#endif
*
 99   Continue
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_real(tmp)
#endif
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine ElRed2(nq,nx,Gmtrx,EVal,EVec,nK,uMtrx,g12K,
     &                 Thr,BM,iBM,nB_Tot,nqB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Gmtrx(nq,nq), EVec(nq,nq),
     &       EVal(nq*(nq+1)/2), uMtrx(nX), BM(nB_Tot)
      Integer iBM(nB_Tot), nqB(nq)
      Logical g12K, Diagonal
      Real*8, Parameter:: Zero_Approx=0.1D-9
      Real*8, Allocatable:: Work(:), W(:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed2: The u matrix','(5e21.12)',umtrx,nx,1)
#endif
      If (nq.eq.0) Then
         nK=0
         Go To 99
      End If
*                             T
*---- Form the G matrix, G=BuB
*
      Call FZero(GMtrx,nq**2)
      iB = 0
      Do i = 1, nq
         nB = nqB(i)
         Do k = 1, nB
            iB = iB + 1
            ik = iBM(iB)
            B_ik=BM(iB)
*
            jB = 0
            Do j = 1, nq
               mB = nqB(j)
               Do l = 1, mB
                  jB = jB + 1
                  jl = iBM(jB)
                  If (ik.eq.jl) Then
                     B_jl=BM(jB)
                     GMtrx(i,j) = GMtrx(i,j)
     &                          + B_ik*umtrx(ik)*B_jl
                  End If
               End Do
            End Do
*
         End Do
      End Do
*
      Diagonal = .True.
      Do i = 1, nq
         Sum = 0.0D0
         Do j = 1, nq
            If (Abs(Gmtrx(i,j)).lt.1.0D-10) Gmtrx(i,j)=Zero
            If (j.ne.i) Sum = Sum + GMtrx(i,j)
         End Do
         Diagonal = Diagonal .and. Sum.eq.0.0D0
      End Do
*
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed2: The G Matrix (nq x nq)',
     &            '(5e21.12)',Gmtrx,nq,nq)
      Write (6,*) 'Diagonal=',Diagonal
#endif
*
*---- Set up a unit matrix
*
      call dcopy_(nq*nq,[Zero],0,EVec,1)
      call dcopy_(nq,[One],0,EVec,nq+1)
*
*---- Set up the Hessian in lower triangular form, the elements
*     are symmetrized.
*
      Do i = 1, nQ
         Do j = 1, i
            ijTri = i*(i-1)/2 + j
            EVal(ijTri) = Half*(Gmtrx(i,j)+Gmtrx(j,i))
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call TriPrt('Eval prediagonalization',' ',EVal,nQ)
#endif
*
*                                                        |g 0|
*---- Compute eigenvalues and eigenvectors G(K L) = (K L)|0 0|
*     K: nonredundant vectors with eigenvalues g
*     L: redundant vectors
*
      If (.NOT.Diagonal) Then
         N=nQ
         LDZ=Max(1,N)
         Call mma_allocate(Work,3*N,Label='Work')
         Work(:)=Zero
         Call mma_allocate(W,N,Label='W')
         W(:)=Zero
         Info=0
         call dspev_('V','U',N,Eval,W,EVec,LDZ,Work,Info)
         If (Info.ne.0) Then
            Write (6,*) 'Info.ne.0'
            Write (6,*) 'Info=',Info
            Call Abend()
         End If
         Call FZero(EVal,N*(N+1)/2)
         Do i = 1, N
            ii = i*(i+1)/2
            EVal(ii)=W(i)
         End Do
         Call mma_deallocate(W)
         Call mma_deallocate(Work)
      End If
      Call DScal_(nQ*(nQ+1)/2,-1.0D0,EVal,1)
      Call JacOrd(EVal,EVec,nQ,nQ)
*     Fix standard direction.
      Do iQ = 1, nQ
         tmp=OrbPhase(EVec(1,iQ),nQ)
      End Do
      Call DScal_(nQ*(nQ+1)/2,-1.0D0,EVal,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed2: Eigenvectors',' ',EVec,nQ,nQ)
      Call TriPrt('ElRed2: Eigenvalues',' ',EVal,nQ)
#endif
*
*                                        -1/2
*---- Remove redundant vectors and form g     K
*
      nK = 0
      Do i = 1, nQ
         ii=i*(i+1)/2
         If (EVal(ii).gt.Thr) Then
            nK = nK + 1
         End If
         EVal(i) = EVal(ii)
c        If (g12K .and. Abs(EVal(i)).gt.Zero)
         If (g12K .and. Abs(EVal(i)).gt.Zero_Approx)
     &      Call DScal_(nQ,One/Sqrt(EVal(i)),EVec(1,i),1)
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('ElRed2: The NonRedundant eigenvectors',
     &            '(5e21.12)',EVec,nQ,nK)
       Call RecPrt('ElRed2: eigenvalues ','(8E12.4)',
     &            EVal,1,nK)
#endif
*
 99   Continue
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_real(tmp)
#endif
      End

