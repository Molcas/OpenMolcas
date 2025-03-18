!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      subroutine CalcAXPzx(AXPzx,GDMat,PUVX,NPUVX,IndTUVX,DDg,zx)
!**************************************************************************
!  Some notes from the author:
!
!  Sum_KL{z_KL * (d^2)Q/(dX_KL)(dP_MLam)}_tuvx =
!
!   2 * Sum_K[ z_MK * <K| E_tu |Lam> * (D^MM_vx - D^KK_vx) ]   !
! - 2 * Sum_K[ z_KM * <K| E_tu |Lam> * (D^MM_vx - D^KK_vx) ]   !       (1)
!
! + <M| E_tu |Lam> * 4*Sum_K[ z_MK * D^MK_vx - z_KM * D^KM_vx ]!       (2)
!
! + Sum_KL z_KL * { c_LLam * [ 2 * (DDg_KMLL - DDg_KLKK) + 4 * DDg_KLLM ]
! + c_KLam * [ 2 * (DDg_LMLL - DDg_LMKK) - 4 * DDg_KLKM ] }            (3)
!
!  Now let's consider some details of the code.
!
!  Loop over M.
!
!  Computing (1) ...
!
!  Loop over K.
!
!  For each K, calculate Ddiff^KM_tu = D^MM_tu - D^KK_tu
!
!  Then construct the operator, Wop
!   if M > K, Wop_tu =  2 * Sum_vx z_MK * Ddiff^KM_vx * g_tuvx
!   if M < K, Wop_tu = -2 * Sum_vx z_KM * Ddiff^KM_vx * g_tuvx
!
!  Then call CISigma, and save the results to the array WSLam
!  (wo...si...le...? an example of a bad variable name)
!  daxpy WSLam((K-1)*nConf1+1) to AXPzx((M-1)*nConf1+1)
!
!  End looping over K
!
!  Computing (2)
!
!  Define a density matrix that accumutates the 4*Sum half in part (2).
!  Note this density matrix as D_acc
!
!  Fzero D_acc
!
!  Loop over K again.
!   if M > K, daxpy  4*z_MK D^MK to D_acc.
!   if M < K, daxpy -4*z_KM D^MK to D_acc.
!
!  End looping over K
!
!  Now my Wop_tu is computed as
!  Wop_tu = sum_vx D_acc_vx * g_tuvx
!
!  Then call CISigma, and daxpy WSLam((M-1)*nConf1+1) to
!  AXPzx((M-1)*nConf1+1)
!
!  Computing (3)
!
!  declare coeff1 and coeff2
!
!  Loop over KL
!
!  In each KL
!  coeff1 = z_KL * [ 2 * (DDg_KMLL - DDg_KLKK) + 4 * DDg_KLLM ]
!  coeff2 = z_KL * [ 2 * (DDg_LMLL - DDg_LMKK) - 4 * DDg_KLKM ]
!
!  Then daxpy coeff1 CIVec((L-1)*nConf1+1) to AXPzx((M-1)...)
!  and  daxpy coeff2 CIVec((K-1)*nConf1+1) to AXPzx((M-1)...)
!
!  End looping over M
!
!  It can be dissembled into three parts to compute the term above, as
!  labeled in the equations.
!
!  In computing part (1), note that z_MK applies to cases when K < M,
!  and z_KM applies to cases when K > M. Computing this term requires
!  calling CISigma for each K != M. After the loop over M is done, there
!  will be N*(N-1) times when it is called, where N is the number of
!  states.
!
!  The rules for z_MK or z_KM are the same for part (2) as for part
!  (1). For each M, CISigma iscalled only once, so CISigma is called
!  only N times for computing part (2).
!
!  In total CISigma will be called N^2 times, but the memory needed to
!  compute AXPzx is a real*8 array sized nRoots*nConf1.
!
!  In CalcAXPzx1, CISigma is called N*(N-1)/2 times, but array is sized
!  N*(N-1)/2*nRoots*nConf1 (note that N = nRoots).  Maybe this is a good
!  reason to keep both algorithms.
!
!  End of the "essay".
!
!***********************************************************************
      use ipPage, only: W
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: nNA, nConf1, ipCI, nDens2
      use MCLR_Data, only: XISPSM
      use input_mclr, only: State_Sym,nSym,nRoots,ntAsh,nAsh
      Implicit None

!*****Input
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
      Integer NPUVX
      Real*8,DIMENSION(NPUVX)::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::DDg
!*****Output
      Real*8,DIMENSION(NConf1*nRoots)::AXPzx
!*****Auxiliaries
      Real*8,DIMENSION(:),Allocatable::Wop,Ddiff,D_acc
      INTEGER K,L,M
      INTEGER jSym
      INTEGER iKK,iLL,iKL,iKL2,iKM2,iLM,iKM
      Integer,DIMENSION(nSym):: off_Ash
      Real*8 coeff1,coeff2,Coeff,dRoots
      Integer tempi1,ipwslam,nconf3
      Real*8,DIMENSION(1)::tempda
      Real*8, External:: DDot_
      Integer, External:: ipGet

      INTEGER I
      Real*8,DIMENSION(:),Allocatable:: ovrlp
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface
       SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,&
     &                       Int2a,nInt2a,ipCI1,ipCI2, Have_2_el)
       Integer iispin, iCsym, iSSym
       Integer nInt1, nInt2s, nInt2a
       Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
       Integer ipCI1, ipCI2
       Logical Have_2_el
       End SubRoutine CISigma_sa
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *


      CALL FZero(AXPzx,nConf1*nRoots)
!*****Covnerting nRoots to double prec type
      dRoots=Real(nRoots,8)
!*****Symmetry off-set
      tempi1=0
      DO jSym=1, nSym
       off_Ash(jSym)=tempi1
       tempi1=tempi1+nAsh(jSym)
      END DO

!*****memory allocation
      nConf3=nint(Max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
      ipwslam=ipGet(nConf3*nRoots)

      CALL mma_allocate(Wop,nDens2)
      CALL mma_allocate(Ddiff,nnA**2)
      CALL mma_allocate(D_acc,nnA**2)

      CALL FZero(Wop,nDens2)
!     only a subset, and the same one, of Wop is overwritten
!     so it is ok to fzero it here only once

!**** Starting the procedure in the "essay"
!***  looping over M


      DO M=1,nRoots
!***  Computing (1)
       Do K=1,nRoots
        IF(K.eq.M) Cycle
        IF(M.gt.K) THEN
         iKM2=(M-2)*(M-1)/2+K
        ELSE
         iKM2=(K-2)*(K-1)/2+M
        END IF
        CALL CalcDdiff(Ddiff,GDMat,M,K,nnA,nRoots)
        Coeff=2.0d0*zx(IKM2)
        IF(K.gt.M) Coeff=-Coeff
        CALL CalcWop(Wop,Ddiff,PUVX,NPUVX,IndTUVX,Coeff,off_Ash)
        CALL CISigma_SA(0,State_Sym,State_Sym,Wop,nDens2,tempda,1,      &
     &  tempda,1,ipci,ipwslam,.false.)
!        irc=ipin(ipwslam)
        CALL dAXpY_(nConf1,dRoots,W(ipwslam)%Vec((K-1)*nConf1+1),1,     &
     &  AXPzx((M-1)*nConf1+1),1)
       End Do
!***  Computing (2)
       CALL FZero(D_acc,nnA**2)
       CALL CalcDacc(D_acc,GDMat,M,nnA,nRoots,zx)
       CALL CalcWop(Wop,D_acc,PUVX,NPUVX,IndTUVX,1.0d0,off_Ash)
       CALL CISigma_SA(0,State_Sym,State_Sym,Wop,nDens2,tempda,1,       &
     & tempda,1,ipci,ipwslam,.false.)
       CALL dAXpY_(nConf1,dRoots,W(ipwslam)%Vec((M-1)*nConf1+1),1,      &
     & AXPzx((M-1)*nConf1+1),1)
!***  Computing (3)
       Do K=2,nRoots
        IKK=(K+1)*K/2
        IF(K.lt.M) IKM=(M-1)*M/2+K
        IF(K.ge.M) IKM=(K-1)*K/2+M
        do L=1,K-1
         ILL=(L+1)*L/2
         IKL=(K-1)*K/2+L
         IKL2=(K-1)*(K-2)/2+L
         IF(L.lt.M) ILM=(M-1)*M/2+L
         IF(L.ge.M) ILM=(L-1)*L/2+M
         Coeff1=zx(IKL2)*                                               &
     &   (2.0d0*(DDg(IKM,ILL)-DDg(IKM,IKK))+4.0d0*DDg(IKL,ILM))
         Coeff2=zx(IKL2)*                                               &
     &   (2.0d0*(DDg(ILM,ILL)-DDg(ILM,IKK))-4.0d0*DDg(IKL,IKM))

         CALL DAXpY_(nConf1,Coeff1,W(ipCI)%Vec((L-1)*nConf1+1),1,       &
     &                                  AXPzx((M-1)*nConf1+1),1)
         CALL DAXpY_(nConf1,Coeff2,W(ipCI)%Vec((K-1)*nConf1+1),1,       &
     &                                  AXPzx((M-1)*nConf1+1),1)
        end do
       End Do
      END DO

!**** memory deallocation
      CALL mma_deallocate(Wop)
      CALL mma_deallocate(Ddiff)
      CALL mma_deallocate(D_acc)

!**** Now making AXPzx orthogonal to intermediate states
!**** This part is not mentioned in the equation above because
!**** the method is straightforward.

!**** AXPzx_orthonal = AXPzx * (1- Sum_I |I><I|)

      CALL mma_allocate(ovrlp,nRoots**2)

      DO M=1,nRoots
       Do I=1,nRoots
        ovrlp((M-1)*nRoots+I)=                                          &
     & ddot_(nConf1,W(ipCI)%Vec((I-1)*nConf1+1),1,                      &
     &                   AXPzx((M-1)*nConf1+1),1)
       End Do
      END DO

      DO M=1,nRoots
       Do I=1,nRoots
        CALL daxpy_(nConf1,-ovrlp((M-1)*nRoots+I),                      &
     &                W(ipCI)%Vec((I-1)*nConf1+1),1,                    &
     &                     AXPzx((M-1)*nConf1+1),1)
       End Do
      END DO

      CALL DScal_(nRoots*nConf1,-1.0d0,AXPzx,1)

      CALL mma_deallocate(ovrlp)

      END SUBROUTINE CalcAXPzx
