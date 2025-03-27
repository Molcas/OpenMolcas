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
!***********************************************************************
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
!  (wo...si... <= ...? an example of a bad variable name)
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
!***********************************************************************

use ipPage, only: W
use MCLR_Data, only: nNA, nConf1, ipCI, nDens2
use MCLR_Data, only: XISPSM
use input_mclr, only: State_Sym, nSym, nRoots, ntAsh, nAsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two, Four

implicit none
! Input
real*8, dimension((nRoots-1)*nRoots/2) :: zX
real*8, dimension(nRoots*(nRoots+1)/2,nnA,nnA) :: GDMat
integer NPUVX
real*8, dimension(NPUVX) :: PUVX
integer, dimension(ntAsh,ntAsh,ntAsh,ntAsh) :: IndTUVX
real*8, dimension((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2) :: DDg
! Output
real*8, dimension(NConf1*nRoots) :: AXPzx
! Auxiliaries
real*8, dimension(:), allocatable :: Wop, Ddiff, D_acc
integer K, L, M
integer jSym
integer iKK, iLL, iKL, iKL2, iKM2, iLM, iKM
integer, dimension(nSym) :: off_Ash
real*8 coeff1, coeff2, Coeff, dRoots
integer tempi1, ipwslam, nconf3
real*8, dimension(1) :: tempda
real*8, external :: DDot_
integer, external :: ipGet
integer I
real*8, dimension(:), allocatable :: ovrlp
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
    integer iispin, iCsym, iSSym
    integer nInt1, nInt2s, nInt2a
    real*8, target :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
    integer ipCI1, ipCI2
    logical Have_2_el
  end subroutine CISigma_sa
end interface
!                                                                      *
!***********************************************************************
!                                                                      *

call FZero(AXPzx,nConf1*nRoots)
! Converting nRoots to double prec type
dRoots = real(nRoots,8)
! Symmetry offset
tempi1 = 0
do jSym=1,nSym
  off_Ash(jSym) = tempi1
  tempi1 = tempi1+nAsh(jSym)
end do

! Memory allocation
nConf3 = nint(max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
ipwslam = ipGet(nConf3*nRoots)

call mma_allocate(Wop,nDens2)
call mma_allocate(Ddiff,nnA**2)
call mma_allocate(D_acc,nnA**2)

call FZero(Wop,nDens2)
! only a subset, and the same one, of Wop is overwritten
! so it is ok to fzero it here only once

! Starting the procedure in the "essay"
! looping over M

do M=1,nRoots
  ! Computing (1)
  do K=1,nRoots
    if (K == M) cycle
    if (M > K) then
      iKM2 = (M-2)*(M-1)/2+K
    else
      iKM2 = (K-2)*(K-1)/2+M
    end if
    call CalcDdiff(Ddiff,GDMat,M,K,nnA,nRoots)
    Coeff = Two*zx(IKM2)
    if (K > M) Coeff = -Coeff
    call CalcWop(Wop,Ddiff,PUVX,NPUVX,IndTUVX,Coeff,off_Ash)
    call CISigma_SA(0,State_Sym,State_Sym,Wop,nDens2,tempda,1,tempda,1,ipci,ipwslam,.false.)
    !irc=ipin(ipwslam)
    call dAXpY_(nConf1,dRoots,W(ipwslam)%Vec((K-1)*nConf1+1),1,AXPzx((M-1)*nConf1+1),1)
  end do
  ! Computing (2)
  call FZero(D_acc,nnA**2)
  call CalcDacc(D_acc,GDMat,M,nnA,nRoots,zx)
  call CalcWop(Wop,D_acc,PUVX,NPUVX,IndTUVX,One,off_Ash)
  call CISigma_SA(0,State_Sym,State_Sym,Wop,nDens2,tempda,1,tempda,1,ipci,ipwslam,.false.)
  call dAXpY_(nConf1,dRoots,W(ipwslam)%Vec((M-1)*nConf1+1),1,AXPzx((M-1)*nConf1+1),1)
  ! Computing (3)
  do K=2,nRoots
    IKK = (K+1)*K/2
    if (K < M) IKM = (M-1)*M/2+K
    if (K >= M) IKM = (K-1)*K/2+M
    do L=1,K-1
      ILL = (L+1)*L/2
      IKL = (K-1)*K/2+L
      IKL2 = (K-1)*(K-2)/2+L
      if (L < M) ILM = (M-1)*M/2+L
      if (L >= M) ILM = (L-1)*L/2+M
      Coeff1 = zx(IKL2)*(Two*(DDg(IKM,ILL)-DDg(IKM,IKK))+Four*DDg(IKL,ILM))
      Coeff2 = zx(IKL2)*(Two*(DDg(ILM,ILL)-DDg(ILM,IKK))-Four*DDg(IKL,IKM))

      call DAXpY_(nConf1,Coeff1,W(ipCI)%Vec((L-1)*nConf1+1),1,AXPzx((M-1)*nConf1+1),1)
      call DAXpY_(nConf1,Coeff2,W(ipCI)%Vec((K-1)*nConf1+1),1,AXPzx((M-1)*nConf1+1),1)
    end do
  end do
end do

! Memory deallocation
call mma_deallocate(Wop)
call mma_deallocate(Ddiff)
call mma_deallocate(D_acc)

! Now making AXPzx orthogonal to intermediate states
! This part is not mentioned in the equation above because
! the method is straightforward.

! AXPzx_orthonal = AXPzx * (1- Sum_I |I><I|)

call mma_allocate(ovrlp,nRoots**2)

do M=1,nRoots
  do I=1,nRoots
    ovrlp((M-1)*nRoots+I) = ddot_(nConf1,W(ipCI)%Vec((I-1)*nConf1+1),1,AXPzx((M-1)*nConf1+1),1)
  end do
end do

do M=1,nRoots
  do I=1,nRoots
    call daxpy_(nConf1,-ovrlp((M-1)*nRoots+I),W(ipCI)%Vec((I-1)*nConf1+1),1,AXPzx((M-1)*nConf1+1),1)
  end do
end do

call DScal_(nRoots*nConf1,-One,AXPzx,1)

call mma_deallocate(ovrlp)

end subroutine CalcAXPzx
