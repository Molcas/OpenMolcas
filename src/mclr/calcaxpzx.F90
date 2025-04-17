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
!  Then call CISigma_sa, and save the results to the array WSLam
!  (wo...si... <= ...? an example of a bad variable name)
!  daxpy WSLam((K-1)*nConf1+1) to AXPzx(1,M)
!
!  End looping over K
!
!  Computing (2)
!
!  Define a density matrix that accumutates the 4*Sum half in part (2).
!  Note this density matrix as D_acc
!
!  zero out D_acc
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
!  Then call CISigma_sa, and daxpy WSLam((M-1)*nConf1+1) to
!  AXPzx(1,M)
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
!  Then daxpy coeff1 CIVec((L-1)*nConf1+1) to AXPzx(1,M)
!  and  daxpy coeff2 CIVec((K-1)*nConf1+1) to AXPzx(1,M)
!
!  End looping over M
!
!  It can be dissembled into three parts to compute the term above, as
!  labeled in the equations.
!
!  In computing part (1), note that z_MK applies to cases when K < M,
!  and z_KM applies to cases when K > M. Computing this term requires
!  calling CISigma_sa for each K != M. After the loop over M is done, there
!  will be N*(N-1) times when it is called, where N is the number of
!  states.
!
!  The rules for z_MK or z_KM are the same for part (2) as for part
!  (1). For each M, CISigma_sa is called only once, so CISigma_sa is called
!  only N times for computing part (2).
!
!  In total CISigma_sa will be called N^2 times, but the memory needed to
!  compute AXPzx is a real(kind=wp) array sized nConf1,nRoots.
!
!  In CalcAXPzx1, CISigma_sa is called N*(N-1)/2 times, but array is sized
!  N*(N-1)/2*nRoots*nConf1 (note that N = nRoots).  Maybe this is a good
!  reason to keep both algorithms.
!
!  End of the "essay".
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: ipget, W
use MCLR_Data, only: ipCI, nConf1, nDens, nNA, XISPSM
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: State_Sym, nSym, nRoots, ntAsh, nAsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: AXPzx(nConf1,nRoots)
integer(kind=iwp), intent(in) :: NPUVX, IndTUVX(ntAsh,ntAsh,ntAsh,ntAsh)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(nRoots),nnA,nnA), PUVX(NPUVX), DDg(nTri_Elem(nRoots),nTri_Elem(nRoots)), &
                             zX(nTri_Elem(nRoots-1))
integer(kind=iwp) :: I, iKK, iKL, iKL2, iKM, iKM2, iLL, iLM, ipwslam, jSym, K, L, M, nconf3, off_Ash(nSym)
real(kind=wp) :: Coeff, coeff1, coeff2, dRoots, tempda(1)
real(kind=wp), allocatable :: D_acc(:), Ddiff(:), ovrlp(:,:), Wop(:)
real(kind=wp), external :: DDot_

AXPzx(:,:) = Zero
! Converting nRoots to double prec type
dRoots = real(nRoots,8)
! Symmetry offset
off_Ash(1) = 0
do jSym=2,nSym
  off_Ash(jSym) = off_Ash(jSym-1)+nAsh(jSym-1)
end do

! Memory allocation
nConf3 = nint(max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
ipwslam = ipGet(nConf3*nRoots)

call mma_allocate(Wop,nDens)
call mma_allocate(Ddiff,nnA**2)
call mma_allocate(D_acc,nnA**2)

Wop(:) = Zero
! only a subset, and the same one, of Wop is overwritten
! so it is ok to zero it here only once

! Starting the procedure in the "essay"
! looping over M

do M=1,nRoots
  ! Computing (1)
  do K=1,nRoots
    if (K == M) cycle
    if (M > K) then
      iKM2 = nTri_Elem(M-2)+K
    else
      iKM2 = nTri_Elem(K-2)+M
    end if
    call CalcDdiff(Ddiff,GDMat,M,K,nnA,nRoots)
    Coeff = Two*zx(IKM2)
    if (K > M) Coeff = -Coeff
    call CalcWop(Wop,Ddiff,PUVX,NPUVX,IndTUVX,Coeff,off_Ash)
    call CISigma_SA(0,State_Sym,State_Sym,Wop,nDens,tempda,1,tempda,1,ipci,ipwslam,.false.)
    !call ipin(ipwslam)
    AXPzx(:,M) = AXPzx(:,M)+dRoots*W(ipwslam)%A((K-1)*nConf1+1:K*nConf1)
  end do
  ! Computing (2)
  D_acc(:) = Zero
  call CalcDacc(D_acc,GDMat,M,nnA,nRoots,zx)
  call CalcWop(Wop,D_acc,PUVX,NPUVX,IndTUVX,One,off_Ash)
  call CISigma_SA(0,State_Sym,State_Sym,Wop,nDens,tempda,1,tempda,1,ipci,ipwslam,.false.)
  AXPzx(:,M) = AXPzx(:,M)+dRoots*W(ipwslam)%A((M-1)*nConf1+1:M*nConf1)
  ! Computing (3)
  do K=2,nRoots
    IKK = nTri_Elem(K)
    IKM = iTri(K,M)
    do L=1,K-1
      ILL = nTri_Elem(L)
      IKL = iTri(K,L)
      IKL2 = nTri_Elem(K-2)+L
      ILM = iTri(L,M)
      Coeff1 = zx(IKL2)*(Two*(DDg(IKM,ILL)-DDg(IKM,IKK))+Four*DDg(IKL,ILM))
      Coeff2 = zx(IKL2)*(Two*(DDg(ILM,ILL)-DDg(ILM,IKK))-Four*DDg(IKL,IKM))

      AXPzx(:,M) = AXPzx(:,M)+Coeff1*W(ipCI)%A((L-1)*nConf1+1:L*nConf1)+Coeff2*W(ipCI)%A((K-1)*nConf1+1:K*nConf1)
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

call mma_allocate(ovrlp,nRoots,nRoots)

do M=1,nRoots
  do I=1,nRoots
    ovrlp(I,M) = ddot_(nConf1,W(ipCI)%A((I-1)*nConf1+1),1,AXPzx(:,M),1)
  end do
end do

do M=1,nRoots
  do I=1,nRoots
    AXPzx(:,M) = AXPzx(:,M)-ovrlp(I,M)*W(ipCI)%A((I-1)*nConf1+1:I*nConf1)
  end do
end do

AXPzx(:,:) = -AXPzx(:,:)

call mma_deallocate(ovrlp)

end subroutine CalcAXPzx
