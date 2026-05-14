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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _FULL_SPACE_
subroutine S_GEK_Optimizer(dq,mOV,dqdq,UpMeth,Step_Trunc,SOrange)
!***********************************************************************
!                                                                      *
!     Object: subspace gradient-enhanced kriging optimization.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry -- BMC,                 *
!             University of Uppsala, SWEDEN                            *
!             May '22, November-December '22                           *
!***********************************************************************

use InfSCF, only: Energy, HDiag, iter, IterGEK, Loosen, TimFld
use LnkLst, only: Init_LLs, LLGrad, LLx, LstPtr, SCF_V
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
#ifndef _FULL_SPACE_
use Constants, only: One
#endif
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mOV
real(kind=wp), intent(inout) :: dq(mOV)
real(kind=wp), intent(out) :: dqdq
character(len=6), intent(inout) :: UpMeth
character, intent(inout) :: Step_Trunc
logical(kind=iwp), intent(in) :: SOrange
integer(kind=iwp) :: i, iFirst, ipg, ipq, j, k, l, mDIIS, nDIIS, nExplicit
real(kind=wp) :: Cpu1, Cpu2, gg, Tim1, Tim2, Tim3
real(kind=wp), allocatable :: D(:,:), dq_diis(:), e_diis(:,:), g(:,:), g_diis(:,:), H_Diis(:,:), q(:,:), q_diis(:,:), w(:,:)
integer(kind=iwp), parameter :: Max_Iter = 50, nWindow = 20
real(kind=wp), parameter :: Beta_Disp_Min = 5.0e-3_wp, Beta_Disp_Seed = 0.05_wp, StepMax_Seed = 0.1_wp, Thr_RS = 1.0e-7_wp, &
                            ThrGrd = 1.0e-7_wp
#ifndef _FULL_SPACE_
real(kind=wp), allocatable :: aux_a(:), aux_b(:)
#endif
real(kind=wp), external :: DDot_

call Timing(Cpu1,Tim1,Tim2,Tim3)

#ifdef _DEBUGPRINT_
write(u6,*) 'Enter S-GEK Optimizer'
#endif
if (.not. Init_LLs) then
  write(u6,*) 'Link list not initiated'
  call Abend()
end if

! define first iteration considered in the subspace
! the last nDIIS iterations, of which the first is iFirst
nDIIS = min(IterGEK,nWindow)
iFirst = Iter-nDIIS+1
!if (nDIIS == 1) then
!# ifdef _DEBUGPRINT_
!  write(u6,*) 'Exit S-GEK Optimizer'
!# endif
!  return
!end if

call mma_allocate(q,mOV,nDIIS,Label='q')
call mma_allocate(g,mOV,nDIIS,Label='g')

if (nDIIS == 1) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit S-GEK Optimizer'
# endif
  call mma_deallocate(g)
  call mma_deallocate(q)
  return
end if

! Pick up coordinates and gradients in full space
j = 0
do i=iFirst,iter
  j = i-iFirst+1
  !write(u6,*) 'i,j,iter=',i,j,iter

  ! Coordinates
  ipq = LstPtr(i,LLx)
  q(:,j) = SCF_V(ipq)%A(:)

  ! Gradients
  ipg = LstPtr(i,LLGrad)
  g(:,j) = SCF_V(ipg)%A(:)

end do

#ifdef _DEBUGPRINT_
write(u6,*) 'nWindow=',nWindow
write(u6,*) 'nDIIS=',nDIIS
write(u6,*) 'IterGEK=',IterGEK
call RecPrt('q',' ',q,mOV,nDIIS)
call RecPrt('g',' ',g,mOV,nDIIS)
call RecPrt('g(:,nDIIS)',' ',g(:,nDIIS),mOV,1)
#endif

!=======================================================================
! Select the subspace

#ifdef _FULL_SPACE_

! Set up the full space
nExplicit = mOV
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')
e_diis(:,:) = Zero
do k=1,nExplicit
  e_diis(k,k) = One
end do

#else

!nExplicit = 2 * (nDIIS - 1) + mOV + 2
nExplicit = 2*(nDIIS-1)+2
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')

call mma_allocate(Aux_a,mOV,Label='Aux_a')
call mma_allocate(Aux_b,mOV,Label='Aux_b')

j = 0
do k=1,nDIIS-1
  j = j+1
  Aux_a(:) = g(:,k+1)-g(:,k)
  e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

  j = j+1
  Aux_a(:) = q(:,k+1)-q(:,k)
  !call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  Aux_b(:) = Aux_a(:)
  e_diis(:,j) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

end do
call mma_deallocate(Aux_b)

! Add some unit vectors corresponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j+1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j+1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
call mma_deallocate(Aux_a)

#endif

! Now orthogonalize all unit vectors
#ifdef _DEBUGPRINT_
if (allocated(e_diis)) call RecPrt('e_diis(unnorm)',' ',e_diis,mOV,nExplicit)
#endif
do l=1,2
  j = 1
  do i=2,nExplicit
    do k=1,j
      gg = DDot_(mOV,e_diis(:,i),1,e_diis(:,k),1)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'i,k,gg=',i,k,gg
#     endif
      e_diis(:,i) = e_diis(:,i)-gg*e_diis(:,k)
    end do
    gg = DDot_(mOV,e_diis(:,i),1,e_diis(:,i),1) ! renormalize
#   ifdef _DEBUGPRINT_
    write(u6,*) 'j,i,gg=',j,i,gg
#   endif
    if (gg > 1.0e-17_wp) then   ! Skip vector if linear dependent.
      j = j+1
      e_diis(:,j) = e_diis(:,i)/sqrt(gg)
    end if
  end do
end do
! normally mDIIS=2*nDIIS, but it can happen that not all unit vectors are linear independent (mDIIS<=2*nDIIS).
! mDIIS is then the number of linear independent e_diis column vectors that span the subspace
mDIIS = j

#ifdef _DEBUGPRINT_
write(u6,*) '      mOV:',mOV
write(u6,*) 'nExplicit:',nExplicit
write(u6,*) 'IterGEK   :',IterGEK
write(u6,*) '    nDIIS:',nDIIS
write(u6,*) '    mDIIS:',mDIIS

write(u6,*) 'Check the orthonormality'
do i=1,mDIIS
  do j=1,i
    write(u6,*) i,j,DDot_(mOV,e_diis(:,i),1,e_diis(:,j),1)
  end do
  write(u6,*)
end do
if (allocated(e_diis)) call RecPrt('e_diis',' ',e_diis,mOV,mDIIS)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the projected displacement coordinates. Note that the displacements are relative to the last coordinate, q(:,nDIIS). !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mma_allocate(q_diis,mDIIS,nDIIS+Max_Iter,Label='q_diis')
q_diis(:,:) = Zero
do i=1,nDIIS
  do k=1,mDIIS
    q_diis(k,i) = sum((q(:,i)-q(:,nDIIS))*e_diis(:,k))
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected gradients !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mma_allocate(g_diis,mDIIS,nDIIS+Max_Iter,Label='g_diis')
g_diis(:,:) = Zero
do i=1,nDIIS
  do k=1,mDIIS
    g_diis(k,i) = sum(g(:,i)*e_diis(:,k))
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('q_diis',' ',q_diis,mDIIS,nDIIS)
call RecPrt('g_diis',' ',g_diis,mDIIS,nDIIS)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Project the approximate Hessian to the subspace !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_allocate(H_diis,mDIIS,mDIIS,Label='H_diis')
!call mma_allocate(HDiag_diis,mDIIS,Label='HDiag_diis')

do i=1,mDiis
  do j=1,mDiis
    H_diis(i,j) = sum(e_diis(:,i)*HDiag(:)*e_diis(:,j))
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Undershoot avoidance: Scale along dq !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (Loosen%Factor /= One) then
  ! Components of dq in the subspace
  call mma_allocate(w,mDIIS,mDIIS,Label='w')
  gg = sqrt(DDot_(mOV,dq,1,dq,1))
  do i=1,mDIIS
    w(i,1) = DDot_(mOV,dq,1,e_diis(:,i),1)/gg
  end do
  ! D = I + (f-1) * w w^T
  ! H' = D^T H D
  call mma_allocate(D,mDIIS,mDIIS,Label='D')
  do i=1,mDIIS
    D(:,i) = (One/Loosen%Factor-One)*w(i,1)*w(:,1)
    D(i,i) = D(i,i)+One
  end do
  call dgemm_('N','N',mDIIS,mDIIS,mDIIS,One,H_diis,mDIIS,D,mDIIS,Zero,w,mDIIS)
  call dgemm_('N','N',mDIIS,mDIIS,mDIIS,One,D,mDIIS,w,mDIIS,Zero,H_diis,mDIIS)
  call mma_deallocate(D)
  call mma_deallocate(w)
end if

!do i=1,mDiis
!  HDiag_Diis(i) = H_Diis(i,i)
!end do
#ifdef _DEBUGPRINT_
call RecPrt('H_diis(HDiag)',' ',H_diis,mDIIS,mDIIS)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_allocate(dq_diis,mDiis,Label='dq_Diis')
dq_diis(:) = Zero

!=======================================================================
! Start the optimization

Call GEK_Optimizer(mDiis,nDiis,Max_Iter,q_diis,g_diis,dq_diis,Energy(iFirst:),H_diis,dqdq,Step_Trunc,UpMeth,SORange)

!=======================================================================

! Compute the displacement in the full space.
dq(:) = Zero
do i=1,mDIIS
  dq(:) = dq(:)+dq_diis(i)*e_diis(:,i)
end do
dqdq = sqrt(DDot_(size(dq),dq(:),1,dq(:),1))

#ifdef _DEBUGPRINT_
write(u6,*) '||dq||=',sqrt(DDot_(size(dq),dq(:),1,dq(:),1))
call RecPrt('dq',' ',dq(:),size(dq),1)
#endif

call Finish_Kriging()
call mma_deallocate(dq_diis)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_deallocate(H_diis)
call mma_deallocate(q_diis)
call mma_deallocate(g_diis)
call mma_deallocate(e_diis,safe='*')

call mma_deallocate(g)
call mma_deallocate(q)

#ifdef _DEBUGPRINT_
write(u6,*) 'Exit S-GEK Optimizer'
#endif
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(12) = TimFld(12)+(Cpu2-Cpu1)

end subroutine S_GEK_Optimizer
