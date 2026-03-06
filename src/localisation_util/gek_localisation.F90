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
! Copyright (C) 2026, Lila Zapp                                        *
!                                                                      *
! Based on the S_GEK_Optimizer for SCF by R. Lindh.                    *
!***********************************************************************


subroutine GEK_localisation(nIter, Functionallist,GradientList,displacements,hdiag,fsdim,dqdq,dq,SGEKdebug,UpMeth)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero,One
use Definitions, only: iwp,wp,u6
use Localisation_globals, only: nMxIter

implicit none

integer(kind=iwp), intent(in) :: nIter,fsdim
logical, intent(in) :: SGEKdebug
real(kind=wp),intent(in) :: GradientList(fsdim,nMxIter),displacements(fsdim,nMxIter),Hdiag(fsdim)
real(kind=wp),intent(inout) :: FunctionalList(nMxIter)
real(kind=wp), intent(inout) :: dqdq,dq(fsdim)
integer(kind=iwp) :: nDiis,iFirst,i,j,k,l,nExplicit,mDiis
real(kind=wp) :: gg,Cpu1,Cpu2, Tim1, Tim2, Tim3, norm,thr
real(kind=wp), allocatable :: q(:,:),g(:,:),Aux_a(:),Aux_b(:),e_diis(:,:),q_diis(:,:),g_diis(:,:),H_diis(:,:),dq_diis(:)
integer(kind=iwp), parameter :: nWindow = 20, Max_Iter_GEK = 50
real(kind=wp), External :: DDot_
character(len=6),intent(out) :: UpMeth
logical :: SORange
character :: Step_Trunc

SORange=.false.

Functionallist(:) =-Functionallist(:)

call Timing(Cpu1,Tim1,Tim2,Tim3)

if (SGEKdebug) write(u6,*) 'Enter S-GEK Optimizer'

! Pick up coordinates and gradients in full space
! -----------------------------------------------
! number of iterations used to build the subspace
nDIIS = min(nIter,nWindow) !1 for first iteration; 2


! index of the first iteration to consider for the subspace
iFirst = nIter-nDIIS+1 !1 for first iteration; 1

call mma_Allocate(q,fsdim, nDiis,Label="q")
call mma_Allocate(g,fsdim, nDiis,Label="g")


j = 0
do i=iFirst,nIter
    j = i-iFirst+1
    !write(u6,*) 'i,j,iter=',i,j,nIter

    ! Coordinates
    q(:,j) = displacements(:,i)

    ! Gradients
    g(:,j) = GradientList(:,i)

end do

if (SGEKdebug) then
    write(u6,*) 'nWindow =',nWindow
    write(u6,*) '  nDIIS =',nDIIS
    write(u6,*) '  nIter =',nIter
    call RecPrt("g(:,:)",' ',g,fsdim, nDiis)
    call RecPrt("q(:,:)",' ',q,fsdim, nDiis)
    call RecPrt("g(:,nDiis)",' ',g(:,nDiis),fsdim, 1)
    call RecPrt("dq(:) before projecting in",' ',dq,fsdim, 1)
end if

! select subspace basis vectors; construct normalized e_diis
! Set up the full space
nExplicit = fsdim
call mma_allocate(e_diis,fsdim,nExplicit,Label='e_diis')
e_diis(:,:) = Zero
do k = 1,nExplicit
    e_diis(k,k) = One
end do

! normally mDIIS=2*nDIIS, but it can happen that not all unit vectors are linear independent (mDIIS<=2*nDIIS).
! mDIIS is then the number of linear independent e_diis column vectors that span the subspace
!mDIIS = j
mDiis = fsdim

if (SGEKdebug) then
write(u6,*) '    fsdim:',fsdim
write(u6,*) 'nExplicit:',nExplicit
write(u6,*) '    nIter:',nIter
write(u6,*) '    nDIIS:',nDIIS
write(u6,*) '    mDIIS:',mDIIS

write(u6,*) 'Check the orthonormality'
!do i=1,mDIIS
!    do j=1,i
!        write(u6,*) i,j,DDot_(fsdim,e_diis(:,i),1,e_diis(:,j),1)
!    end do
!    write(u6,*)
!end do
if (allocated(e_diis)) call RecPrt('e_diis',' ',e_diis,fsdim,nExplicit)
end if




! Compute the projected displacement coordinates
! ----------------------------------------------
!Note that the displacements are relative to the last coordinate, q(:,nDIIS).
! q_diis(u) = e_diis(Kxu)^T * q(KxK)
! where u is the subspace dimension mDiis; and K is the fullspace dimension fsdim
call mma_Allocate(q_diis,mDiis,nDiis+Max_Iter_GEK,Label='q_diis')
q_diis(:,:) = Zero
do i=1,nDiis ! we project only those q vectors that were used to build the subspace, so that they are fully expressed within it
    do k=1,mDiis
        q_diis(k,i) = sum( (q(:,i)-q(:,nDIIS)) * e_diis(:,k))
    end do
end do


! Compute projected gradients
! ---------------------------
! g_diis(u) = e_diis(Kxu)^T * g(K)
call mma_Allocate(g_diis,mDiis,nDiis+Max_Iter_GEK,Label='g_diis')
g_diis(:,:) = Zero
do i=1,nDIIS
    do k=1,mDIIS
        g_diis(k,i) = sum(g(:,i)*e_diis(:,k))
    end do
end do


! project also the Hessian (diagonal) onto the subspace
! -----------------------------------------------------
! H_diis(uxu) = e_diis(Kxu)^T * H(KxK) * e_diis(Kxu)
call mma_allocate(H_diis,mDIIS,mDIIS,Label='H_diis')

do i=1,mDiis
  do j=1,mDiis
    H_diis(i,j) = sum(e_diis(:,i)*HDiag(:)*e_diis(:,j))
  end do
end do


! define dq as null vector in the subspace
! ----------------------------------------
call mma_allocate(dq_diis,mDiis,Label='dq_diis')
dq_diis(:) = Zero

if (SGEKdebug) then
    call RecPrt('q_diis',' ',q_diis,mDIIS,nDIIS)
    call RecPrt('g_diis',' ',g_diis,mDIIS,nDIIS)
    call RecPrt('H_diis(HDiag)',' ',H_diis,mDIIS,mDIIS)
end if

! build the surrogate model & perform the optimization
! ----------------------------------------------------
Call GEK_Optimizer(mDiis,nDiis,Max_Iter_GEK,q_diis,g_diis,dq_diis,Functionallist(iFirst:),H_diis,dqdq,Step_Trunc,UpMeth,SORange)


! project the resulting displacement dq_diis back into the fullspace
! ------------------------------------------------------------------
if (SGEKdebug) call RecPrt('dq(:) before projecting out',' ',dq_diis(:),size(dq_diis),1)

dq(:) = Zero
do i=1,mDIIS
  dq(:) = dq(:)+dq_diis(i)*e_diis(:,i)
end do
dqdq = sqrt(DDot_(size(dq),dq(:),1,dq(:),1))

if (SGEKdebug) then
    write(u6,*) '||dq||=',dqdq
    call RecPrt('dq(:) after projecting out',' ',dq(:),size(dq),1)
end if


! deallocations
! -------------
call mma_Deallocate(q)
call mma_Deallocate(g)

call mma_Deallocate(e_diis,safe='*')
call mma_Deallocate(q_diis)
call mma_Deallocate(g_diis)
call mma_Deallocate(H_diis)
call mma_Deallocate(dq_diis)

! print timing & finalize GEK
! ---------------------------
if (SGEKdebug) write(u6,*) 'Exit S-GEK Optimizer'
call Timing(Cpu2,Tim1,Tim2,Tim3)

if (SGEKdebug) write(u6,*) 'CPU Time for GEK iteration',Cpu2-Cpu1

Functionallist(:) =-Functionallist(:)
end subroutine GEK_localisation
