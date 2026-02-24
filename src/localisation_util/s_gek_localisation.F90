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

subroutine S_GEK_localisation(nIter, Functionallist,GradientList,displacements,hdiag,fsdim,dqdq,dq,SGEKdebug)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: iwp,wp,u6
use Localisation_globals, only: nMxIter

implicit none

integer(kind=iwp), intent(in) :: nIter,fsdim
logical, intent(in) :: SGEKdebug
real(kind=wp),intent(inout) :: FunctionalList(nMxIter),GradientList(fsdim,nMxIter),displacements(fsdim,nMxIter),Hdiag(fsdim)
real(kind=wp), intent(inout) :: dqdq,dq(fsdim)
integer(kind=iwp) :: nDiis,iFirst,i,j,k,l,nExplicit,mDiis
real(kind=wp) :: gg,Cpu1,Cpu2, Tim1, Tim2, Tim3
real(kind=wp), allocatable :: q(:,:),g(:,:),Aux_a(:),Aux_b(:),e_diis(:,:),q_diis(:,:),g_diis(:,:),H_diis(:,:),dq_diis(:)
integer(kind=iwp), parameter :: nWindow = 20, Max_Iter_GEK = 50
real(kind=wp), External :: DDot_
character(len=6) :: UpMeth
logical :: SORange
character :: Step_Trunc

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
    call RecPrt("dq(:)",' ',dq,fsdim, 1)
end if

! select subspace basis vectors; construct normalized e_diis
! -----------------------------------------------------------

!number of subspace basis vectors, potentially linear dependent => difference vecs of ndiis displacements and gradients +2 additional vecs (see below)
nExplicit = 2*(nDIIS-1)+2

call mma_allocate(e_diis,fsdim,nExplicit,Label='e_diis')
call mma_allocate(Aux_a,fsdim,Label='Aux_a')
call mma_allocate(Aux_b,fsdim,Label='Aux_b')

j = 0
do k=1,nDIIS-1
    !n-th column of e_diis
    j = j+1
    ! gradient difference vector
    Aux_a(:) = g(:,k+1)-g(:,k)
    !normalize
    e_diis(:,j) = Aux_a(:)/sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))

    !(n+1)-th column of e_diis
    j = j+1
    ! displacement difference vector
    Aux_a(:) = q(:,k+1)-q(:,k)
    Aux_b(:) = Aux_a(:)
    !normalize
    e_diis(:,j) = Aux_b(:)/sqrt(DDot_(fsdim,Aux_b(:),1,Aux_b(:),1))

end do
call mma_deallocate(Aux_b)

! Add some unit vectors corresponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j+1
!current gradient
Aux_a(:) = g(:,nDIIS)
!normalize
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))

j = j+1
!second order method's displacement suggestion
Aux_a(:) = dq(:)
!normalize
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))
call mma_deallocate(Aux_a)


if (SGEKdebug) then
    if (allocated(e_diis)) call RecPrt('e_diis(unorth)',' ',e_diis,fsdim,nExplicit)
end if

! orthogonalize e_diis; remove redundancies from linear dependences
! -----------------------------------------------------------------
do l=1,2
    j = 1
    do i=2,nExplicit
        do k=1,j
            gg = DDot_(fsdim,e_diis(:,i),1,e_diis(:,k),1)
            if (SGEKdebug) write(u6,*) 'i,k,gg=',i,k,gg
            e_diis(:,i) = e_diis(:,i)-gg*e_diis(:,k)
        end do
        gg = DDot_(fsdim,e_diis(:,i),1,e_diis(:,i),1) ! renormalize
        if (SGEKdebug) write(u6,*) 'j,i,gg=',j,i,gg

        if (gg > 1.0e-17_wp) then   ! Skip vector if linear dependent.
            j = j+1
            e_diis(:,j) = e_diis(:,i)/sqrt(gg)
        end if
    end do
end do

! normally mDIIS=2*nDIIS, but it can happen that not all unit vectors are linear independent (mDIIS<=2*nDIIS).
! mDIIS is then the number of linear independent e_diis column vectors that span the subspace
mDIIS = j

if (SGEKdebug) then
write(u6,*) '    fsdim:',fsdim
write(u6,*) 'nExplicit:',nExplicit
write(u6,*) '    nIter:',nIter
write(u6,*) '    nDIIS:',nDIIS
write(u6,*) '    mDIIS:',mDIIS

write(u6,*) 'Check the orthonormality'
do i=1,mDIIS
    do j=1,i
        write(u6,*) i,j,DDot_(fsdim,e_diis(:,i),1,e_diis(:,j),1)
    end do
    write(u6,*)
end do
if (allocated(e_diis)) call RecPrt('e_diis',' ',e_diis,fsdim,mDIIS)
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
dq(:) = Zero
do i=1,mDIIS
  dq(:) = dq(:)+dq_diis(i)*e_diis(:,i)
end do
dqdq = sqrt(DDot_(size(dq),dq(:),1,dq(:),1))

if (SGEKdebug) then
    write(u6,*) '||dq||=',dqdq
    call RecPrt('dq',' ',dq(:),size(dq),1)
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

end subroutine S_GEK_localisation
