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
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUG2_
#define _DEBUGPRINT_

subroutine S_GEK_localisation(nIter,IterGEK,fsdim,dqdq,dq,UpMeth,SORange,nOrb2Loc,nDIIS,HDiag,useFH,CMO,nBasis,PA,nAtoms)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero,One
use Definitions, only: iwp,wp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
use Constants, only: Pi
#endif

use Definitions, only: u6
use Localisation_globals, only: OptMeth,FuncList,GradList,DispList,bias,SOFact

implicit none

integer(kind=iwp), intent(in) :: nIter,fsdim,nOrb2Loc,nBasis,nAtoms
integer(kind=iwp),intent(out) :: nDIIS
integer(kind=iwp), intent(inout) :: IterGEK
real(kind=wp),intent(in) :: Hdiag(fsdim),CMO(nBasis,nOrb2Loc),PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: dqdq,dq(fsdim)
integer(kind=iwp) :: iFirst,i,j,k,l,nExplicit,mDiis, iLast
real(kind=wp) :: gg,Cpu1,Cpu2, Tim1, Tim2, Tim3, norm,thr, dq_NR(fsdim)
real(kind=wp), allocatable :: coords(:,:),grads(:,:),Aux_1(:),Aux_2(:),e_diis(:,:),q_diis(:,:),g_diis(:,:),H_diis(:,:),&
                              dq_diis(:), FHrow_k(:)
integer(kind=iwp), parameter :: nWindow = 20, Max_IterGEK = 50, minDP = 1
real(kind=wp), External :: DDot_
character(len=6),intent(out) :: UpMeth
logical(kind=iwp), intent(in) :: SORange,useFH
character(LEN=1) :: Step_Trunc
real(kind=wp), allocatable :: CoordsAbs(:,:)

If (nOrb2Loc==0) Call abend() ! Dummy statement

call Timing(Cpu1,Tim1,Tim2,Tim3)


#ifdef _DEBUGPRINT_
write(u6,*) 'Enter S-GEK Optimizer'
#endif


! number of iterations used to build the subspace
nDIIS = min(IterGEK,nWindow)


! leave this subroutine if not enough data points were collected
if (nDIIS < mindp) then
    write(u6,*) 'Exit S-GEK Optimizer (not enough sampling points)', ndiis
    return
end if

! find the list indeces corresponding to the first and last data point
iFirst = nIter-nDIIS+1
iLast = nIter


# ifdef _DEBUG2_
    write(u6,'(A,6(I4))') "Iter,IterGEK,nDIIS,iFirst,iLast =",nIter,IterGEK,nDIIS,iFirst,iLast
    write(u6,*) "Iter    =",nIter
    write(u6,*) "IterGEK =",IterGEK
    write(u6,*) "nDIIS   =",nDIIS
    write(u6,*) "iFirst  =",iFirst
    write(u6,*) "iLast   =",iLast
# endif

! read in coordinates and gradients to build the GEK model
call mma_Allocate(coords,fsdim, nDiis,Label="coords")
call mma_Allocate(grads,fsdim, nDiis,Label="grads")
j = 0
do i=iFirst,iLast
    j = i-iFirst+1
    coords(:,j) = DispList(:,i)
    grads(:,j) = GradList(:,i)
end do

call moveref()



#ifdef _DEBUGPRINT_
    write(u6,*) 'iFirst =',iFirst
    write(u6,*) 'iLast =',iLast
    write(u6,*) 'nWindow =',nWindow
    write(u6,*) '  nDIIS =',nDIIS
    write(u6,*) 'IterGEK =',IterGEK
    call RecPrt("coords(:,:)",' ',coords,fsdim, nDiis)
    call RecPrt("grads(:,:)",' ',grads,fsdim, nDiis)
    call RecPrt("grads(:,nDiis)",' ',grads(:,nDiis),fsdim, 1)
    call RecPrt("dq(:) = NR suggestion",' ',dq,fsdim, 1)
#endif


! select subspace basis vectors; construct normalized e_diis
! -----------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FULL SPACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (OptMeth == 4) then

    ! Set up the full space
    nExplicit = fsdim
    call mma_allocate(e_diis,fsdim,nExplicit,Label='e_diis')
    e_diis(:,:) = Zero
    do k = 1,nExplicit
        e_diis(k,k) = One
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBSPACE VERSION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (OptMeth == 5) then

    if (nDIIS == 1) then
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Exit S-GEK Optimizer'
#       endif
        call mma_deallocate(grads)
        call mma_deallocate(coords)
        return
    end if

    !number of subspace basis vectors, potentially linear dependent => difference vecs of ndiis displacements and gradients +2 additional vecs (see below)
    nExplicit = 2*(nDIIS-1)+2

    call mma_allocate(e_diis,fsdim,nExplicit,Label='e_diis')
    call mma_allocate(Aux_1,fsdim,Label='Aux_1')
    call mma_allocate(Aux_2,fsdim,Label='Aux_2')

    j = 0
    thr = 1E-18_wp
    do k=1,nDIIS-1
        !n-th column of e_diis
        j = j+1
        ! gradient difference vector
        Aux_1(:) = grads(:,k+1)-grads(:,k)
        !normalize
        norm = sqrt(DDot_(fsdim,Aux_1(:),1,Aux_1(:),1))
        if (norm < thr) then
            e_diis(:,j) = Zero
        else
            e_diis(:,j) = Aux_1(:)/norm
        end if

        !(n+1)-th column of e_diis
        j = j+1
        ! displacement difference vector
        Aux_1(:) = coords(:,k+1)-coords(:,k)
        Aux_2(:) = Aux_1(:)
        !normalize
        norm = sqrt(DDot_(fsdim,Aux_2(:),1,Aux_2(:),1))
        if (norm < thr) then
            e_diis(:,j) = Zero
        else
            e_diis(:,j) = Aux_2(:)/norm
        end if


    end do
    call mma_deallocate(Aux_2)


    ! Add some unit vectors corresponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
    j = j+1
    !current gradient
    Aux_1(:) = grads(:,nDIIS)
    !normalize
    e_diis(:,j) = Aux_1(:)/sqrt(DDot_(fsdim,Aux_1(:),1,Aux_1(:),1))

    j = j+1
    !second order method's displacement suggestion
    Aux_1(:) = dq(:)
    !normalize
    e_diis(:,j) = Aux_1(:)/sqrt(DDot_(fsdim,Aux_1(:),1,Aux_1(:),1))

    call mma_deallocate(Aux_1)

#   ifdef _DEBUGPRINT_
    if (allocated(e_diis)) call RecPrt('e_diis(unorth)',' ',e_diis,fsdim,nExplicit)
#   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else
    nExplicit = 0
    write(u6,*) "ERROR: entered S_GEK_localisation even though neither GEK, nor SGEK were specified as opt method"
    call Abend()
end if !framework: fullspace/subspace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! orthogonalize e_diis; remove redundancies from linear dependences
! -----------------------------------------------------------------
do l=1,2
    j = 1
    do i=2,nExplicit
        do k=1,j
            gg = DDot_(fsdim,e_diis(:,i),1,e_diis(:,k),1)
!            if (SGEKdebug) write(u6,*) 'i,k,gg=',i,k,gg
            e_diis(:,i) = e_diis(:,i)-gg*e_diis(:,k)
        end do
        gg = DDot_(fsdim,e_diis(:,i),1,e_diis(:,i),1) ! renormalize
!        if (SGEKdebug) write(u6,*) 'j,i,gg=',j,i,gg

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
write(u6,*) '    fsdim:',fsdim
write(u6,*) 'nExplicit:',nExplicit
write(u6,*) '  IterGEK:',IterGEK
write(u6,*) '    nDIIS:',nDIIS
write(u6,*) '    mDIIS:',mDIIS

!write(u6,*) 'Check the orthonormality'
!do i=1,mDIIS
!    do j=1,i
!        write(u6,*) i,j,DDot_(fsdim,e_diis(:,i),1,e_diis(:,j),1)
!    end do
!    write(u6,*)
!end do
if (allocated(e_diis)) call RecPrt('e_diis',' ',e_diis,fsdim,nExplicit)
#endif




! Compute the projected displacement coordinates
! ----------------------------------------------
!Note that the displacements are relative to the last coordinate, coords(:,nDIIS).
! coords_diis(u) = e_diis(Kxu)^T * coords(KxK)
! where u is the subspace dimension mDiis; and K is the fullspace dimension fsdim
call mma_Allocate(q_diis,mDiis,nDiis+Max_IterGEK,Label='q_diis')
q_diis(:,:) = Zero
do i=1,nDiis ! we project only those q vectors that were used to build the subspace, so that they are fully expressed within it
    do k=1,mDiis
        q_diis(k,i) = sum( (coords(:,i)-coords(:,nDIIS)) * e_diis(:,k))
        !q_diis(k,i) = sum( (coords(:,nDIIS)-coords(:,i)) * e_diis(:,k))
    end do
end do




! Compute projected gradients
! ---------------------------
! g_diis(u) = e_diis(Kxu)^T * g(K)
call mma_Allocate(g_diis,mDiis,nDiis+Max_IterGEK,Label='g_diis')
g_diis(:,:) = Zero
do i=1,nDIIS
    do k=1,mDIIS
        g_diis(k,i) = sum(grads(:,i)*e_diis(:,k))
    end do
end do




! project also the Hessian (diagonal) onto the subspace
! -----------------------------------------------------
! H_diis(uxu) = e_diis(Kxu)^T * H(KxK) * e_diis(Kxu)
call mma_allocate(H_diis,mDIIS,mDIIS,Label='H_diis')
call mma_allocate(FHrow_k,fsdim,Label='FHrow_k')
H_diis(:,:) = Zero

! use full hessian or just Hdiag like in NR
if (useFH) then
    do k=1,fsdim
        call GetFHrow_PM(nAtoms,nOrb2Loc,PA,fsdim,FHrow_k,k,CMO,nBasis)
        do i=1,mDiis
            do j=1,mDiis
                !H_diis(i,j) = sum(e_diis(:,i)*HDiag(:)*e_diis(:,j))
                H_diis(i,j) = H_diis(i,j) + e_diis(i,k)*sum(FHrow_k(:)*e_diis(:,j))
            end do
        end do
    end do
else
    do i=1,mDiis
        do j=1,mDiis
            H_diis(i,j) = sum(e_diis(:,i)*HDiag(:)*e_diis(:,j))
        end do
    end do
end if

call mma_deallocate(FHrow_k)

! define dq as null vector in the subspace
! ----------------------------------------
call mma_allocate(dq_diis,mDiis,Label='dq_diis')
dq_diis(:) = Zero

#ifdef _DEBUGPRINT_
    call RecPrt('q_diis',' ',q_diis,mDIIS,nDIIS)
    call RecPrt('g_diis',' ',g_diis,mDIIS,nDIIS)
    call RecPrt('H_diis',' ',H_diis,mDIIS,mDIIS)
    write(u6,'(/A)') "H_diis diagonal elements:"
    do i=1,mDiis
        write(u6,'(F26.16)') H_diis(i,i)
    end do
    write(u6,*) ''
#endif




! build the surrogate model & perform the optimization
! ----------------------------------------------------
dq_NR(:) = dq(:)

if (SORANGE .and. .False.) SOFact = One * SOFact ! Dummy statement

!if (SORange) then
!  SOFact = One
!else
!  SOFact = 10000000.0_wp
!end if
!write(u6,*) "call GEK_Optimizer"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call GEK_Optimizer(mDiis,nDiis,Max_IterGEK,q_diis(:,:),g_diis(:,:),dq_diis(:),FuncList(iFirst:),H_diis(:,:),dqdq,&
                   Step_Trunc,UpMeth,SOFact,bias)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! project the resulting displacement dq_diis back into the fullspace
! ------------------------------------------------------------------
#ifdef _DEBUGPRINT_
call RecPrt('dq(:) before projecting out',' ',dq_diis(:),size(dq_diis),1)
#endif

dq(:) = Zero
do i=1,mDIIS
  dq(:) = dq(:)+dq_diis(i)*e_diis(:,i)
end do
dqdq = sqrt(DDot_(size(dq),dq(:),1,dq(:),1))



! compute angle between GEK result and NR suggestion
! --------------------------------------------------
norm = sqrt(DDot_(fsdim,dq_NR,1,dq_NR,1))
#ifdef _DEBUGPRINT_
write(u6,'(A,F12.6,2X,A,F12.3,2x,A,I4)') "Angle(dq_NR,dq) (deg) =", acos(DDot_(fsdim,dq_NR,1,dq,1)/(norm*dqdq))/Pi*180.0_wp,&
                                         "norm(dq)/norm(dq_NR) = ",dqdq/norm, "IterGEK=",IterGEK

    write(u6,*) '||dq||=',dqdq
    call RecPrt('dq(:) after projecting out',' ',dq(:),size(dq),1)
#endif



! deallocations
! -------------
call mma_Deallocate(coords)
call mma_Deallocate(grads)

call mma_Deallocate(e_diis,safe='*')
call mma_Deallocate(q_diis)
call mma_Deallocate(g_diis)
call mma_Deallocate(H_diis)
call mma_Deallocate(dq_diis)


! print timing & finalize GEK
! ---------------------------
call Timing(Cpu2,Tim1,Tim2,Tim3)

#ifdef _DEBUGPRINT_
write(u6,*) 'CPU Time for GEK iteration',Cpu2-Cpu1
write(u6,*) 'Exit S-GEK Optimizer'
#endif

contains

subroutine moveref()
    integer(kind=iwp) :: I, J

    call mma_Allocate(CoordsAbs,fsdim, nDiis,Label="CoordsAbs")
    CoordsAbs(:,:) = Zero

    ! DispList contains displacements relative to the CMO of each iteration
    !we have to switch from relative to absolute coords for the GEK model:

    CoordsAbs(:,:) = coords(:,:)
    do i = 1, ndiis
        do j = i,ndiis
            CoordsAbs(:,i) = CoordsAbs(:,i) - coords(:,j)
        end do
    end do
    coords(:,:ndiis) = CoordsAbs(:,:ndiis)

call mma_Deallocate(CoordsAbs)

end subroutine moveref



#ifdef _NOT_USED_
subroutine sizecontrol()
integer(kind=iwp) :: I, J

do i=1,fsdim
    do j=1,ndiis
        if (coords(i,j) > 0.1_wp) then
            write(u6,*) "coords(i,j) too big i,j =",i,j,coords(i,j)
        end if
    end do
end do

end subroutine sizecontrol




subroutine getUtot()
use Localisation_globals, only: UmatList
integer(kind=iwp) :: I
real(kind=wp), allocatable :: disp_summed(:) UmatProd(:,:),xUmatProd(:,:),Umat_i(:,:),kappa_summed(:,:),UmatKsum(:,:)
! compute product matrix U_1...n = U_1 * ... * U_n
! -------------------------------------------------
! look later if this can be done in pipekmezey_iter, to save comp
call mma_allocate(UmatProd,nOrb2Loc,nOrb2Loc,Label='UmatProd')
call mma_allocate(xUmatProd,nOrb2Loc,nOrb2Loc,Label='xUmatProd')
call mma_allocate(Umat_i,nOrb2Loc,nOrb2Loc,Label='Umat_i')
call mma_allocate(kappa_summed,nOrb2Loc,nOrb2Loc,Label="kappa_summed")
call mma_allocate(disp_summed,fsdim,Label="disp_summed")
call mma_allocate(UmatKsum,nOrb2Loc,nOrb2Loc,Label="UmatKsum")

xUmatProd(:,:) = Zero
Umat_i(:,:) = Zero
call unitmat(xUmatProd,nOrb2Loc)

do i=iFirst,iLast
    Umat_i(:,:) = UmatList(:,:,i)
    !call RecPrt("UmatList(:,:,i) = ",' ',UmatList(:,:,i),nOrb2Loc,nOrb2Loc)

    !call RecPrt("U_i = "," ",Umat_i,nOrb2Loc,nOrb2Loc)
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,One,xUmatProd,nOrb2Loc,Umat_i,nOrb2Loc,Zero,UmatProd,norb2Loc)

    !call RecPrt("U_1...i = "," ",UmatProd,nOrb2Loc,nOrb2Loc)

    xUmatProd(:,:) = UmatProd(:,:)
end do
!call RecPrt("U_1...n = "," ",UmatProd,nOrb2Loc,nOrb2Loc)


disp_summed(:) = Zero
UmatKsum(:,:) = Zero

do i=iFirst,iLast
    disp_summed(:) = disp_summed(:) + DispList(:,i)
end do


call vec2upper_triag(kappa_summed,nOrb2Loc,disp_summed,fsdim,.true.)
!call expkap_localisation(kappa_summed,nOrb2Loc,Umat_i,xUmatProd,UmatKsum)
call expkap_localisation(kappa_summed,nOrb2Loc,UmatKsum)

norm = sqrt(DDot_(nOrb2Loc*nOrb2Loc,UmatProd-UmatKsum,1,UmatProd-UmatKsum,1))
# ifdef _DEBUGPRINT_
call RecPrt("-K_1-K_2-...-K_n = "," ",kappa_summed,nOrb2Loc,nOrb2Loc)
call RecPrt("exp(-K_1-K_2-...-K_n) = "," ",UmatKsum,nOrb2Loc,nOrb2Loc)
call RecPrt("U_1...n - exp(-K_1-K_2-...-K_n) = "," ",UmatProd-UmatKsum,nOrb2Loc,nOrb2Loc)
write(u6,*) "norm of U_1...n - exp(-K_1-K_2-...-K_n):", norm, "ndiis =",ndiis
# endif

call mma_Deallocate(xUmatProd)
call mma_Deallocate(Umat_i)
call mma_Deallocate(UmatProd)
call mma_Deallocate(UmatKsum)
call mma_Deallocate(disp_summed)
call mma_Deallocate(kappa_summed)



end subroutine getUtot
#endif


end subroutine S_GEK_localisation
