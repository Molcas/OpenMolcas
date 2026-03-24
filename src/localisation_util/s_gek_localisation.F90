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

!#define _DEBUGPRINT_

subroutine S_GEK_localisation(Iter_GEK,hdiag,fsdim,dqdq,dq,UpMeth,SORange,nOrb2Loc,usmitigation)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero,One,Pi
use Definitions, only: iwp,wp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif
use Localisation_globals, only: Loosen,OptMeth,FuncList,GradList,DispList,UmatList
use Definitions, only: u6

implicit none

integer(kind=iwp), intent(in) :: Iter_GEK,fsdim,nOrb2Loc
real(kind=wp),intent(in) :: Hdiag(fsdim)
real(kind=wp), intent(inout) :: dqdq,dq(fsdim)
integer(kind=iwp) :: nDiis,iFirst,i,j,k,l,nExplicit=0,mDiis
real(kind=wp) :: gg,Cpu1,Cpu2, Tim1, Tim2, Tim3, norm,thr, SOFact
real(kind=wp), allocatable :: q(:,:),g(:,:),Aux_a(:),Aux_b(:),e_diis(:,:),q_diis(:,:),g_diis(:,:),H_diis(:,:),dq_diis(:),&
                              w(:,:),D(:,:),dq_NR(:),UmatProd(:,:),xUmatProd(:,:),Umat_i(:,:),disp_summed(:),kappa_summed(:,:),&
                              UmatKsum(:,:)
integer(kind=iwp), parameter :: nWindow =22, Max_Iter_GEK = 50
real(kind=wp), External :: DDot_
character(len=6),intent(out) :: UpMeth
logical, intent(in) :: SORange,usmitigation
character :: Step_Trunc

call Timing(Cpu1,Tim1,Tim2,Tim3)

call mma_allocate(dq_NR,fsdim,Label='dq_NR')
dq_NR(:) = dq(:)


#ifdef _DEBUGPRINT_
write(u6,*) 'Enter S-GEK Optimizer'
#endif

! number of iterations used to build the subspace
if (Iter_GEK > 3) then
    nDIIS = min(Iter_GEK,nWindow)-2 !skip the pure NR data
else
    nDIIS = min(Iter_GEK,nWindow) !1 for first iteration; 2
end if

! index of the first iteration to consider for the subspace
iFirst = Iter_GEK-nDIIS+1 !1 for first iteration; 1

call mma_Allocate(q,fsdim, nDiis,Label="q")
call mma_Allocate(g,fsdim, nDiis,Label="g")

! compute product matrix U_1...n = U_1 * ... * U_n
! -------------------------------------------------
! look later if this can be done in pipekmezey_iter, to save comp
call mma_allocate(UmatProd,nOrb2Loc,nOrb2Loc,Label='UmatProd')
call mma_allocate(xUmatProd,nOrb2Loc,nOrb2Loc,Label='xUmatProd')
call mma_allocate(Umat_i,nOrb2Loc,nOrb2Loc,Label='Umat_i')

xUmatProd(:,:) = Zero
Umat_i(:,:) = Zero
call unitmat(xUmatProd,nOrb2Loc)

do i=iFirst,Iter_GEK
    Umat_i(:,:) = UmatList(:,:,i)
    !call RecPrt("UmatList(:,:,i) = ",' ',UmatList(:,:,i),nOrb2Loc,nOrb2Loc)

    !call RecPrt("U_i = "," ",Umat_i,nOrb2Loc,nOrb2Loc)
    call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,One,xUmatProd,nOrb2Loc,Umat_i,nOrb2Loc,Zero,UmatProd,norb2Loc)

    !call RecPrt("U_1...i = "," ",UmatProd,nOrb2Loc,nOrb2Loc)

    xUmatProd(:,:) = UmatProd(:,:)
end do
!call RecPrt("U_1...n = "," ",UmatProd,nOrb2Loc,nOrb2Loc)
write(u6,*) "ndiis=",ndiis

! -------------------------------------------------
call mma_allocate(disp_summed,fsdim,Label="disp_summed")
call mma_allocate(UmatKsum,nOrb2Loc,nOrb2Loc,Label="UmatKsum")
disp_summed(:) = Zero
UmatKsum(:,:) = Zero

j = 0
do i=iFirst,Iter_GEK
    j = i-iFirst+1
    !write(u6,*) 'i,j,iter=',i,j,Iter_GEK

    ! Coordinates
    q(:,j) = DispList(:,i)
    disp_summed(:) = disp_summed(:) + DispList(:,i)

    ! Gradients
    g(:,j) = GradList(:,i)

end do
write(u6,*) "n =",j


call mma_allocate(kappa_summed,nOrb2Loc,nOrb2Loc,Label="kappa_summed")


call vec2upper_triag(kappa_summed,nOrb2Loc,disp_summed,fsdim,.true.)
call expkap_localisation(kappa_summed,nOrb2Loc,Umat_i,xUmatProd,UmatKsum)


!call RecPrt("exp(-K_1-K_2-...-K_n) = "," ",UmatKsum,nOrb2Loc,nOrb2Loc)
!call RecPrt("U_1...n - exp(-K_1-K_2-...-K_n) = "," ",UmatProd-UmatKsum,nOrb2Loc,nOrb2Loc)

norm = 0
norm = sqrt(DDot_(nOrb2Loc*nOrb2Loc,UmatProd-UmatKsum,1,UmatProd-UmatKsum,1))
write(u6,*) "norm of U_1...n - exp(-K_1-K_2-...-K_n):", norm
norm = 0

call mma_Deallocate(xUmatProd)
call mma_Deallocate(Umat_i)
call mma_Deallocate(UmatProd)
call mma_Deallocate(UmatKsum)
call mma_Deallocate(disp_summed)
call mma_Deallocate(kappa_summed)

if (nDIIS < 3) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit S-GEK Optimizer'
# endif
  call mma_deallocate(g)
  call mma_deallocate(q)
  return
end if

#ifdef _DEBUGPRINT_
    write(u6,*) 'nWindow =',nWindow
    write(u6,*) '  nDIIS =',nDIIS
    write(u6,*) '  Iter_GEK =',Iter_GEK
    call RecPrt("q(:,:)",' ',q,fsdim, nDiis)
    call RecPrt("g(:,:)",' ',g,fsdim, nDiis)
    call RecPrt("g(:,nDiis)",' ',g(:,nDiis),fsdim, 1)
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

!number of subspace basis vectors, potentially linear dependent => difference vecs of ndiis displacements and gradients +2 additional vecs (see below)
nExplicit = 2*(nDIIS-1)+2

call mma_allocate(e_diis,fsdim,nExplicit,Label='e_diis')
call mma_allocate(Aux_a,fsdim,Label='Aux_a')
call mma_allocate(Aux_b,fsdim,Label='Aux_b')

j = 0
thr = 1E-18_wp
do k=1,nDIIS-1
    !n-th column of e_diis
    j = j+1
    ! gradient difference vector
    Aux_a(:) = g(:,k+1)-g(:,k)
    !normalize
    norm = sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))
    if (norm < thr) then
        e_diis(:,j) = Zero
    else
        e_diis(:,j) = Aux_a(:)/norm
    end if

    !(n+1)-th column of e_diis
    j = j+1
    ! displacement difference vector
    Aux_a(:) = q(:,k+1)-q(:,k)
    Aux_b(:) = Aux_a(:)
    !normalize
    norm = sqrt(DDot_(fsdim,Aux_b(:),1,Aux_b(:),1))
    if (norm < thr) then
        e_diis(:,j) = Zero
    else
        e_diis(:,j) = Aux_b(:)/norm
    end if


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


#ifdef _DEBUGPRINT_
    if (allocated(e_diis)) call RecPrt('e_diis(unorth)',' ',e_diis,fsdim,nExplicit)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
write(u6,*) '    Iter_GEK:',Iter_GEK
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

#ifdef _DEBUGPRINT_
    call RecPrt('q_diis',' ',q_diis,mDIIS,nDIIS)
    call RecPrt('g_diis',' ',g_diis,mDIIS,nDIIS)
    !call RecPrt('H_diis(HDiag)',' ',H_diis,mDIIS,mDIIS)
    write(u6,'(/A)') "H_diis diagonal elements:"
    do i=1,mDiis
        write(u6,'(F26.16)') H_diis(i,i)
    end do
    write(u6,*) ''
#endif


if (usmitigation) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Undershoot avoidance: Scale along dq !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (Loosen%Factor /= One) then
  ! Components of dq in the subspace
  call mma_allocate(w,mDIIS,mDIIS,Label='w')
  gg = sqrt(DDot_(fsdim,dq,1,dq,1))
  do i=1,mDIIS
    w(i,1) = DDot_(fsdim,dq,1,e_diis(:,i),1)/gg
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

# ifdef _DEBUGPRINT_
  write(u6,*) "Undershoot mitigation, scale along dq"
  call RecPrt('H_diis(after scaling)',' ',H_diis,mDIIS,mDIIS)
# endif

end if
end if



! build the surrogate model & perform the optimization
! ----------------------------------------------------


if (SORange) then
  SOFact = One
else
  SOFact = 10000000.0_wp
end if
Call GEK_Optimizer(mDiis,nDiis,Max_Iter_GEK,q_diis(:,:),g_diis(:,:),dq_diis(:),FuncList(iFirst:),H_diis(:,:),dqdq,&
                   Step_Trunc,UpMeth,SOFact,10.0_wp)

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



!compute angle
norm = sqrt(DDot_(fsdim,dq_NR,1,dq_NR,1))
write(u6,'(A,F12.6,2X,A,F12.3,2x,A,I4)') "Angle(dq_NR,dq) (deg):", acos(DDot_(fsdim,dq_NR,1,dq,1)/(norm*dqdq))/Pi*180.0_wp,&
                                         "norm(dq)/norm(dq_NR)",dqdq/norm, "Iter_GEK=",Iter_GEK

#ifdef _DEBUGPRINT_
    write(u6,*) '||dq||=',dqdq
    call RecPrt('dq(:) after projecting out',' ',dq(:),size(dq),1)
#endif

! deallocations
! -------------
call mma_Deallocate(q)
call mma_Deallocate(g)
call mma_Deallocate(dq_NR)

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

end subroutine S_GEK_localisation
