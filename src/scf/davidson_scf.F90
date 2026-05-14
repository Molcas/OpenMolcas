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
! Copyright (C) 2014, Ignacio Fdez. Galvan                             *
!               2017, Roland Lindh                                     *
!***********************************************************************
!  Davidson_SCF
!
!> @brief
!>   Compute the lowest \p k eigenvalues of a symmetric matrix.
!> @author Ignacio Fdez. Galv&aacute;n
!> @modified_by Roland Lindh
!>
!> @details
!> Simple application of the Davidson procedure to obtain the lowest \p k eigenvalues
!> and corresponding eigenvectors of a symmetric matrix.
!> On input, \p Vec can contain an initial guess for the eigenvectors (from a previous
!> run with smaller \p k, for example), only the non-zero vectors are used.
!> This routine is adapted to an augmented Hessian, which is not explicitly expressed,
!> rather the original Hessian is implicitly there, via a diagonal and an on-the-fly
!> update when multiplied by a vector, and the gradient is explicit there.
!>
!> @param[in]     HDiag  Diagonal of the Hessian matrix
!> @param[in]     g      Gradient vector
!> @param[in]     m      Size of diagonal Hessian and gradient
!> @param[in]     k      Number of lowest eigenvalues to compute
!> @param[in]     Fact   Scaling factor
!> @param[out]    Eig    Lowest eigenvalues
!> @param[in,out] Vec    Lowest eigenvectors
!> @param[out]    iRC    Return code (0 if converged)
!***********************************************************************

!#define _DEBUGPRINT_
!#define _DEBUGCode_
subroutine Davidson_SCF(g,m,k,Fact,Eig,Vec,iRC)

#ifdef _DEBUGCode_
use Index_Functions, only: iTri, nTri_Elem
#endif
use InfSCF, only: HDiag
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: m, k
real(kind=wp), intent(in) :: g(m), Fact
real(kind=wp), intent(out) :: Eig(k)
real(kind=wp), intent(inout) :: Vec(m+1,k)
integer(kind=iwp), intent(out) :: iRC
integer(kind=iwp) :: i, ig, ii, info, iter, j, jj, maxk, mink, mk, n, nTmp, old_mk
real(kind=wp) :: Alpha, Aux, Conv, Dum(1) = Zero, tmp
logical(kind=iwp) :: Augmented, Last, Reduced
integer(kind=iwp), allocatable :: Index_D(:)
real(kind=wp), allocatable :: Ab(:,:), Diag(:), Eig_old(:), EVal(:), EVec(:), Proj(:), Sub(:,:), TAV(:), TmpVec(:), TRes(:), TVec(:)
integer(kind=iwp), parameter :: maxiter = 100
real(kind=wp), parameter :: Thr = 1.0e-6_wp, Thr2 = 1.0e-12_wp, Thr3 = 1.0e-16_wp
real(kind=wp), external :: ddot_
#ifdef _DEBUGCode_
integer(kind=iwp) :: ij
real(kind=wp), allocatable :: EVal(:), EVec(:), HAug(:,:), HM(:,:), Vec_(:)
#endif
n = m+1

#ifdef _DEBUGCode_
call mma_allocate(Vec_,m,Label='Vec_')
call mma_allocate(HM,m,m,Label='HM')
HM(:,:) = Zero
call mma_allocate(HAug,n,n,Label='HAug')

do i=1,m
  Vec_(:) = Zero
  Vec_(i) = One
  call SOrUpV(Vec_(:),m,HM(:,i),'GRAD','BFGS')
end do
!call RecPrt('HM',' ',HM,m,m)

call mma_allocate(EVal,nTri_Elem(m),Label='EVal')
call mma_allocate(EVec,m**2,Label='EVec')
do i=1,m
  do j=1,i
    EVal(iTri(i,j)) = HM(i,j)
  end do
end do

! Set up a unit matrix

call unitmat(EVec,m)

! Compute eigenvalues and eigenvectors

call NIDiag_new(EVal,EVec,m,m)
call Jacord(EVal,EVec,m,m)

!do i=1,m
!  ij = iTri(i,j)
!  write(u6,*) 'Eval0(ij)=',EVal(ij)
!end do

call mma_deallocate(EVal)
call mma_deallocate(EVec)

HAug(n,1:m) = g(:)
HAug(1:m,n) = g(:)
HAug(1:m,1:m) = HM(:,:)

call mma_allocate(EVal,nTri_Elem(n),Label='EVal')
call mma_allocate(EVec,n**2,Label='EVec')

do i=1,n
  do j=1,i
    EVal(iTri(i,j)) = HAug(i,j)
  end do
end do

! Set up a unit matrix

call unitmat(EVec,n)

! Compute eigenvalues and eigenvectors

call NIDiag_new(EVal,EVec,n,n)
call Jacord(EVal,EVec,n,n)

do i=1,n
  ij = iTri(i,j)
  if (EVal(ij) < Zero) write(u6,*) 'Eval(ij)=',EVal(ij)
end do

call mma_deallocate(EVal)
call mma_deallocate(EVec)
call mma_deallocate(HAug)
call mma_deallocate(HM)
call mma_deallocate(Vec_)
#endif

#ifdef _DEBUGPRINT_
call NrmClc(HDiag,m,'Davidson_SCF','HDiag')
call NrmClc(g,m,'Davidson_SCF','g')
!call RecPrt('HDiag',' ',HDiag,1,m)
!call RecPrt('g',' ',g,1,m)
#endif

! Initialize some parameters
!  mk   = subspace size (initially k)
!  maxk = maximum subspace size (25 if k=1)
!  mink = subspace size to reduce to when the maximum is exceeded (5 if k=1)

if (k > n) call SysAbendMsg('Davidson_SCF','Wrong k value.','')

mink = min(max(k+2,5),n)
maxk = min(5*mink,n)
mk = k
iRC = 0

! Allocate matrices
!  Sub  = Vectors (columns) defining the subspace (maximum maxk vectors)
!  Ab   = A*b vectors (A * Sub)
!  EVal = Set of computed eigenvalues (maximum maxk elements)
!  EVec = Set of computed eigenvectors, in the subspace (maximum maxk*maxk)

call mma_allocate(Eig_old,k,label='Eig_old')
call mma_allocate(Sub,n,maxk,Label='Sub')
call mma_allocate(Ab,n,maxk,Label='Ab ')
call mma_allocate(Proj,maxk*maxk,Label='Proj')
call mma_allocate(EVal,maxk,Label='EVal')
call mma_allocate(EVec,maxk*maxk,Label='EVec')
Ab(:,:) = Zero
EVal(:) = Zero
EVec(:) = Zero

! Build an index of sorted diagonal elements in A

call mma_allocate(Index_D,n,Label='Index_D')
do i=1,n
  Index_D(i) = i
end do

do i=1,n
  ii = Index_D(i)
  if (ii == n) then
    Aux = Zero
  else
    Aux = HDiag(ii)
  end if
  ii = i
  do j=i,n
    jj = Index_D(j)
    if (jj == n) then
      Tmp = Zero
    else
      Tmp = HDiag(jj)
    end if
    if (Tmp < Aux) then
      Aux = Tmp
      ii = j
    end if
  end do
  if (ii /= i) then
    jj = Index_D(ii)
    Index_D(ii) = Index_D(i)
    Index_D(i) = jj
  end if
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'Index_D=',Index_D
#endif

! Setup the initial subspace
!  Read the non-linear-dependent columns from the initial eigenvector matrix
!  Fill up to mk with selected base vectors from the initial matrix
!   (those corresponding to the lowest diagonal elements)
!  The rest is set to zero, just in case

nTmp = 0
call mma_allocate(TmpVec,n,Label='TmpVec')
do i=1,k
  TmpVec(:) = Vec(1:n,i)
  call Add_Vector(n,nTmp,Sub,TmpVec,Thr3)
end do

ii = 0
TmpVec(:) = Zero

do while ((nTmp < mk) .and. (ii < n))
  ii = ii+1
  jj = Index_D(ii)

  ! A large value indicates a forbidden rotation.
  ! A negative value indicates large rotation to another
  ! global minimum. Avoid these!!!

  if (jj == n) then
    Aux = Zero
  else
    Aux = HDiag(jj)
  end if
  if ((Aux < 1.0e10_wp) .and. (Aux > -0.1_wp)) then
    TmpVec(jj) = One
    call Add_Vector(n,nTmp,Sub,TmpVec,Thr3)
    TmpVec(jj) = Zero
  end if
end do

! ig will be a global counter to loop across all n base vectors
ig = ii
call mma_deallocate(TmpVec)
Sub(:,mk+1:) = Zero

! Iterative procedure starts here
!  mk     = subspace size at each iteration
!  old_mk = value of mk at the previous iteration

Augmented = .false.
Reduced = .false.
Last = .false.
old_mk = 0
iter = 0
call mma_allocate(Diag,n,Label='Diag')
call mma_allocate(TVec,n,Label='TVec')
call mma_allocate(TAV,n,Label='TAV ')
call mma_allocate(TRes,n,Label='TRes')
do while (.not. Last)
  iter = iter+1
  if (iter > 1) Eig_old(:) = Eig(:)
# ifdef _DEBUGPRINT_
  if (.not. Reduced) then
    write(u6,'(A)') '---------------'
    write(u6,'(A,1X,I5)') 'Iteration',iter
  end if
  !call RecPrt('Orthonormalized subspace',' ',Sub,n,mk)
# endif

  ! Compute the matrix product
  !  Ab = A * Sub
  !  Only the new vectors since the last iterations need to be calculated

  ! Note that the A-matrix is the augmented Hessian of a rs-rfo
  ! approach. The A matrix is not explicitly stored but rather only
  ! the associated gradient is. The original Hessian is implicitly
  ! there and a vector corresponding to the contraction of the updated
  ! Hessian and a trial vector can be computed on-the-fly.

  do j=old_mk,mk-1
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Davidson_SCF: j,Fact=',j,Fact
    call NrmClc(Sub(1,j+1),n,'Davidson_SCF','Sub(1,j+1)')
    !call RecPrt('Sub',' ',Sub(1,j+1),1,n)
#   endif

    ! Pick up the contribution for the updated Hessian (BFGS update)
    ! Add contribution from the gradient

    call SOrUpV(Sub(1,j+1),m,Ab(1,j+1),'GRAD','BFGS')
    !call RecPrt('Ab(0)',' ',Ab(1,j+1),1,n)
    tmp = Sub(n,j+1)/sqrt(Fact)
    Ab(1:m,j+1) = Ab(1:m,j+1)/Fact+tmp*g(:)

    Ab(n,j+1) = DDot_(m,g,1,Sub(1,j+1),1)/sqrt(Fact)
#   ifdef _DEBUGPRINT_
    call NrmClc(Ab(1,j+1),n,'Davidson_SCF','Ab(1,j+1)')
    !call RecPrt('Ab',' ',Ab(1,j+1),1,n)
#   endif

  end do

  ! Compute the matrix to diagonalize (symmetric)
  !  Proj = Sub^t * Ab
  !  Again, only the new rows/columns are needed

  if (old_mk == 0) then
    call DGeMM_('T','N',mk,mk,n,One,Sub,n,Ab,n,Zero,Proj,maxk)
  else
    do i=0,mk-1
      do j=max(old_mk,i),mk-1
        Proj(1+i*maxk+j) = DDot_(n,Sub(:,j+1),1,Ab(:,i+1),1)
      end do
    end do
  end if

  ! Compute the eigenvalues of the projected matrix
  !  Make sure the eigenpairs are sorted
  !  If the subspace has been reduced, no need to compute new eigenpairs

  if (.not. Reduced) then
#   ifdef _DEBUGPRINT_
    write(u6,'(2X,A,1X,I5)') 'Solving for subspace size:',mk
#   endif
    EVec(:) = Proj(:)
    call dsyev_('V','L',mk,EVec,maxk,EVal,Dum,-1,info)
    if (info /= 0) then
      write(u6,*) 'info(2)/=0',info
      call Abend()
    end if
    nTmp = max(1,int(Dum(1)))
    call mma_allocate(TmpVec,nTmp,Label='TmpVec')
    call dsyev_('V','L',mk,EVec,maxk,EVal,TmpVec,nTmp,info)
    if (info /= 0) then
      write(u6,*) 'info(2)/=0',info
      call Abend()
    end if
    call mma_deallocate(TmpVec)
    call SortEig(EVal,EVec,mk,maxk,1,.false.)
    Eig(:) = EVal(1:k)
#   ifdef _DEBUGPRINT_
    !call RecPrt('Current guess',' ',Eig,1,k)
#   endif
  end if
# ifdef _DEBUGPRINT_
  !call RecPrt('Eigenvalues',' ',EVal,1,mk)
  !call SubRecPrt('Subspace Eigenvectors',' ',EVec,maxk,mk,mk)
  write(u6,*)
# endif
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  ! Check for convergence                                              *
  !
  ! Converge if the change in the eigenvalues is small
  ! (but if a mink size has been reached)
  ! Converge if the full system has been solved
  ! Stop if the number of iterations exceeds the maximum
  ! Stop if no new vectors to add are found
  !
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute Conv

  if (iter > 1) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    Conv = Zero
    do i=1,k
      if (abs(Eig(i)) > Thr2) then
        Conv = max(Conv,abs((Eig(i)-Eig_old(i))/Eig(i)))
      else
        Conv = max(Conv,abs(Eig(i)-Eig_old(i)))
      end if
    end do
#   ifdef _DEBUGPRINT_
    if (Augmented) write(u6,'(2X,A,1X,G12.6)') 'Maximum relative eigenvalue change:',Conv
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    Conv = Ten*Thr
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now check for convergence

  old_mk = mk
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Augmented .and. (Conv <= Thr) .and. (mk >= mink)) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Converged due to small change'
#   endif
    Last = .true.
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (mk == n) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Complete system solved'
#   endif
    Last = .true.
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (iter >= maxiter) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Not converged'
#   endif
    Last = .true.
    iRC = 1

    ! Reduce the subspace size if it exceeds the maximum (maxk)
    !  Sub' = Sub * Vec(1:mink)
    !  Sub' should be orthonormal if Sub is orthonormal
    ! (A reduction does not consume an iteration)
    ! There is also a reduction if the process is stagnated

    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (mk > mink .and. (min(mk+k,n) > maxk) .or. (iRC == 2)) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (iRC == 2) iRC = 0
#   ifdef _DEBUGPRINT_
    write(u6,'(2X,A,1X,I5)') 'Reducing search space to',mink
#   endif
    call mma_allocate(TmpVec,mink*n,Label='TmpVec')
    call DGeMM_('N','N',n,mink,mk,One,Sub,n,EVec,maxk,Zero,TmpVec,n)
    Sub(:,1:mink) = reshape(TmpVec,[n,mink])
    call mma_deallocate(TmpVec)

    ! To make sure Sub' is orthonormal, add the vectors one by one

    j = 0
    i = 0
    do while ((j < mink) .and. (i < mk))
      i = i+1
      call Add_Vector(n,j,Sub,Sub(1,i),Thr3)
    end do

    ! j should be mink, but who knows...

#   ifdef _DEBUGPRINT_
    if (j < mink) write(u6,'(2X,A,1X,I5)') 'Fewer vectors found:',j
#   endif
    Sub(:,j+1:) = Zero
    Ab(:,j+1:) = Zero
    EVec(:) = Zero
    do i=0,j-1
      EVec(1+i*(maxk+1)) = One
    end do
    mk = j
    old_mk = 0
    Augmented = .false.
    Reduced = .true.
    iter = iter-1
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Expand the subspace
    !  For each eigenpair i of the first k,
    !  check convergence for the residuals r:
    !   r = Ab * Vec(i) - Val(i) * Sub * Vec(i)
    !  Add a new vector, orthonormalized with the previous vectors,
    !  computed from r and the eigenpair
    !  (different possible variants)

    call mma_allocate(TmpVec,n,Label='TmpVec')
    Conv = Zero

    jj = 0
    do i=0,k-1
      ! Vector in full space: Sub*Vec(i)
      call dGeMV_('N',n,mk,One,Sub,n,EVec(1+i*maxk),1,Zero,TVec,1)
      ! Product of matrix and vector: Ab*Vec(i)
      call dGeMV_('N',n,mk,One,Ab,n,EVec(1+i*maxk),1,Zero,TAV,1)
      ! Residual: (A-Val(i))*Vec(i) = Ab*Vec(i) - Val(i)*Sub*Vec(i)
      TRes(:) = TAV(:)-EVal(1+i)*TVec(:)
      Conv = max(Conv,DDot_(n,TRes,1,TRes,1))

      ! Scale vector, orthonormalize, and add to subspace

      ! Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
      do j=0,n-1
        !Aux = A(nTri_Elem(j+1))-EVal(1+i)
        if (j == n-1) then
          Aux = -Eval(1+i)
        else
          Aux = HDiag(j+1)-Eval(1+i)
        end if
        if (j == n-1) then
          Diag(1+j) = One/sign(max(abs(Aux),Thr2),Aux)
        else
          if (HDiag(j+1) < 1.0e20_wp) then
            Diag(1+j) = One/sign(max(abs(Aux),Thr2),Aux)
          else
            Diag(1+j) = 1.0e20_wp
          end if
        end if
      end do

      ! scale
      do j=0,n-1
        if (Diag(1+j) < Ten**2) then
          TmpVec(1+j) = TRes(1+j)*Diag(1+j)
        else
          TmpVec(1+j) = Zero
        end if
      end do

      Alpha = Zero
      do j=0,n-1
        if (Diag(1+j) < Ten**2) Alpha = Alpha+Diag(1+j)*TVec(1+j)**2
      end do
      Alpha = DDot_(n,TVec,1,TmpVec,1)/Alpha
      ! subtract
      do j=0,n-1
        if (Diag(1+j) < Ten**2) then
          TVec(1+j) = TVec(1+j)*Diag(1+j)
        else
          TVec(1+j) = Zero
        end if
      end do
      TmpVec(:) = TmpVec(:)-Alpha*TVec(:)

      if (mk+jj <= n-1) then
        jj = mk+jj
        call Add_Vector(n,jj,Sub,TmpVec,Thr3)
        jj = jj-mk
      end if
    end do

#   ifdef _DEBUGPRINT_
    write(u6,'(2X,A,1X,G12.6)') 'Maximum residual:',Conv
#   endif
    !                                                                  *
    !------------------------------------------------------------------*
    !                                                                  *
    if ((Conv < Thr3) .and. (mk >= mink)) then
      !                                                                *
      !----------------------------------------------------------------*
      !                                                                *
#     ifdef _DEBUGPRINT_
      write(u6,'(A)') 'Converged due to small residual'
#     endif
      Last = .true.
      !                                                                *
      !----------------------------------------------------------------*
      !                                                                *
    else
      !                                                                *
      !----------------------------------------------------------------*
      !                                                                *
      mk = min(mk+jj,n)

      ! If no new vector is found to add to the subspace, we are in
      ! trouble. Try to find a non-linear-dependent base vector in
      ! the original matrix

      if (jj == 0) then
#       ifdef _DEBUGPRINT_
        write(u6,'(A)') 'Process stagnated'
#       endif
        if (mk < maxk) then
          TmpVec(:n) = Zero
          i = 0

          do while ((jj < 1) .and. (i < n))
            i = i+1
            ig = mod(ig,n)+1
            ii = Index_D(ig)

            ! Avoid explicitly rotations between fermions of
            ! different types. Avoid rotations which will be
            ! large,

            if (ii == n) then
              Aux = Zero
            else
              Aux = HDiag(ii)
            end if
            if ((Aux < 1.0e20_wp) .and. (Aux > -0.1_wp)) then
              TmpVec(ii) = One
              jj = mk+jj
              call Add_Vector(n,jj,Sub,TmpVec,Thr3)
              TmpVec(ii) = Zero
              jj = jj-mk
            end if

          end do

          mk = min(mk+jj,n)
          if (jj > 0) Augmented = .true.
        end if
        if (jj == 0) then
          Augmented = .false.
          iRC = 2
        end if
      else
        Augmented = .true.
      end if
      !                                                                *
      !----------------------------------------------------------------*
      !                                                                *
    end if
    !                                                                  *
    !------------------------------------------------------------------*
    !                                                                  *
    call mma_deallocate(TmpVec)
    Reduced = .false.
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do

call mma_deallocate(Diag)
call mma_deallocate(TVec)
call mma_deallocate(TAV)
call mma_deallocate(TRes)
call mma_deallocate(Index_D)

! Store the current lowest k eigenvectors (in the full space)
! Vec' = Sub * Vec(1:k)

call DGeMM_('N','N',n,k,mk,One,Sub,n,EVec,maxk,Zero,Vec,n)

#ifdef _DEBUGPRINT_
call NrmClc(Vec(1,1),m,'Davidson_SCF','Vec(1-m)')
call NrmClc(Vec(n,1),1,'Davidson_SCF','Vec(n)')
#endif
call mma_deallocate(Sub)
call mma_deallocate(Ab)
call mma_deallocate(Proj)
call mma_deallocate(EVal)
call mma_deallocate(EVec)
call mma_deallocate(Eig_old)

end subroutine Davidson_SCF
