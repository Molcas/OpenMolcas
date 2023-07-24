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
!***********************************************************************
!  Davidson
!
!> @brief
!>   Compute the lowest \p k eigenvalues of a symmetric matrix.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Simple application of the Davidson procedure to obtain the lowest \p k eigenvalues
!> and corresponding eigenvectors of a symmetric matrix.
!> On input, \p Vec can contain an initial guess for the eigenvectors (from a previous
!> run with smaller \p k, for example), only the non-zero vectors are used.
!>
!> @param[in]     A   Symmetric matrix, in upper triangular packed format
!> @param[in]     n   Size of the matrix
!> @param[in]     k   Number of lowest eigenvalues to compute
!> @param[out]    Eig Lowest eigenvalues
!> @param[in,out] Vec Lowest eigenvectors
!> @param[out]    iRC Return code (0 if converged)
!***********************************************************************

subroutine Davidson(A,n,k,Eig,Vec,iRC)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k
real(kind=wp), intent(in) :: A(nTri_Elem(n))
real(kind=wp), intent(out) :: Eig(k)
real(kind=wp), intent(inout) :: Vec(n,k)
integer(kind=iwp), intent(out) :: iRC
integer(kind=iwp) :: i, ig, ii, info, iter, j, jj, maxk, mink, mk, nTmp, old_mk
real(kind=wp) :: Alpha, Aux, Conv, rDum(1)
logical(kind=iwp) :: Augmented, Last, Reduced
integer(kind=iwp), allocatable :: Indx(:)
real(kind=wp), allocatable :: Ab(:), Diag(:), Eig_old(:), EVal(:), EVec(:), Proj(:), Sub(:), TAV(:), Tmp(:), TRes(:), TVec(:), &
                              Val(:), Vec2(:)
integer(kind=iwp), parameter :: maxiter = 300
real(kind=wp), parameter :: Thr = 1.0e-7_wp, Thr2 = 1.0e-16_wp, Thr3 = 1.0e-16_wp
real(kind=wp), external :: ddot_

!#define _DEBUGPRINT_

! Diagonal preconditioned residue (Davidson)
#define DAV_DPR 1
! Inverse-iteration generalized Davidson (Olsen et al.)
#define DAV_IIGD 2
! Generalized Jacobi-Davidson (Sleijpen et al.) (DO NOT USE IT)
#define DAV_GJD 3

#ifdef _DEBUGPRINT_
call TriPrt('Initial matrix','',A,n)
#endif

! Initialize some parameters
! mk   = subspace size (initially k)
! maxk = maximum subspace size (25 if k=1)
! mink = subspace size to reduce to when the maximum is exceeded (5 if k=1)

if (k > n) call SysAbendMsg('Davidson','Wrong k value.','')
mink = min(max(k+2,5),n)
maxk = min(5*mink,n)
mk = k
iRC = 0

! If all the eigenvalues are wanted, better solve the system directly
! and return

if (mk >= n) then
  call mma_allocate(Val,n,Label='Val')
  call mma_allocate(Vec2,n*n,Label='Vec2')
  do j=1,n
    jj = (j-1)*n
    Vec2(jj+1:jj+j) = A(nTri_Elem(j-1)+1:nTri_Elem(j))
  end do
  call dsyev_('V','U',n,Vec2,n,Val,rDum,-1,info)
  nTmp = int(rDum(1))

  call mma_allocate(Tmp,nTmp,Label='Tmp')
  call dsyev_('V','U',n,Vec2,n,Val,Tmp,nTmp,info)
  call SortEig(Val,Vec2,n,n,1,.false.)
  call mma_deallocate(Tmp)

  Eig(:) = Val(1:k)
  Vec(:,:) = reshape(Vec2(1:n*k),[n,k])

  call mma_deallocate(Val)
  call mma_deallocate(Vec2)
# ifdef _DEBUGPRINT_
  call RecPrt('Eigenvalues',' ',Eig,1,n)
  write(u6,*)
  write(u6,'(A)') 'Complete system solved'
# endif
  return
end if

! Allocate matrices
! Sub  = Vectors (columns) defining the subspace (maximum maxk vectors)
! Ab   = A*b vectors (A * Sub)
! EVal  = Set of computed eigenvalues (maximum maxk elements)
! EVec  = Set of computed eigenvectors, in the subspace (maximum maxk*maxk)

call mma_allocate(Eig_old,k,label='Eig_old')
call mma_allocate(Sub,n*maxk,Label='Sub')
call mma_allocate(Ab,n*maxk,Label='Ab ')
call mma_allocate(Proj,maxk*maxk,Label='Proj')
call mma_allocate(EVal,maxk,Label='EVal')
call mma_allocate(EVec,maxk*maxk,Label='EVec')
AB(:) = Zero
EVal(:) = Zero
EVec(:) = Zero

! Build an index of sorted diagonal elements in A

call mma_allocate(Indx,n,Label='Indx')
do i=1,n
  Indx(i) = i
end do
do i=1,n
  ii = Indx(i)
  Aux = A(nTri_Elem(ii))
  ii = i
  do j=i,n
    jj = Indx(j)
    if (A(nTri_Elem(jj)) < Aux) then
      Aux = A(nTri_Elem(jj))
      ii = j
    end if
  end do
  if (ii /= i) then
    jj = Indx(ii)
    Indx(ii) = Indx(i)
    Indx(i) = jj
  end if
end do

! Setup the initial subspace
! Read the non-linear-dependent columns from the initial eigenvector matrix
! Fill up to mk with selected base vectors from the initial matrix
!  (those corresponding to the lowest diagonal elements)
! The rest is set to zero, just in case

nTmp = 0
call mma_allocate(Tmp,n,Label='Tmp')
do i=1,k
  Tmp(:) = Vec(:,i)
  call Add_Vector(n,nTmp,Sub,Tmp,Thr3)
end do

ii = 0
Tmp(:) = Zero
do while ((nTmp < mk) .and. (ii < n))
  ii = ii+1
  jj = Indx(ii)
  Tmp(jj) = One
  call Add_Vector(n,nTmp,Sub,Tmp,Thr3)
  Tmp(jj) = Zero
end do
! ig will be a global counter to loop across all n base vectors
ig = ii
call mma_deallocate(Tmp)
Sub(n*mk+1:) = Zero

! Iterative procedure starts here
! mk     = subspace size at each iteration
! old_mk = value of mk at the previous iteration

Augmented = .false.
Reduced = .false.
Last = .false.
old_mk = 0
iter = 0
call mma_allocate(Diag,n,Label='Diag')
call mma_allocate(TVec,n,Label='TVec')
call mma_allocate(TAV,n,Label='TAV')
call mma_allocate(TRes,n,Label='TRes')
do while (.not. Last)
  iter = iter+1
  if (iter > 1) Eig_old(:) = Eig(:)
# ifdef _DEBUGPRINT_
  if (.not. Reduced) then
    write(u6,'(A)') '---------------'
    write(u6,'(A,1X,I5)') 'Iteration',iter
  end if
  call RecPrt('Orthonormalized subspace',' ',Sub,n,mk)
# endif

  ! Compute the matrix product
  ! Ab = A * Sub
  ! Only the new vectors since the last iterations need to be calculated

  call mma_allocate(Tmp,n,Label='Tmp')
  do i=1,n
    ! Reconstruct a row (or column) of the matrix A
    Tmp(1:i) = A(nTri_Elem(i-1)+1:nTri_Elem(i))
    do j=i+1,n
      Tmp(j) = A(nTri_Elem(j-1)+i)
    end do
    ! Compute the i-th element of each new column
    do j=old_mk,mk-1
      Ab(j*n+i) = DDot_(n,Tmp,1,Sub(1+j*n),1)
    end do
  end do
  call mma_deallocate(Tmp)

  ! Compute the matrix to diagonalize (symmetric)
  ! Proj = Sub^t * Ab
  ! Again, only the new rows/columns are needed

  if (old_mk == 0) then
    call DGeMM_('T','N',mk,mk,n,One,Sub,n,Ab,n,Zero,Proj,maxk)
  else
    do i=0,mk-1
      do j=max(old_mk,i),mk-1
        Proj(1+i*maxk+j) = DDot_(n,Sub(1+j*n),1,Ab(1+i*n),1)
      end do
    end do
  end if

  ! Compute the eigenvalues of the projected matrix
  ! Make sure the eigenpairs are sorted
  ! If the subspace has been reduced, no need to compute new eigenpairs

  if (.not. Reduced) then
#   ifdef _DEBUGPRINT_
    write(u6,'(2X,A,1X,I5)') 'Solving for subspace size:',mk
#   endif
    EVec(:) = Proj(:)
    call dsyev_('V','L',mk,EVec,maxk,EVal,rDum,-1,info)
    nTmp = int(rDum(1))
    call mma_allocate(Tmp,nTmp,Label='Tmp')
    call dsyev_('V','L',mk,EVec,maxk,EVal,Tmp,nTmp,info)
    call mma_deallocate(Tmp)
    call SortEig(EVal,EVec,mk,maxk,1,.false.)
    Eig(:) = EVal(1:k)
#   ifdef _DEBUGPRINT_
    call RecPrt('Current guess',' ',Eig,1,k)
#   endif
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('Eigenvalues',' ',EVal,1,mk)
  call SubRecPrt('Subspace Eigenvectors',' ',EVec,maxk,mk,mk)
  write(u6,*)
# endif

  ! Check for convergence
  ! Converge if the change in the eigenvalues is small
  !  (but if a mink size has been reached)
  ! Converge if the full system has been solved
  ! Stop if the number of iterations exceeds the maximum
  ! Stop if no new vectors to add are found

  if (iter > 1) then
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
  else
    Conv = Ten*Thr
  end if
  old_mk = mk
  if (Augmented .and. (Conv <= Thr) .and. (mk >= mink)) then
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Converged due to small change'
#   endif
    Last = .true.
  else if (mk == n) then
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Complete system solved'
#   endif
    Last = .true.
  else if (iter >= maxiter) then
#   ifdef _DEBUGPRINT_
    write(u6,'(A)') 'Not converged'
#   endif
    Last = .true.
    iRC = 1

    ! Reduce the subspace size if it exceeds the maximum (maxk)
    ! Sub' = Sub * Vec(1:mink)
    ! Sub' should be orthonormal if Sub is orthonormal
    ! (A reduction does not consume an iteration)
    ! There is also a reduction if the process is stagnated

  else if ((min(mk+k,n) > maxk) .or. (iRC == 2)) then
    if (iRC == 2) iRC = 0
#   ifdef _DEBUGPRINT_
    write(u6,'(2X,A,1X,I5)') 'Reducing search space to',mink
#   endif
    call mma_allocate(Tmp,mink*n,Label='Tmp')
    call DGeMM_('N','N',n,mink,mk,One,Sub,n,EVec,maxk,Zero,Tmp,n)
    Sub(1:mink*n) = Tmp(:)
    call mma_deallocate(Tmp)

    ! To make sure Sub' is orthonormal, add the vectors one by one

    j = 0
    i = 0
    do while ((j < mink) .and. (i < mk))
      i = i+1
      call Add_Vector(n,j,Sub,Sub(1+(i-1)*n),Thr3)
    end do

    ! j should be mink, but who knows...

#   ifdef _DEBUGPRINT_
    if (j < mink) write(u6,'(2X,A,1X,I5)') 'Fewer vectors found:',j
#   endif
    Sub(j*n+1:) = Zero
    Ab(j*n+1:) = Zero
    call unitmat(EVec,maxk)
    mk = j
    old_mk = 0
    Augmented = .false.
    Reduced = .true.
    iter = iter-1

    ! Expand the subspace
    ! For each eigenpair i of the first k,
    ! check convergence for the residuals r:
    !  r = Ab * Vec(i) - Val(i) * Sub * Vec(i)
    ! Add a new vector, orthonormalized with the previous vectors,
    ! computed from r and the eigenpair
    ! (different possible variants)

  else
    call mma_allocate(Tmp,n,Label='Tmp')
    Conv = Zero
    jj = 0
    do i=0,k-1
      ! Vector in full space: Sub*Vec(i)
      call dGeMV_('N',n,mk,One,Sub,n,EVec(1+i*maxk),1,Zero,TVec,1)
      ! Product of matrix and vector: Ab*Vec(i)
      call dGeMV_('N',n,mk,One,Ab,n,EVec(1+i*maxk),1,Zero,TAV,1)
      ! Residual: (A-Val(i))*Vec(i) = Ab*Vec(i) - Val(i)*Sub*Vec(i)
      TRes(:) = TAV(:)-EVal(i+1)*TVec(:)
      Conv = max(Conv,DDot_(n,TRes,1,TRes,1))

      ! Scale vector, orthonormalize, and add to subspace

#     define DAV_METH DAV_IIGD
#     if DAV_METH == DAV_DPR
      ! Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
      do j=1,n
        Aux = A(nTri_Elem(j))-Eval(1+i)
        Diag(j) = One/sign(max(abs(Aux),Thr2),Aux)
      end do
      ! scale
      Tmp(:) = TRes(:)*Diag(:)
#     elif DAV_METH == DAV_IIGD
      ! Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
      do j=1,n
        Aux = A(nTri_Elem(j))-EVal(1+i)
        Diag(j) = One/sign(max(abs(Aux),Thr2),Aux)
      end do
      ! scale
      Tmp(:) = TRes(:)*Diag(:)
      Alpha = Zero
      do j=1,n
        Alpha = Alpha+Diag(j)*TVec(j)**2
      end do
      Alpha = DDot_(n,TVec,1,Tmp,1)/Alpha
      ! subtract
      Tmp(:) = Tmp(:)-Alpha*TVec(:)*Diag(:)
#     elif DAV_METH == DAV_GJD
      ! DO NOT USE THIS VARIANT!
      ! This is not practical as it stands, the equation should be
      ! solved only approximately, and it is not efficient to do this
      ! for each of the possibly many eigenpairs
      block
        integer(kind=iwp) :: iDum(1)
        real(kind=wp), allocatable :: P1(:)
        call mma_allocate(P1,n*n,Label='P1')
        ! project: (I-|v><v|) (A-e*I) (I-|v><v|) =
        !          A - e*I + (<v|P>+e)*|v><v| - |v><P| - |P><v|
        ! e = Val(i); |P> = A|v>
        Aux = DDot_(n,TVec,1,TAV,1)
        do kk=0,n-1
          ll = nTri_Elem(kk)
          do ii=0,kk
            P1(1+ii*n+kk) = A(ll+ii+1)+(Aux+EVal(1+i))*TVec(1+ii)*TVec(1+kk)-TAV(1+ii)*TVec(1+kk)-TAV(1+kk)*TVec(1+ii)
          end do
          P1(1+kk*n+kk) = P1(1+kk*n+kk)-EVal(1+i)
        end do
        iDum(1) = 0
        ! solve the equation
        call CG_Solver(n,n*n,P1,iDum,TRes,Tmp,info,5)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'CG iterations',info
#       endif
        call mma_deallocate(P1)
      end block
#     endif
      if (mk+jj <= n-1) then
        jj = mk+jj
        call Add_Vector(n,jj,Sub,Tmp,Thr3)
        jj = jj-mk
      end if
    end do
#   ifdef _DEBUGPRINT_
    write(u6,'(2X,A,1X,G12.6)') 'Maximum residual:',Conv
#   endif
    if ((Conv < Thr3) .and. (mk >= mink)) then
#     ifdef _DEBUGPRINT_
      write(u6,'(A)') 'Converged due to small residual'
#     endif
      Last = .true.
    else
      mk = min(mk+jj,n)

      ! If no new vector is found to add to the subspace, we are in trouble
      ! -Try to find a non-linear-dependent base vector in the original matrix

      if (jj == 0) then
#       ifdef _DEBUGPRINT_
        write(u6,'(A)') 'Process stagnated'
#       endif
        if (mk < maxk) then
          Tmp(:) = Zero
          i = 0
          do while ((jj < 1) .and. (i < n))
            i = i+1
            ig = mod(ig,n)+1
            ii = Indx(ig)
            Tmp(ii) = One
            jj = mk+jj
            call Add_Vector(n,jj,Sub,Tmp,Thr3)
            Tmp(ii) = Zero
            jj = jj-mk
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
    end if
    call mma_deallocate(Tmp)
    Reduced = .false.
  end if
end do
call mma_deallocate(Diag)
call mma_deallocate(TVec)
call mma_deallocate(TAV)
call mma_deallocate(TRes)
call mma_deallocate(Indx)

! Store the current lowest k eigenvectors (in the full space)
! Vec' = Sub * Vec(1:k)

call DGeMM_('N','N',n,k,mk,One,Sub,n,EVec,maxk,Zero,Vec,n)

call mma_deallocate(Sub)
call mma_deallocate(Ab)
call mma_deallocate(Proj)
call mma_deallocate(EVal)
call mma_deallocate(EVec)
call mma_deallocate(Eig_old)

end subroutine Davidson
