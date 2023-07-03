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

implicit none
integer n, k, iRC
real*8 A(n*(n+1)/2), Eig(k), Vec(n,k)
real*8, dimension(:), allocatable :: Eig_old, EVec, Sub, Ab, Proj, EVal
real*8 Aux, Thr, Thr2, Thr3, Conv, Alpha
real*8 ddot_
integer mk, old_mk, mink, maxk, ig, info, nTmp, iter, maxiter
integer i, j, ii, jj
logical Last, Augmented, Reduced
external ddot_
parameter(Thr=1.0D-7,maxiter=300,Thr2=1.0D-16,Thr3=1.0D-16)
real*8 rDum(1)
real*8, allocatable :: Vec2(:), Val(:), Tmp(:), Diag(:), TVec(:), TAV(:), TRes(:)
integer, allocatable :: index(:)
#include "stdalloc.fh"
#include "real.fh"

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

if (k > n) then
  call SysAbendMsg('Davidson','Wrong k value.','')
end if
mink = min(max(k+2,5),n)
maxk = min(5*mink,n)
mk = k
iRC = 0

! If all the eigenvalues are wanted, better solve the system directly
! and return

if (mk >= n) then
  call mma_allocate(Val,n,Label='Val')
  call mma_allocate(Vec2,n*n,Label='Vec2')
  call FZero(Vec,n*n)
  do j=1,n
    jj = (j-1)*n
    do i=1,j
      Vec2(jj+i) = A(j*(j-1)/2+i)
    end do
  end do
  call dsyev_('V','U',n,Vec2,n,Val,rDum,-1,info)
  nTmp = int(rDum(1))

  call mma_allocate(Tmp,nTmp,Label='Tmp')
  call dsyev_('V','U',n,Vec2,n,Val,Tmp,nTmp,info)
  call SortEig(Val,Vec2,n,n,1,.false.)
  call mma_deallocate(Tmp)

  call dcopy_(k,Val,1,Eig,1)
  call dcopy_(n*k,Vec2,1,Vec,1)

  call mma_deallocate(Val)
  call mma_deallocate(Vec2)
# ifdef _DEBUGPRINT_
  call RecPrt('Eigenvalues',' ',Eig,1,n)
  write(6,*)
  write(6,'(A)') 'Complete system solved'
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

call mma_allocate(Index,n,Label='Index')
do i=1,n
  index(i) = i
end do
do i=1,n
  ii = index(i)
  Aux = A(ii*(ii+1)/2)
  ii = i
  do j=i,n
    jj = index(j)
    if (A(jj*(jj+1)/2) < Aux) then
      Aux = A(jj*(jj+1)/2)
      ii = j
    end if
  end do
  if (ii /= i) then
    jj = index(ii)
    index(ii) = index(i)
    index(i) = jj
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
  call dcopy_(n,Vec(1,i),1,Tmp,1)
  call Add_Vector(n,nTmp,Sub,Tmp,Thr3)
end do

ii = 0
Tmp(:) = Zero
do while ((nTmp < mk) .and. (ii < n))
  ii = ii+1
  jj = index(ii)
  Tmp(jj) = One
  call Add_Vector(n,nTmp,Sub,Tmp,Thr3)
  Tmp(jj) = Zero
end do
! ig will be a global counter to loop across all n base vectors
ig = ii
call mma_deallocate(Tmp)
call FZero(Sub(1+mk*n),(maxk-mk)*n)

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
  if (iter > 1) call dcopy_(k,Eig,1,Eig_old,1)
# ifdef _DEBUGPRINT_
  if (.not. Reduced) then
    write(6,'(A)') '---------------'
    write(6,'(A,1X,I5)') 'Iteration',iter
  end if
  call RecPrt('Orthonormalized subspace',' ',Sub,n,mk)
# endif

  ! Compute the matrix product
  ! Ab = A * Sub
  ! Only the new vectors since the last iterations need to be calculated

  call mma_allocate(Tmp,n,Label='Tmp')
  do i=1,n
    ! Reconstruct a row (or column) of the matrix A
    do j=1,i
      Tmp(j) = A(i*(i-1)/2+j)
    end do
    do j=i+1,n
      Tmp(j) = A(j*(j-1)/2+i)
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
    write(6,'(2X,A,1X,I5)') 'Solving for subspace size:',mk
#   endif
    call dcopy_(maxk*maxk,Proj,1,EVec,1)
    call dsyev_('V','L',mk,EVec,maxk,EVal,rDum,-1,info)
    nTmp = int(rDum(1))
    call mma_allocate(Tmp,nTmp,Label='Tmp')
    call dsyev_('V','L',mk,EVec,maxk,EVal,Tmp,nTmp,info)
    call mma_deallocate(Tmp)
    call SortEig(EVal,EVec,mk,maxk,1,.false.)
    call dcopy_(k,EVal,1,Eig,1)
#   ifdef _DEBUGPRINT_
    call RecPrt('Current guess',' ',Eig,1,k)
#   endif
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('Eigenvalues',' ',EVal,1,mk)
  call SubRecPrt('Subspace Eigenvectors',' ',EVec,maxk,mk,mk)
  write(6,*)
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
    if (Augmented) write(6,'(2X,A,1X,G12.6)') 'Maximum relative eigenvalue change:',Conv
#   endif
  else
    Conv = Ten*Thr
  end if
  old_mk = mk
  if (Augmented .and. (Conv <= Thr) .and. (mk >= mink)) then
#   ifdef _DEBUGPRINT_
    write(6,'(A)') 'Converged due to small change'
#   endif
    Last = .true.
  else if (mk == n) then
#   ifdef _DEBUGPRINT_
    write(6,'(A)') 'Complete system solved'
#   endif
    Last = .true.
  else if (iter >= maxiter) then
#   ifdef _DEBUGPRINT_
    write(6,'(A)') 'Not converged'
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
    write(6,'(2X,A,1X,I5)') 'Reducing search space to',mink
#   endif
    call mma_allocate(Tmp,mink*n,Label='Tmp')
    call DGeMM_('N','N',n,mink,mk,One,Sub,n,EVec,maxk,Zero,Tmp,n)
    call dcopy_(mink*n,Tmp,1,Sub,1)
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
    if (j < mink) then
      write(6,'(2X,A,1X,I5)') 'Fewer vectors found:',j
    end if
#   endif
    call FZero(Sub(1+j*n),(maxk-j)*n)
    call FZero(Ab(1+j*n),(maxk-j)*n)
    call FZero(EVec,maxk*maxk)
    do i=0,j-1
      EVec(1+i*(maxk+1)) = One
    end do
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
      call dcopy_(n,TAV,1,TRes,1)
      call daxpy_(n,-EVal(1+i),TVec,1,TRes,1)
      Conv = max(Conv,DDot_(n,TRes,1,TRes,1))

      ! Scale vector, orthonormalize, and add to subspace

#     define DAV_METH DAV_IIGD
#     if DAV_METH == DAV_DPR
      ! Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
      do j=1,n
        Aux = A(j*(j+1)/2)-Eval(1+i)
        Diag(j) = One/sign(max(abs(Aux),Thr2),Aux)
      end do
      ! scale
      do j=1,n
        Tmp(j) = TRes(j)*Diag(j)
      end do
#     elif DAV_METH == DAV_IIGD
      ! Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
      do j=1,n
        Aux = A(j*(j+1)/2)-EVal(1+i)
        Diag(j) = One/sign(max(abs(Aux),Thr2),Aux)
      end do
      ! scale
      do j=1,n
        Tmp(j) = TRes(j)*Diag(j)
      end do
      Alpha = Zero
      do j=1,n
        Alpha = Alpha+Diag(j)*TVec(j)**2
      end do
      Alpha = DDot_(n,TVec,1,Tmp,1)/Alpha
      ! subtract
      do j=1,n
        TVec(j) = TVec(j)*Diag(j)
      end do
      call daxpy_(n,-Alpha,TVec,1,Tmp,1)
#     elif DAV_METH == DAV_GJD
      ! DO NOT USE THIS VARIANT!
      ! This is not practical as it stands, the equation should be
      ! solved only approximately, and it is not efficient to do this
      ! for each of the possibly many eigenpairs
      block
        integer iDum(1)
        real*8, allocatable :: P1(:)
        call mma_allocate(P1,n*n,Label='P1')
        ! project: (I-|v><v|) (A-e*I) (I-|v><v|) =
        !          A - e*I + (<v|P>+e)*|v><v| - |v><P| - |P><v|
        ! e = Val(i); |P> = A|v>
        Aux = DDot_(n,TVec,1,TAV,1)
        do kk=0,n-1
          ll = kk*(kk+1)/2
          do ii=0,kk
            P1(1+ii*n+kk) = A(ll+ii+1)+(Aux+EVal(1+i))*TVec(1+ii)*TVec(1+kk)-TAV(1+ii)*TVec(1+kk)-TAV(1+kk)*TVec(1+ii)
          end do
          P1(1+kk*n+kk) = P1(1+kk*n+kk)-EVal(1+i)
        end do
        iDum(1) = 0
        ! solve the equation
        call CG_Solver(n,n*n,P1,iDum,TRes,Tmp,info,5)
#       ifdef _DEBUGPRINT_
        write(6,*) 'CG iterations',info
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
    write(6,'(2X,A,1X,G12.6)') 'Maximum residual:',Conv
#   endif
    if ((Conv < Thr3) .and. (mk >= mink)) then
#     ifdef _DEBUGPRINT_
      write(6,'(A)') 'Converged due to small residual'
#     endif
      Last = .true.
    else
      mk = min(mk+jj,n)

      ! If no new vector is found to add to the subspace, we are in trouble
      ! -Try to find a non-linear-dependent base vector in the original matrix

      if (jj == 0) then
#       ifdef _DEBUGPRINT_
        write(6,'(A)') 'Process stagnated'
#       endif
        if (mk < maxk) then
          call FZero(Tmp,n)
          i = 0
          do while ((jj < 1) .and. (i < n))
            i = i+1
            ig = mod(ig,n)+1
            ii = index(ig)
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
call mma_deallocate(Index)

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
