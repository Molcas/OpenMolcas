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
! Copyright (C) 2020, Oskar Weser                                      *
!***********************************************************************

module linalg_mod

#include "compiler_features.h"
#include "macros.fh"
#include "intent.fh"

use stdalloc, only: mma_allocate, mma_deallocate
use constants, only: Zero, One, cZero, cOne
use definitions, only: wp, iwp
#ifdef _ADDITIONAL_RUNTIME_CHECK_
use definitions, only: r8
#endif
use sorting, only: sort, argsort
use sorting_funcs, only: leq_i, leq_r, geq_r

implicit none
private

public :: mult, sym_diagonalize, is_close, operator(.isclose.), dot_product_, norm, canonicalize, Gram_Schmidt, Canonical, Lowdin, &
          symmetric
! TODO Move to different module
public :: abort_, verify_

type, abstract :: ParentCanonicalize_t
contains
  procedure :: canonicalize => base_canonicalize
  procedure(get_projections_t), deferred :: get_projections
end type

abstract interface
  subroutine get_projections_t(self,M,P)
    import :: wp, ParentCanonicalize_t
    class(ParentCanonicalize_t), intent(in) :: self
    real(kind=wp), intent(in) :: M(:,:)
    real(kind=wp), intent(out) :: P(:,:)
  end subroutine get_projections_t
end interface

type, extends(ParentCanonicalize_t) :: CanonicalBasisCanonicalize_t
contains
  procedure :: get_projections => project_canonical_unit_vectors
end type

type, extends(ParentCanonicalize_t) :: GeneralBasisCanonicalize_t
  real(kind=wp), pointer :: ref(:,:) => null()
contains
  procedure :: get_projections => project_general_unit_vectors
end type

!>  @brief
!>    Wrapper around dgemm.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  This is an overloaded wrapper around DGEMM which
!>  uses assumed shape arrays to deduce as much information
!>  as possible about the DGEMM calls automatically.
!>  The preferred way is to use arrays of actual rank 2
!>  instead of rank 1 arrays which are interpreted to be rank 2
!>  just by convention.
!>
!>  There are some run time checks to test for matching shapes.
interface mult
  module procedure :: mult_2D, multZ_2D, mult_2D_1D, mult_2d_raw
end interface mult

!>  @brief
!>  Canonicalize Eigenvectors and Eigenvalues
!>
!>  @details
!>  This procedure can be used to make diagonalizations deterministic,
!>  by forcing the Eigenvectors to be as similar as possible to proj_B.
!>  The algorithm is as following:
!>
!>  1. Sort the Eigenvalues (and the Eigenvectors accordingly).
!>  2. Find degenerate Eigenvalues (up to float equality) and put
!>      corresponding Eigenvectors into separate Eigenspaces.
!>  For each Eigenspace
!>      1. Let $b_i$ be the basis vectors in `proj_B` and $E$ your Eigenspace of dimension $d$.
!>      2. Calculate all projections of $b_i$ onto $E$. Let's call them $p_i$
!>      3. Stable-sort the $p_i$ by their norm.
!>      4. Take the first $d$ projections and normalize them. These are your canonical Eigenvectors.
!>
!>  @param[in,out] V 2D matrix which contains the Eigenvectors.
!>      The j-th column corresponds to the j-th Eigenvalue.
!>  @param[in,out] lambda 1D vector of Eigenvalues.
!>  @param[in] proj_B Optional and overloaded argument.
!>      If it is ommited, the canonical unit vector basis is assumed.
!>      Otherwise it can be 2D orthogonal matrix that represents the reference basis for
!>      the canonicalization or a function pointer of type `get_matrix_t`
!>      that returns matrix elements (trading speed for memory).
!>  @param[out] info Optional error code. If calling code does not check the error
!>      code, then this routine crashes upon errors.
interface canonicalize
  module procedure :: canonicalize_canonical_basis, canonicalize_general
end interface canonicalize

interface operator(.isclose.)
  module procedure :: isclose_for_operator
end interface

contains

!>  @brief
!>    Wrapper around dgemm for matrix-matrix multiplication.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @param[in] A
!>  @param[in] B
!>  @param[out] C The shape of the output array is usually
!>      [size(A, 1), size(B, 2)] which changes of course, if
!>      A or B are transposed.
!>  @param[in] transpA Optional argument to specify that A
!>      should be transposed.
!>  @param[in] transpB Optional argument to specify that B
!>      should be transposed.
subroutine mult_2D(A,B,C,transpA,transpB)
  real(kind=wp), intent(in) :: A(:,:), B(:,:)
  real(kind=wp), intent(out) :: C(:,:)
  logical(kind=iwp), intent(in), optional :: transpA, transpB
  logical(kind=iwp) :: transpA_, transpB_
  integer(kind=iwp) :: M, N, K_1, K
# ifdef _ADDITIONAL_RUNTIME_CHECK_
  integer(kind=iwp) :: K_2
# endif
  debug_function_name('mult_2D')

  if (present(transpA)) then
    transpA_ = transpA
  else
    transpA_ = .false.
  end if
  if (present(transpB)) then
    transpB_ = transpB
  else
    transpB_ = .false.
  end if

  M = size(A,merge(2,1,transpA_))
  ASSERT(M == size(C,1))
  N = size(B,merge(1,2,transpB_))
  ASSERT(N == size(C,2))
  K_1 = size(A,merge(1,2,transpA_))
# ifdef _ADDITIONAL_RUNTIME_CHECK_
  K_2 = size(B,merge(2,1,transpB_))
# endif
  ASSERT(K_1 == K_2)
  K = K_1

  ASSERT(wp == r8)
  call dgemm_(merge('T','N',transpA_),merge('T','N',transpB_),M,N,K,One,A,size(A,1),B,size(B,1),Zero,C,size(C,1))
end subroutine mult_2D

!>  @brief
!>    Wrapper around zgemm for matrix-matrix multiplication.
!>
!>  @author Vladislav Kochetov
!>
!>  @details
!>
!>  @param[in] A
!>  @param[in] B
!>  @param[out] C The shape of the output array is usually
!>      [size(A, 1), size(B, 2)] which changes of course, if
!>      A or B are conjugate transposed.
!>  @param[in] transpA Optional argument to specify that A
!>      should be conjugate transposed.
!>  @param[in] transpB Optional argument to specify that B
!>      should be conjugate transposed.
subroutine multZ_2D(A,B,C,transpA,transpB)
  complex(kind=wp), intent(in) :: A(:,:), B(:,:)
  complex(kind=wp), intent(out) :: C(:,:)
  logical(kind=iwp), intent(in), optional :: transpA, transpB
  logical(kind=iwp) :: transpA_, transpB_
  integer(kind=iwp) :: M, N, K_1, K
# ifdef _ADDITIONAL_RUNTIME_CHECK_
  integer(kind=iwp) :: K_2
# endif

  if (present(transpA)) then
    transpA_ = transpA
  else
    transpA_ = .false.
  end if
  if (present(transpB)) then
    transpB_ = transpB
  else
    transpB_ = .false.
  end if

  M = size(A,merge(2,1,transpA_))
  ASSERT(m == size(C,1))
  N = size(B,merge(1,2,transpB_))
  ASSERT(n == size(C,2))
  K_1 = size(A,merge(1,2,transpA_))
# ifdef _ADDITIONAL_RUNTIME_CHECK_
  K_2 = size(B,merge(2,1,transpB_))
# endif
  ASSERT(K_1 == K_2)
  K = K_1

  ASSERT(wp == r8)
  call zgemm_(merge('C','N',transpA_),merge('C','N',transpB_),M,N,K,cOne,A,size(A,1),B,size(B,1),cZero,C,size(C,1))
end subroutine multZ_2D

!>  @brief
!>    Wrapper around dgemm for matrix-vector multiplication.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @param[in] A
!>  @param[in] x
!>  @param[out] y The shape of the output array is size(x)
!>  @param[in] transpA Optional argument to specify that A
!>      should be transposed.
subroutine mult_2D_1D(A,x,y,transpA)
  real(kind=wp), intent(in) :: A(:,:), x(:)
  real(kind=wp), intent(out) :: y(:)
  logical(kind=iwp), intent(in), optional :: transpA
  logical(kind=iwp) :: transpA_
  integer(kind=iwp) :: M, N, K

  if (present(transpA)) then
    transpA_ = transpA
  else
    transpA_ = .false.
  end if

  M = size(A,merge(1,2,.not. transpA_))
  ASSERT(M == size(y,1))
  N = 1
  K = size(A,merge(2,1,.not. transpA_))
  ASSERT(M == size(x,1))

  ASSERT(wp == r8)
  call dgemm_(merge('T','N',transpA_),'N',M,N,K,One,A,size(A,1),x,size(x,1),Zero,y,size(y,1))
end subroutine mult_2D_1D

!>  @brief
!>    Wrapper around dgemm for matrix-matrix multiplication with raw memory.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is important to note, that it is expected to pass in the
!>  actual vectors and not just the first element.
!>  So if `A_ptr` is the pointer to A in the work array
!>  it is necessary to call it with `Work(A_ptr : )`.
!>
!>  @param[in] A
!>  @param[in] shapeA The shape of A.
!>  @param[in] B
!>  @param[in] shapeB The shape of B.
!>  @param[out] C The shape of the output array is usually
!>      (shapeA(1) * shapeB(2)) which changes of course, if
!>      A or B are transposed.
!>  @param[in] transpA Optional argument to specify that A
!>      should be transposed.
!>  @param[in] transpB Optional argument to specify that B
!>      should be transposed.
subroutine mult_2D_raw(A,shapeA,B,shapeB,C,transpA,transpB)
  real(kind=wp), intent(in), target :: A(*)
  integer(kind=iwp), intent(in) :: shapeA(2)
  real(kind=wp), intent(in), target :: B(*)
  integer(kind=iwp), intent(in) :: shapeB(2)
  real(kind=wp), intent(_OUT_), target :: C(*)
  logical(kind=iwp), intent(in), optional :: transpA, transpB
  logical(kind=iwp) :: transpA_, transpB_
  integer(kind=iwp) :: shapeC(2)
  real(kind=wp), pointer :: ptr_A(:,:), ptr_B(:,:), ptr_C(:,:)

  if (present(transpA)) then
    transpA_ = transpA
  else
    transpA_ = .false.
  end if
  if (present(transpB)) then
    transpB_ = transpB
  else
    transpB_ = .false.
  end if

  shapeC = [merge(shapeA(1),shapeA(2),.not. transpA_),merge(shapeB(2),shapeB(1),.not. transpB_)]

  ptr_A(1:shapeA(1),1:shapeA(2)) => A(1:product(shapeA))
  ptr_B(1:shapeB(1),1:shapeB(2)) => B(1:product(shapeB))
  ptr_C(1:shapeC(1),1:shapeC(2)) => C(1:product(shapeC))

  call mult(ptr_A,ptr_B,ptr_C,transpA_,transpB_)
end subroutine mult_2D_raw

!>  @brief
!>  Diagonalize symmetric A.
!>
!>  @details
!>  Diagonalize A to fullfill `A V(:, j) = lambda(j) V(:, j)`.
!>  Wrapper around BLAS DSYEV.
!>
!>  @param[in] A Symmetric 2D matrix to be diagonalized.
!>  @param[out] V 2D matrix which contains the Eigenvectors.
!>      The j-th column corresponds to the j-th Eigenvalue.
!>  @param[out] lambda 1D vector of Eigenvalues.
!>  @param[in] info Error code. If calling code does not check the error
!>      code, then this routine crashes upon errors.
subroutine sym_diagonalize(A,V,lambda,info)
  real(kind=wp), intent(in) :: A(:,:)
  real(kind=wp), intent(out) :: V(:,:), lambda(:)
  integer(kind=iwp), optional, intent(out) :: info
  integer(kind=iwp), parameter :: do_worksize_query = -1
  integer(kind=iwp) :: info_
  real(kind=wp), allocatable :: work(:)
  real(kind=wp) :: dummy(2), query_result(2)

  ASSERT(symmetric(A))

  V(:,:) = A(:,:)
  info_ = 0
  call dsyev_('V','L',size(V,2),dummy,size(V,1),dummy,query_result,do_worksize_query,info_)

  if (info_ /= 0) then
    if (present(info)) then
      info = info_
    else
      call abort_('Error in diagonalize')
    end if
  end if

  call mma_allocate(work,int(query_result(1)))
  info_ = 0
  call dsyev_('V','L',size(V,2),V,size(V,1),lambda,work,size(work),info_)
  call mma_deallocate(work)

  if (info_ /= 0) then
    if (present(info)) then
      info = info_
    else
      call abort_('Error in diagonalize')
    end if
  end if
end subroutine sym_diagonalize

!>  @brief
!>  Order Eigenvectors by ascending eigenvalues.
subroutine order_eigenvectors(V,lambda)
  real(kind=wp), intent(inout) :: V(:,:), lambda(:)
  integer(kind=iwp), allocatable :: idx(:)
  call mma_allocate(idx,size(lambda))
  idx(:) = argsort(lambda,leq_r)
  V(:,:) = V(:,idx)
  lambda(:) = lambda(idx)
  call mma_deallocate(idx)
end subroutine order_eigenvectors

!>  @brief
!>  Determine the Eigenspaces.
!>
!>  @details
!>  Sort the eigenvalues and collect "equal" eigenvalues together.
!>  Equality of eigenvalues considers floating point error and is determined
!>  by `operator(.isclose.)`.
!>  For this reason all Eigenvalues of the same Eigenspace are replaced
!>  with their mean.
!>
!>  @param[in,out] lambda Eigenvalues are sorted ascendingly and
!>      Eigenvalues of the same Eigenspace (up to floating point error)
!>      are replaced with their mean.
!>  @param[out] dimensions The dimension of each Eigenspace.
!>      size(dimensions) is the number of Eigenspaces, i.e. the number
!>      of distinct Eigenvalues and sum(dimensions) is the number of Eigenvalues.
subroutine determine_eigenspaces(lambda,dimensions)
  real(kind=wp), intent(inout) :: lambda(:)
  integer(kind=iwp), allocatable, intent(out) :: dimensions(:)
  integer(kind=iwp), allocatable :: d_buffer(:)
  integer(kind=iwp) :: i, low
  integer(kind=iwp) :: n_spaces

# ifdef _ADDITIONAL_RUNTIME_CHECK_
  if (any(lambda(2:) < lambda(:size(lambda)-1))) then
    call abort_('Eigenvalues not sorted in'//__FILE__)
  end if
# endif

  call mma_allocate(d_buffer,size(lambda))
  d_buffer(:) = 0

  low = 1
  n_spaces = 1
  do i=1,size(lambda)
    d_buffer(n_spaces) = d_buffer(n_spaces)+1
    if (i+1 <= size(lambda)) then
      if (.not. is_close(lambda(i),lambda(i+1),epsilon(lambda)*1.0e3_wp,1.0e-8_wp)) then
        n_spaces = n_spaces+1
        lambda(low:i) = mean(lambda(low:i))
        low = i+1
      end if
    else
      lambda(low:i) = mean(lambda(low:i))
    end if
  end do

  call mma_allocate(dimensions,n_spaces)
  dimensions(:) = d_buffer(:n_spaces)
  call mma_deallocate(d_buffer)

contains
  pure function mean(X)
    real(kind=wp) :: mean
    real(kind=wp), intent(in) :: X(:)
    mean = sum(X)/size(X)
  end function mean
end subroutine determine_eigenspaces

subroutine base_canonicalize(self,V,lambda)
  class(ParentCanonicalize_t), intent(in) :: self
  real(kind=wp), intent(inout) :: V(:,:), lambda(:)
  !> The norms of the projections
  real(kind=wp), allocatable :: norms_projections(:)
  integer(kind=iwp) :: low, i, j, d, n_new
  integer(kind=iwp), allocatable :: idx(:), dimensions(:)
  real(kind=wp), allocatable :: projections(:,:), ONB(:,:)
  debug_function_name('canonicalize_factory')

  ASSERT(size(V,1) == size(V,2))
  call mma_allocate(projections,size(V,1),size(V,2))
  call mma_allocate(norms_projections,size(projections,2))
  call mma_allocate(idx,size(V,1))

  call order_eigenvectors(V,lambda)
  do j=1,size(V,2)
    V(:,j) = V(:,j)/norm(V(:,j))
  end do

  call determine_eigenspaces(lambda,dimensions)

  call mma_allocate(ONB,size(lambda),maxval(dimensions))

  low = 1
  do j=1,size(dimensions)
    d = dimensions(j)
    ! `self` "knows" which vectors are to be projected
    !  into V(:, low : low + d - 1)
    call self%get_projections(V(:,low:low+d-1),projections)

    do i=1,size(norms_projections)
      norms_projections(i) = norm(projections(:,i))
    end do
    idx(:) = argsort(norms_projections,geq_r)
    call sort(idx(:d),leq_i)

    ! Get the first d projections with the largest overlap.
    do i=1,d
      V(:,low+i-1) = projections(:,idx(i))/norms_projections(idx(i))
    end do

    ! reorthogonalize
    call Gram_Schmidt(V(:,low:low+d-1),d,ONB(:,:d),n_new)
    ASSERT(d == n_new)
    V(:,low:low+d-1) = ONB(:,:d)

    low = low+d
  end do

  call mma_deallocate(norms_projections)
  call mma_deallocate(idx)
  call mma_deallocate(ONB)
  call mma_deallocate(projections)
  call mma_deallocate(dimensions)
end subroutine base_canonicalize

subroutine project_canonical_unit_vectors(self,M,P)
  class(CanonicalBasisCanonicalize_t), intent(in) :: self
  real(kind=wp), intent(in) :: M(:,:)
  real(kind=wp), intent(out) :: P(:,:)
  integer(kind=iwp) :: i, j
  unused_var(self)

  P(:,:) = Zero
  ! calculate projections of basis b_j into eigenvectors v_i:
  ! p_j = sum_i < b_j | v_i > v_i
  ! if b_j are canonical unit vectors, the dot product is the j-th
  ! component of v_i
  do j=1,size(P,2)
    do i=1,size(M,2)
      P(:,j) = P(:,j)+M(j,i)*M(:,i)
    end do
  end do
end subroutine project_canonical_unit_vectors

subroutine project_general_unit_vectors(self,M,P)
  class(GeneralBasisCanonicalize_t), intent(in) :: self
  real(kind=wp), intent(in) :: M(:,:)
  real(kind=wp), intent(out) :: P(:,:)
  integer(kind=iwp) :: i, j

  P(:,:) = Zero
  ! calculate projections of basis b_j into eigenvectors v_i:
  ! p_j = sum_i < b_j | M_i > M_i
  do j=1,size(P,2)
    do i=1,size(M,2)
      P(:,j) = P(:,j)+dot_product_(self%ref(:,j),M(:,i))*M(:,i)
    end do
  end do
end subroutine project_general_unit_vectors

subroutine canonicalize_canonical_basis(V,lambda)
  real(kind=wp), intent(inout) :: V(:,:), lambda(:)
  type(CanonicalBasisCanonicalize_t) :: canonicalizer

  call canonicalizer%canonicalize(V,lambda)
end subroutine canonicalize_canonical_basis

subroutine canonicalize_general(V,lambda,ref)
  real(kind=wp), intent(inout) :: V(:,:), lambda(:)
  real(kind=wp), intent(in), target :: ref(:,:)
  type(GeneralBasisCanonicalize_t) :: canonicalizer

  canonicalizer%ref => ref
  call canonicalizer%canonicalize(V,lambda)
end subroutine canonicalize_general

!>  @brief
!>  Check for floating point equality.
!>
!>  @details
!>  Checks if \f[ | a - b | \leq max(rtol * max(|a|, |b|), atol)  \f].
!>
!>  @param[in] a A real number.
!>  @param[in] b A real number.
!>  @param[in] atol Absolute tolerance.
!>  @param[in] rtol Relative tolerance.
!>
!>  @author
!>    Oskar Weser
! Change name from isclose to is_close, for nvfortran's sake
! (https://forums.developer.nvidia.com/t/nvfortran-wrong-error-compiler-is-confused-by-similar-names/186146)
elemental function is_close(a,b,atol,rtol) result(res)
  real(kind=wp), intent(in) :: a, b
  real(kind=wp), intent(in) :: atol, rtol
  logical(kind=iwp) :: res
  res = abs(a-b) <= max(rtol*max(abs(a),abs(b)),atol)
end function is_close

! Operator functions may only have two arguments.
elemental function isclose_for_operator(a,b) result(res)
  real(kind=wp), intent(in) :: a, b
  logical(kind=iwp) :: res

  ! For real64, the epsilon is 2.3e-16 so that atol becomes roughly 2.3e-14 which is very tight.
  res = is_close(a,b,atol=epsilon(a)*1.0e2_wp,rtol=1.0e-9_wp)
end function isclose_for_operator

!>  @brief
!>    Calculates v1^T S v2.
!>
!>  @details
!>  Calculates \f[ v1^T S v2 \f]
!>  S has to be a symmetric positive definite matrix, which is not
!>  tested.
!>  If S is ommited, it defaults to the unit matrix,
!>  i.e. the Euclidean dot-product.
!>
!>  @param[in] v1 A vector.
!>  @param[in] v2 A vector.
!>  @param[in] S Optional positive definite overlap matrix. If ommited
!>      it is assumed to be the unit matrix.
!>
!>  @author
!>    Oskar Weser
function dot_product_(v1,v2,S) result(dot)
  real(kind=wp), intent(in) :: v1(:), v2(:)
  real(kind=wp), intent(in), optional :: S(:,:)
  real(kind=wp), allocatable :: tmp(:)
  real(kind=wp) :: dot

  if (present(S)) then
    call mma_allocate(tmp,size(v2))
    call mult(S,v2,tmp)
    dot = dot_product(v1,tmp)
    call mma_deallocate(tmp)
  else
    dot = dot_product(v1,v2)
  end if
end function dot_product_

!>  @brief
!>    The induced norm from the dot product given by S.
!>
!>  @details
!>    It is assumed, that S is symmetric and positive definite.
!>
!>  @param[in] v A vector.
!>  @param[in] S Optional positive definite overlap matrix. If ommited
!>      it is assumed to be the unit matrix.
!>
!>  @author
!>    Oskar Weser
function norm(v,S) result(L)
  real(kind=wp), intent(in) :: v(:)
  real(kind=wp), intent(in), optional :: S(:,:)
  real(kind=wp) :: L
  if (present(S)) then
    L = sqrt(dot_product_(v,v,S))
  else
    ! One could use norm2 here, but Sun and PGI compilers don't know this.
    L = sqrt(sum(v**2))
  end if
end function norm

!>  @brief
!>  Check if matrix is symmetric.
!>
!>  @param[in] M A matrix.
pure function symmetric(M)
  logical(kind=iwp) :: symmetric
  real(kind=wp), intent(in) :: M(:,:)
  integer(kind=iwp) :: i, j

  symmetric = .true.
  do j=1,size(M,2)
    do i=j,size(M,1)
      if (.not. (M(i,j) .isclose. M(j,i))) then
        symmetric = .false.
        return
      end if
    end do
  end do
end function symmetric

!>  @brief
!>  Loewdin-orthonormalize basis to get ONB using the overlap matrix S.
!>
!>  @details
!>  If the given basis is linear dependent, then this routine will abort.
!>    For a detailed explanation see \cite szabo_ostlund (p. 143).
!>
!>  @author
!>    Oskar Weser
!>
!>  @param[in] basis A matrix of full rank.
!>  @param[out] ONB The orthonormalized version of basis.
!>  @param[in] S Optional positive definite overlap matrix. If ommited
!>      it is assumed to be the unit matrix.
subroutine Lowdin(basis,ONB,S)
  real(kind=wp), intent(in) :: basis(:,:)
  real(kind=wp), intent(out) :: ONB(:,:)
  real(kind=wp), intent(in), optional :: S(:,:)
  integer(kind=iwp) :: i
  real(kind=wp), allocatable :: U(:,:), s_diag(:), X(:,:), S_transf(:,:), tmp(:,:)

  call mma_allocate(S_transf,size(S,1),size(S,2))
  call mma_allocate(U,size(S,1),size(S,2))
  call mma_allocate(X,size(S,1),size(S,2))
  call mma_allocate(tmp,size(S,1),size(S,2))
  call mma_allocate(s_diag,size(S,2))

  ! Transform AO-overlap matrix S to the overlap matrix of basis.
  ! S_transf = basis^T S basis
  ! We search X that diagonalizes S_transf
  if (present(S)) then
    call mult(S,basis,tmp)
    call mult(basis,tmp,S_transf,transpA=.true.)
  else
    call mult(basis,basis,S_transf,transpA=.true.)
  end if

  call sym_diagonalize(S_transf,U,s_diag)
  call canonicalize(U,s_diag)

  call verify_(all(s_diag > 1.0e-10_wp), 'Linear dependency detected. Lowdin can''t cure it. '// &
                                         'Please use Gram_Schmidt or Canonical orthonormalization.')

  ! X = U s_diag^{-1/2} U^T
  do i=1,size(tmp,2)
    tmp(:,i) = U(:,i)/sqrt(s_diag(i))
  end do
  call mult(tmp,U,X,transpB=.true.)
  ! With this X the overlap matrix S_transf has diagonal form.
  ! X^T basis^T S basis X = 1
  ! We finally have to convert to get the form:
  ! ONB^T S ONB = 1
  call mult(basis,X,ONB)

  call mma_deallocate(tmp)
  call mma_deallocate(s_diag)
  call mma_deallocate(X)
  call mma_deallocate(U)
  call mma_deallocate(S_transf)
end subroutine Lowdin

!>  @brief
!>  Canonical-orthonormalize basis to get ONB using the overlap matrix S.
!>
!>  @details
!>  If the given basis is linear dependent, then this routine will still
!>  work. The out parameter `n_new` returns the number of valid
!>  column vectors. A valid orthonormalized basis is given by ONB(:, : n_new).
!>
!>  For a detailed explanation see \cite szabo_ostlund (p. 143).
!>
!>  @author
!>    Oskar Weser
!>
!>  @param[in] basis A matrix of full rank.
!>  @param[in] n_to_ON Number of columns that should be orthonormalized.
!>  @param[out] ONB The orthonormalized version of basis.
!>  @param[out] n_new The number of columns that could be orthonormalized.
!>     Note that `n_new <= n_to_ON`.
!>  @param[in] S Optional positive definite overlap matrix. If ommited
!>      it is assumed to be the unit matrix.
subroutine Canonical(basis,n_to_ON,ONB,n_new,S)
  real(kind=wp), intent(in) :: basis(:,:)
  integer(kind=iwp), intent(in) :: n_to_ON
  real(kind=wp), intent(out) :: ONB(:,:)
  integer(kind=iwp), intent(out) :: n_new
  real(kind=wp), intent(in), optional :: S(:,:)
  logical(kind=iwp) :: lin_dep_detected
  integer(kind=iwp) :: i
  integer(kind=iwp), allocatable :: idx(:)
  real(kind=wp), allocatable :: U(:,:), s_diag(:), S_transf(:,:), X(:,:), tmp(:,:)

  call mma_allocate(S_transf,size(S,1),size(S,2))
  call mma_allocate(U,size(S,1),size(S,2))
  call mma_allocate(s_diag,size(S,2))
  call mma_allocate(X,size(S,1),size(S,2))
  call mma_allocate(idx,size(S,1))
  call mma_allocate(tmp,size(S,1),size(S,2))

  ! Transform AO-overlap matrix S to the overlap matrix of basis.
  ! We search X that diagonalizes basis^T S basis
  if (present(S)) then
    call mult(S,basis,tmp)
    call mult(basis,tmp,S_transf,transpA=.true.)
  else
    call mult(basis,basis,S_transf,transpA=.true.)
  end if

  call sym_diagonalize(S_transf,U,s_diag)
  call canonicalize(U,s_diag)

  idx(:) = argsort(s_diag,geq_r)
  U(:,:) = U(:,idx)
  s_diag(:) = s_diag(idx)

  i = 0
  lin_dep_detected = .false.
  do while ((.not. lin_dep_detected) .and. (i < n_to_ON))
    if (s_diag(i+1) < 1.0e-10_wp) then
      n_new = i
      lin_dep_detected = .true.
    end if
    i = i+1
  end do
  if (.not. lin_dep_detected) n_new = n_to_ON

  ! X = U s_diag^{-1/2}
  do i=1,n_new
    X(:,i) = U(:,i)/sqrt(s_diag(i))
  end do
  ! With this X the overlap matrix S_transf has diagonal form.
  ! X^T basis^T S basis X = 1
  ! We finally have to convert to get the form:
  ! ONB^T S ONB = 1
  ONB(:,n_new+1:) = basis(:,n_new+1:)
  call mult(basis,X(:,:n_new),ONB(:,:n_new))

  call mma_deallocate(tmp)
  call mma_deallocate(X)
  call mma_deallocate(idx)
  call mma_deallocate(s_diag)
  call mma_deallocate(U)
  call mma_deallocate(S_transf)
end subroutine Canonical

!>  @brief
!>  Gram-Schmidt-orthonormalize basis to get ONB using the overlap matrix S.
!>
!>  @details
!>  If the given basis is linear dependent, then this routine will still
!>  work. The out parameter `n_new` returns the number of valid
!>  column vectors. A valid orthonormalized basis is given by ONB(:, : n_new).
!>
!>  For a detailed explanation see \cite szabo_ostlund (p. 143).
!>
!>  @author
!>    Oskar Weser
!>
!>  @param[in] basis A matrix of full rank.
!>  @param[in] n_to_ON Number of columns that should be orthonormalized.
!>  @param[out] ONB The orthonormalized version of basis.
!>  @param[out] n_new The number of columns that could be orthonormalized.
!>     Note that `n_new <= n_to_ON`.
!>  @param[in] S Optional positive definite overlap matrix. If ommited
!>      it is assumed to be the unit matrix.
subroutine Gram_Schmidt(basis,n_to_ON,ONB,n_new,S)
  real(kind=wp), intent(in) :: basis(:,:)
  integer(kind=iwp), intent(in) :: n_to_ON
  real(kind=wp), target, intent(out) :: ONB(:,:)
  integer(kind=iwp), intent(out) :: n_new
  real(kind=wp), intent(in), optional :: S(:,:)
  real(kind=wp) :: L
  integer(kind=iwp) :: i, j
  logical(kind=iwp) :: lin_dep_detected, improve_solution
  real(kind=wp), allocatable :: previous(:), correction(:), v(:)
  real(kind=wp), pointer :: curr(:)

  call mma_allocate(previous,size(basis,1))
  call mma_allocate(correction,size(basis,1))
  call mma_allocate(v,size(basis,1))

  n_new = 0
  ONB(:,n_to_ON+1:) = basis(:,n_to_ON+1:)
  do i=1,n_to_ON
    curr => ONB(:,n_new+1)
    curr = basis(:,i)

    improve_solution = .true.
    lin_dep_detected = .false.
    do while (improve_solution .and. (.not. lin_dep_detected))
      correction = Zero
      if (present(S)) then
        call mult(S,curr,v)
      else
        v(:) = curr(:)
      end if

      do j=1,n_new
        correction(:) = correction(:)+ONB(:,j)*dot_product(ONB(:,j),v)
      end do
      curr = curr-correction
      improve_solution = norm(correction,S=S) > 0.1_wp
      L = norm(curr,S=S)
      lin_dep_detected = L < 1.0e-10_wp
      if (.not. lin_dep_detected) then
        curr = curr/L
      end if
      if (.not. (improve_solution .or. lin_dep_detected)) then
        n_new = n_new+1
      end if
    end do
  end do
  ONB(:,n_new+1:n_to_ON) = basis(:,n_new+1:n_to_ON)

  call mma_deallocate(v)
  call mma_deallocate(correction)
  call mma_deallocate(previous)
end subroutine Gram_Schmidt

!> @brief
!>    Print error message and abort.
subroutine abort_(message)
  character(len=*), intent(in) :: message
  call WarningMessage(2,message)
  call pure_abort()
end subroutine abort_

!> @brief
!>    Abort.
pure subroutine pure_abort()
  ! I know, that these functions are not pure, but we are aborting anyway.
  interface
    pure subroutine Abend()
    end subroutine Abend
  end interface
  call Abend()
end subroutine pure_abort

!> @brief
!>    Runtime check, that is not switched off in Debug mode
subroutine verify_(test_expression,message)
  logical(kind=iwp), intent(in) :: test_expression
  character(len=*), intent(in) :: message
  if (.not. test_expression) then
    call abort_(message)
  end if
end subroutine verify_

end module linalg_mod
