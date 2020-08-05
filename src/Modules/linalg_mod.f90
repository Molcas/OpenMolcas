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
    use stdalloc, only: mma_allocate, mma_deallocate
    use definitions, only: wp, r8
    use sorting, only: sort, argsort
    use sorting_funcs, only : integer_leq => le_i, ge_r
    implicit none
    private
    public :: mult, diagonalize, isclose, operator(.isclose.), &
        dot_product_, norm, canonicalize, Gram_Schmidt, Canonical, Lowdin, &
        symmetric

    ! TODO Move to different module
    public :: assert_, abort_, assume_


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
        module procedure mult_2D, mult_2D_1D, mult_2d_raw
    end interface

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
!>  @param[inout] V 2D matrix which contains the Eigenvectors.
!>      The j-th column corresponds to the j-th Eigenvalue.
!>  @param[inout] lambda 1D vector of Eigenvalues.
!>  @param[in] proj_B Optional and overloaded argument.
!>      If it is ommited, the canonical unit vector basis is assumed.
!>      Otherwise it can be 2D orthogonal matrix that represents the reference basis for
!>      the canonicalization or a function pointer of type `get_matrix_t`
!>      that returns matrix elements (trading speed for memory).
!>  @param[in] info Error code. If calling code does not check the error
!>      code, then this routine crashes upon errors.
    interface canonicalize
        module procedure canonicalize_default
    end interface

    interface operator(.isclose.)
        module procedure isclose_for_operator
    end interface

contains

!>  @brief
!>    Wrapper around dgemm for matrix-matrix multiplication.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @paramin[in] A
!>  @paramin[in] B
!>  @paramin[out] C The shape of the output array is usually
!>      [size(A, 1), size(B, 2)] which changes of course, if
!>      A or B are transposed.
!>  @paramin[in] transpA, Optional argument to specify that A
!>      should be transposed.
!>  @paramin[in] transpB, Optional argument to specify that B
!>      should be transposed.
    subroutine mult_2D(A, B, C, transpA, transpB)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp), intent(out) :: C(:, :)
        logical, intent(in), optional :: transpA, transpB
        logical :: transpA_, transpB_

        integer :: M, N, K_1, K_2, K

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

        M = size(A, merge(1, 2, .not. transpA_))
        call assert_(M == size(C, 1), 'Shape mismatch.')
        N = size(B, merge(2, 1, .not. transpB_))
        call assert_(N == size(C, 2), 'Shape mismatch.')
        K_1 = size(A, merge(2, 1, .not. transpA_))
        K_2 = size(B, merge(1, 2, .not. transpB_))
        call assert_(K_1 == K_2, 'Shape mismatch.')
        K = K_1

        call assert_(wp == r8, 'Precision mismatch for DGEMM')
        call dgemm_(merge('T', 'N', transpA_), merge('T', 'N',transpB_), &
                    M, N, K, 1._wp, A, size(A, 1), B, size(B, 1), &
                    0._wp, C, size(C, 1))
    end subroutine

!>  @brief
!>    Wrapper around dgemm for matrix-vector multiplication.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @paramin[in] A
!>  @paramin[in] x
!>  @paramin[out] y The shape of the output array is size(x)
!>  @paramin[in] transpA, Optional argument to specify that A
!>      should be transposed.
    subroutine mult_2D_1D(A, x, y, transpA)
        real(wp), intent(in) :: A(:, :), x(:)
        real(wp), intent(out) :: y(:)
        logical, intent(in), optional :: transpA
        logical :: transpA_

        integer :: M, N, K

        if (present(transpA)) then
          transpA_ = transpA
        else
          transpA_ = .false.
        end if

        M = size(A, merge(1, 2, .not. transpA_))
        call assert_(M == size(y, 1), 'Shape mismatch.')
        N = 1
        K = size(A, merge(2, 1, .not. transpA_))
        call assert_(K == size(x, 1), 'Shape mismatch.')

        call assert_(wp == r8, 'Precision mismatch for DGEMM')
        call dgemm_(merge('T', 'N', transpA_), 'N', &
                    M, N, K, 1._wp, A, size(A, 1), x, size(x, 1), &
                    0._wp, y, size(y, 1))
    end subroutine


!>  @brief
!>    Wrapper around dgemm for matrix-matrix multiplication with raw memory.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is important to note, that it is expected to pass in the
!>  actual vectors and not just the first element.
!>  So if `A_ptr` is the pointer to A in the work array
!>  and `A_L` is the length of A it is not possible to call
!>  this procedure just with `Work(A_ptr)`. It is necessary to call it
!>  with `Work(A_ptr : A_ptr + A_L - 1)`.
!>  (Which implies that the size of the arrays is automatically passed in.)
!>
!>  @paramin[in] A
!>  @paramin[in] rows_A Number of rows of A. The column number is deduced.
!>  @paramin[in] B
!>  @paramin[in] rows_B Number of rows of B. The column number is deduced.
!>  @paramin[in] C
!>  @paramin[out] C The shape of the output array is usually
!>      (rows_A * (size(B) / rows_B)) which changes of course, if
!>      A or B are transposed.
!>  @paramin[in] transpA, Optional argument to specify that A
!>      should be transposed.
!>  @paramin[in] transpB, Optional argument to specify that B
!>      should be transposed.
    subroutine mult_2D_raw(A, rows_A, B, rows_B, C, transpA, transpB)
        real(wp), intent(in), target :: A(:)
        integer, intent(in) :: rows_A
        real(wp), intent(in), target :: B(:)
        integer, intent(in) :: rows_B
        real(wp), intent(out), target :: C(:)
        logical, intent(in), optional :: transpA, transpB
        logical :: transpA_, transpB_

        integer :: shapeA(2), shapeB(2), shapeC(2)
        integer :: M, N, K
        integer :: K_1, K_2

        real(wp), pointer :: ptr_A(:, :), ptr_B(:, :), ptr_C(:, :)

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

        call assert_(modulo(size(A), rows_A) == 0, 'the number of rows has to be a divisor of size')
        shapeA = [rows_A, size(A) / rows_A]
        call assert_(modulo(size(B), rows_B) == 0, 'the number of rows has to be a divisor of size')
        shapeB = [rows_B, size(B) / rows_B]

        M = merge(shapeA(1), shapeA(2), .not. transpA_)
        N = merge(shapeB(2), shapeB(1), .not. transpB_)
        shapeC = [M, N]

        ptr_A(1 : shapeA(1), 1 : shapeA(2)) => A(:)
        ptr_B(1 : shapeB(1), 1 : shapeB(2)) => B(:)
        ptr_C(1 : shapeC(1), 1 : shapeC(2)) => C(:)

        call mult(ptr_A, ptr_B, ptr_C, transpA_, transpB_)
    end subroutine

!>  @brief
!>  Diagonalize A
!>
!>  @details
!>  Diagonalize A to fullfill `A V(:, j) = lambda(j) V(:, j)`.
!>  Wrapper around BLAS DSYEV.
!>
!>  @param[in] A 2D matrix to be diagonalized.
!>  @param[out] V 2D matrix which contains the Eigenvectors.
!>      The j-th column corresponds to the j-th Eigenvalue.
!>  @param[out] lambda 1D vector of Eigenvalues.
!>  @param[in] info Error code. If calling code does not check the error
!>      code, then this routine crashes upon errors.
    subroutine diagonalize(A, V, lambda, info)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: V(:, :), lambda(:)
        integer, optional, intent(out) :: info

        integer, parameter :: do_worksize_query = -1
        integer :: info_
        real(wp), allocatable :: work(:)
        real(wp) :: dummy(2), query_result(2)

        call assert_(symmetric(A), 'Matrix not symmetric')

        V(:, :) = A(:, :)
        call dsyev_('V', 'L', size(V, 2), dummy, size(V, 1), dummy, &
                    query_result, do_worksize_query, info_)

        if (info_ /= 0) then
            if (present(info)) then
                info = info_
            else
                call abort_('Error in diagonalize')
            end if
        end if

        call mma_allocate(work, int(query_result(1)))
        call dsyev_('V', 'L', size(V, 2), V, size(V, 1), lambda, &
                    work, size(work), info_)
        call mma_deallocate(work)

        if (info_ /= 0) then
            if (present(info)) then
                info = info_
            else
                call abort_('Error in diagonalize')
            end if
        end if
    end subroutine diagonalize

!>  @brief
!>  Order Eigenvectors by ascending eigenvalues.
    subroutine order_eigenvectors(V, lambda)
        real(wp), intent(inout) :: V(:, :), lambda(:)
        integer :: i
        integer, allocatable :: idx(:)
        call mma_allocate(idx, size(lambda))
        do i = 1, size(idx)
            idx(i) = i
        end do
        call sort(idx, leq)
        V(:, :) = V(:, idx)
        lambda(:) = lambda(idx)
        call mma_deallocate(idx)
    contains
        logical pure function leq(i, j)
            integer, intent(in) :: i, j
            leq = lambda(i) <= lambda(j)
        end function
    end subroutine

!> @brief
!> Determine the Eigenspaces.
!>
!>  @details
    subroutine determine_eigenspaces(lambda, dimensions)
        real(wp), intent(inout) :: lambda(:)
        integer, allocatable, intent(out) :: dimensions(:)

        integer, allocatable :: d_buffer(:)

        integer :: i, low
        integer :: n_spaces

        call assert_(all(lambda(2:) >= lambda(:size(lambda) - 1)), 'Not sorted')

        call mma_allocate(d_buffer, size(lambda))
        d_buffer(:) = 0

        low = 1
        n_spaces = 1
        do i = 1, size(lambda)
            d_buffer(n_spaces) = d_buffer(n_spaces) + 1
            if (i + 1 <= size(lambda)) then
                if (.not. (lambda(i) .isclose. lambda(i + 1))) then
                    n_spaces = n_spaces + 1
                    lambda(low : i) = mean(lambda(low : i))
                    low = i + 1
                end if
            else
                lambda(low : i) = mean(lambda(low : i))
            end if
        end do

        call mma_allocate(dimensions, n_spaces)
        dimensions(:) = d_buffer(: n_spaces)
        call mma_deallocate(d_buffer)

        contains
            real(wp) pure function mean(X)
                real(wp), intent(in) :: X(:)
                mean = sum(X) / size(X)
            end function
    end subroutine


    subroutine canonicalize_default(V, lambda)
        real(wp), intent(inout) :: V(:, :), lambda(:)
        character(*), parameter :: this_routine = 'canonicalize_default'


        integer :: info_
        integer :: low, i, j, d
        logical :: converged
        integer, allocatable :: idx(:), dimensions(:)
        real(wp), allocatable :: projections(:, :)

        call assert_(size(V, 1) == size(V, 2), 'non square matrix')
        call mma_allocate(projections, size(V, 1), size(V, 2))
        call mma_allocate(idx, size(V, 1))

        call order_eigenvectors(V, lambda)
        do j = 1, size(V, 2)
            V(:, j) = V(:, j) / norm(V(:, j))
        end do

        call determine_eigenspaces(lambda, dimensions)

        low = 1
        do j = 1, size(dimensions)
            converged = .false.
            do while (.not. converged)
                d = dimensions(j)
                call get_projections(V(:, low : low + d - 1), projections)

                do i = 1, size(idx)
                    idx(i) = i
                end do
                call sort(idx, geq)
                call sort(idx( : d), integer_leq)

                ! Get the first d projections with the largest overlap.
                block
                    real(wp) :: p(size(projections, 1))
                    converged = .true.
                    do i = 1, d
                        p = projections(:, idx(i)) / norm(projections(:, idx(i)))
                        converged = converged .and. all(V(:, low + i - 1) .isclose. p)
                        V(:, low + i - 1) = p
                    end do
                end block
            end do
            low = low + d
        end do

        call mma_deallocate(idx)
        call mma_deallocate(projections)
        call mma_deallocate(dimensions)
    contains

        logical pure function geq(i, j)
            integer, intent(in) :: i, j
            geq = norm(projections(:, i)) >= norm(projections(:, j))
        end function
    end subroutine


    ! Calculate the projections p_j of the e_j canonical unit vectors
    ! into the subspace spanned by M.
    subroutine get_projections(M, P)
        real(wp), intent(in) :: M(:, :)
        real(wp), intent(out) :: P(:, :)
        integer :: i, j

        P(:, :) = 0._wp
        ! calculate projections of basis b_j into eigenvectors v_i:
        ! p_j = sum_i < b_j | v_i > v_i
        ! if b_j are canonical unit vectors, the dot product is the j-th
        ! component of v_i
        do j = 1, size(P, 2)
            do i = 1, size(M, 2)
                P(:, j) = P(:, j) + M(j, i) * M(:, i)
            end do
        end do
    end subroutine

    elemental function isclose(a, b, atol, rtol) result(res)
        real(wp), intent(in) :: a, b
        real(wp), intent(in), optional :: atol, rtol
        logical :: res
        res = abs(a - b) <= (atol + rtol * abs(b))
    end function

! Operator functions may only have two arguments.
    elemental function isclose_for_operator(a, b) result(res)
        real(wp), intent(in) :: a, b
        logical :: res

        ! TODO(@Oskar, @Ignacio): decide on good constants
        res = isclose(a, b, atol=epsilon(a) * 10._wp, rtol=1e-8_wp)
    end function


!>  @brief
!>    Calculates v1^T S v2.
!>
!>  @author
!>    Oskar Weser
!>
!>  @details
!>  Calculates \f[ v1^T S v2 \f]
!>  S has to be a symmetric positive definite matrix, which is not
!>  tested.
!>  If S is ommited, it defaults to the unit matrix,
!>  i.e. the Euclidean dot-product.
      pure function dot_product_(v1, v2, S) result(dot)
        real(wp), intent(in) :: v1(:), v2(:)
        real(wp), intent(in), optional :: S(:, :)
        real(wp) :: dot

        if (present(S)) then
          dot = dot_product(matmul(S, v1), v2)
        else
          dot = dot_product(v1, v2)
        end if
      end function

!>  @brief
!>    The induced norm from the dot product given by S.
!>
!>  @details
!>    It is assumed, that S is symmetric and positive definite.
!>
!>  @author
!>    Oskar Weser
    pure function norm(v, S) result(L)
        real(wp), intent(in) :: v(:)
        real(wp), intent(in), optional :: S(:, :)
        real(wp) :: L
        if (present(S)) then
          L = sqrt(dot_product_(v, v, S))
        else
! One could use norm2 here, but Sun and PGI compilers don't know this.
          L = sqrt(sum(v**2))
        end if
    end function

    logical pure function symmetric(M)
        real(wp), intent(in) :: M(:, :)
        integer :: i, j

        symmetric = .true.
        do j = 1, size(M, 2)
            do i = j, size(M, 1)
                if (.not. (M(i, j) .isclose. M(j, i))) then
                    symmetric = .false.
                    return
                end if
            end do
        end do
    end function

    subroutine Lowdin(basis, ONB, S)
        real(wp), intent(in) :: basis(:, :)
        real(wp), intent(out) :: ONB(:, :)
        real(wp), intent(in), optional :: S(:, :)

        integer :: i
        real(wp), allocatable :: U(:, :), s_diag(:), X(:, :), &
              S_transf(:, :), tmp(:, :)

        call mma_allocate(S_transf, size(S, 1), size(S, 2))
        call mma_allocate(U, size(S, 1), size(S, 2))
        call mma_allocate(X, size(S, 1), size(S, 2))
        call mma_allocate(tmp, size(S, 1), size(S, 2))
        call mma_allocate(s_diag, size(S, 2))

! Transform AO-overlap matrix S to the overlap matrix of basis.
! S_transf = basis^T S basis
! We search X that diagonalizes S_transf
        if (present(S)) then
            call mult(S, basis, tmp)
            call mult(basis, tmp, S_transf, transpA=.true.)
        else
            call mult(basis, basis, S_transf, transpA=.true.)
        end if

        call diagonalize(S_transf, U, s_diag)
        call canonicalize(U, s_diag)

        call assert_(all(s_diag > 1.0d-10), &
           "Linear dependency detected. "// &
           "Lowdin can't cure it. Please use other ORTH keyword "// &
           "from {Gram_Schmidt, Canonical}.")

! X = U s_diag^{-1/2} U^T
        do i = 1, size(tmp, 2)
          tmp(:, i) = U(:, i) / sqrt(s_diag(i))
        end do
        call mult(tmp, U, X, transpB=.true.)
! With this X the overlap matrix S_transf has diagonal form.
! X^T basis^T S basis X = 1
! We finally have to convert to get the form:
! ONB^T S ONB = 1
        call mult(basis, X, ONB)

        call mma_deallocate(tmp)
        call mma_deallocate(s_diag)
        call mma_deallocate(X)
        call mma_deallocate(U)
        call mma_deallocate(S_transf)
    end subroutine Lowdin

    subroutine Canonical(basis, n_to_ON, ONB, n_new, S)
        real(wp), intent(in) :: basis(:, :)
        integer, intent(in) :: n_to_ON
        real(wp), intent(out) :: ONB(:, :)
        integer, intent(out) :: n_new
        real(wp), intent(in), optional :: S(:, :)

        logical :: lin_dep_detected
        integer :: i
        integer, allocatable :: idx(:)
        real(wp), allocatable :: &
            U(:, :), s_diag(:), S_transf(:, :), X(:, :), tmp(:, :)

        call mma_allocate(S_transf, size(S, 1), size(S, 2))
        call mma_allocate(U, size(S, 1), size(S, 2))
        call mma_allocate(s_diag, size(S, 2))
        call mma_allocate(X, size(S, 1), size(S, 2))
        call mma_allocate(idx, size(S, 1))
        call mma_allocate(tmp, size(S, 1), size(S, 2))

! Transform AO-overlap matrix S to the overlap matrix of basis.
! We search X that diagonalizes basis^T S basis
        if (present(S)) then
            call mult(S, basis, tmp)
            call mult(basis, tmp, S_transf, transpA=.true.)
        else
            call mult(basis, basis, S_transf, transpA=.true.)
        end if

        call diagonalize(S_transf, U, s_diag)
        call canonicalize(U, s_diag)

        idx(:) = argsort(s_diag, ge_r)
        U(:, :) = U(:, idx)
        s_diag(:) = s_diag(idx)

        i = 0
        lin_dep_detected = .false.
        do while(.not. lin_dep_detected .and. i < n_to_ON)
          if (s_diag(i + 1) < 1.0d-10) then
            n_new = i
            lin_dep_detected = .true.
          end if
          i = i + 1
        end do
        if (.not. lin_dep_detected) n_new = n_to_ON

! X = U s_diag^{-1/2}
        do i = 1, n_new
          X(:, i) = U(:, i) / sqrt(s_diag(i))
        end do
! With this X the overlap matrix S_transf has diagonal form.
! X^T basis^T S basis X = 1
! We finally have to convert to get the form:
! ONB^T S ONB = 1
        ONB(:, n_new + 1:) = basis(:, n_new + 1 :)
        call mult(basis, X(:, :n_new), ONB(:, :n_new))

        call mma_deallocate(tmp)
        call mma_deallocate(X)
        call mma_deallocate(idx)
        call mma_deallocate(s_diag)
        call mma_deallocate(U)
        call mma_deallocate(S_transf)
    end subroutine Canonical

    subroutine Gram_Schmidt(basis, n_to_ON, ONB, n_new, S)
        real(wp), intent(in) :: basis(:, :)
        integer, intent(in) :: n_to_ON
        real(wp), target, intent(out) :: ONB(:, :)
        integer, intent(out) :: n_new
        real(wp), intent(in), optional :: S(:, :)

        real(wp) :: L
        integer :: i, j
        logical :: lin_dep_detected, improve_solution

        real(wp), allocatable :: previous(:), correction(:), v(:)
        real(wp), pointer :: curr(:)

        call mma_allocate(previous, size(basis, 1))
        call mma_allocate(correction, size(basis, 1))
        call mma_allocate(v, size(basis, 1))

        n_new = 0
        ONB(:, n_to_ON + 1 :) = basis(:, n_to_ON + 1 :)
        do i = 1, n_to_ON
          curr => ONB(:, n_new + 1)
          curr = basis(:, i)

          improve_solution = .true.
          lin_dep_detected = .false.
          do while (improve_solution .and. .not. lin_dep_detected)
            correction = 0._wp
            if (present(S)) then
                call mult(S, curr, v)
            else
                v(:) = curr(:)
            end if

            do j = 1, n_new
              correction(:) = &
                    correction(:) + ONB(:, j) * dot_product(ONB(:, j), v)
            end do
            curr = curr - correction
            improve_solution = norm(correction, S=S) > 0.1_wp
            L = norm(curr, S=S)
            lin_dep_detected = L < 1.0e-10_wp
            if (.not. lin_dep_detected) then
              curr = curr / L
            end if
            if (.not. (improve_solution .or. lin_dep_detected)) then
              n_new = n_new + 1
            end if
          end do
        end do
        ONB(:, n_new + 1 : n_to_ON) = basis(:, n_new + 1 : n_to_ON)

        call mma_deallocate(v)
        call mma_deallocate(correction)
        call mma_deallocate(previous)
    end subroutine Gram_Schmidt



    subroutine abort_(message)
        character(*), intent(in) :: message
        call WarningMessage(2, message)
        call QTrace()
        call Abend()
    end subroutine

    !> @brief
    !>    Runtime check, that is switched off in Debug mode
    subroutine assert_(test_expression, message)
        logical, intent(in) :: test_expression
        character(*), intent(in) :: message
#ifdef _DEBUG_
        if (.not. test_expression) then
            call abort_(message)
        end if
#endif
    end subroutine

    !> @brief
    !>    Runtime check, that is not switched off in Debug mode
    subroutine assume_(test_expression, message)
        logical, intent(in) :: test_expression
        character(*), intent(in) :: message
        if (.not. test_expression) then
            call abort_(message)
        end if
    end subroutine

end module
