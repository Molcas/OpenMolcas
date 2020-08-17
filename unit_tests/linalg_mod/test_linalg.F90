module test_linalg_mod
    use fruit
    use definitions, only: wp
    use linalg_mod, only: mult, operator(.isclose.), Gram_Schmidt, symmetric, &
        diagonalize, canonicalize, norm
    implicit none
    private
    public :: test_mult, test_diagonalization
    real(wp), parameter :: tolerance = 10._wp**2 * epsilon(1._wp)

contains

    subroutine test_mult()
        block
            real(wp), target :: A(20, 10), B(10, 20), C(size(A, 1), size(B, 2))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(A, B)

            call mult(A, B, C)
            call assert_true(all(C .isclose. expected))

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, shape(A), raw_B, shape(B), raw_C)
            call assert_true(all(C .isclose. expected))
        end block

        block
            real(wp), target :: A(20, 10), B(20, 10), C(size(A, 2), size(B, 2))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(transpose(A), B)

            call mult(A, B, C, transpA=.true.)
            call assert_true(all(C .isclose. expected))

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, shape(A), raw_B, shape(B), raw_C, transpA=.true.)
            call assert_true(all(C .isclose. expected))
        end block

        block
            real(wp), target :: A(20, 10), B(20, 10), C(size(A, 1), size(B, 1))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(A, transpose(B))

            call mult(A, B, C, transpB=.true.)
            call assert_true(all(C .isclose. expected))

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, shape(A), raw_B, shape(B), raw_C, transpB=.true.)
            call assert_true(all(C .isclose. expected))
        end block

        block
            real(wp), target :: A(20, 10), B(10, 20), C(size(A, 2), size(B, 1))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(transpose(A), transpose(B))

            call mult(A, B, C, transpA=.true., transpB=.true.)
            call assert_true(all(C .isclose. expected))

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, shape(A), raw_B, shape(B), raw_C, transpA=.true., transpB=.true.)
            call assert_true(all(C .isclose. expected))
        end block
    end subroutine

    subroutine test_diagonalization()
        integer, parameter :: test_size = 10
        ! dimension of the Eigenspaces
        integer, parameter :: dimension_E(4) = [5, 5, 1, 5]
        real(wp) :: lambdas(sum(dimension_E))
        real(wp) :: M(size(lambdas), size(lambdas))
        real(wp) :: V(size(M, 1), size(M, 2))
        real(wp) :: U(size(M, 1), size(M, 2))
        real(wp) :: test_V(size(M, 1), size(M, 2))

        integer :: i, i_test
        integer :: offset

        create_test_matrix : block
            offset = 1
            do i = 1, size(dimension_E)
                lambdas(offset : offset + dimension_E(i) - 1) = real(i, kind=wp)
                offset = offset + dimension_E(i)
            end do

            M = 0._wp
            do i = 1, size(M, 1)
                M(i, i) = lambdas(i)
            end do

            call get_rand_orthogonal(U)
            M = matmul(matmul(transpose(U), M), U)
            call assert_true(symmetric(M))
        end block create_test_matrix

        call diagonalize(M, V, lambdas)

        test_canonical_basis : block

            call canonicalize(V, lambdas)

            do i_test = 1, test_size
                call create_test_V(V, dimension_E, test_V)
                call canonicalize(test_V, lambdas)
                call assert_true(all(test_V .isclose. V))
            end do
        end block test_canonical_basis

        ! Use Roland's contraints
        test_general_basis : block
            real(wp) :: ref(size(V, 1), size(V, 2))

            call assert_true(size(V, 1) == size(V, 2))

            ref(:, :) = 0._wp

            do i = 1, size(V, 2)
                ref(i, i) = 1._wp
                ref(mod(i + 1, size(V, 1)), i) = -1._wp
            end do

            call canonicalize(V, lambdas, ref)

            do i_test = 1, test_size
                call create_test_V(V, dimension_E, test_V)
                call canonicalize(test_V, lambdas, ref)
                call assert_true(all(test_V .isclose. V))
            end do
        end block test_general_basis

    end subroutine

!>  @brief
!>    Create random rotations inside degenerate Eigenspaces
!>
!>  @author
!>    Oskar Weser
!>
!>  @param[in] V The eigenvectors. It is assumed, that Eigenvectors from the same Eigenspace
!>          are neighbouring.
!>  @param[in] dimension_E The dimension for each Eigenspace.
!>  @param[out] test_V The eigenvectors where degenerate Eigenspaces are randomly changed.
    subroutine create_test_V(V, dimension_E, test_V)
        real(wp), intent(in) :: V(:, :)
        integer, intent(in) :: dimension_E(:)
        real(wp), intent(out) :: test_V(:, :)

        integer :: offset, i

        offset = 1
        test_V = 0._wp
        do i = 1, size(dimension_E)
        block
            real(wp) :: U(dimension_E(i), dimension_E(i))
            integer :: j, k
            call get_rand_orthogonal(U)

            do j = 1, dimension_E(i)
                do k = 1, dimension_E(i)
                    test_V(:, offset + j - 1) = test_V(:, offset + j - 1) + U(k, j) * V(:, offset + k - 1)
                end do
            end do
            offset = offset + dimension_E(i)
        end block
        end do
    end subroutine

    subroutine get_rand_orthogonal(U)
        real(wp), intent(out) :: U(:, :)
        integer :: n_new
        real(wp) :: basis(size(U, 1), size(U, 2))
        logical :: linear_independent

        linear_independent = .false.
        do while (.not. linear_independent)
            call random_number(basis)
            call Gram_Schmidt(basis, size(basis, 2), U, n_new)
            linear_independent = n_new == size(basis, 2)
        end do
    end subroutine

end module test_linalg_mod

program test_linalg

    use fruit
    use test_linalg_mod, only: test_mult, test_diagonalization

    implicit none
    integer :: failed_count
    integer, parameter :: seed(20) = &
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

    block

        call init_fruit()
        call random_seed(put=seed)
        call inimem()

        call test_linalg_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0)  error stop
    end block

contains

    subroutine test_linalg_driver()
        call run_test_case(test_mult, "test_mult")
        call run_test_case(test_diagonalization, "test_diagonalization")
    end subroutine
end program test_linalg
