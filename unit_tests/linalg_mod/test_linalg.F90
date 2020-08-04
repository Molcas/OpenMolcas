module test_linalg_mod
    use fruit
    use definitions, only: wp
    use linalg_mod, only: mult
    implicit none
    private
    public :: test_mult
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
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, size(A, 1), raw_B, size(B, 1), raw_C)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)
        end block

        block
            real(wp), target :: A(20, 10), B(20, 10), C(size(A, 2), size(B, 2))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(transpose(A), B)

            call mult(A, B, C, transpA=.true.)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, size(A, 1), raw_B, size(B, 1), raw_C, transpA=.true.)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)
        end block

        block
            real(wp), target :: A(20, 10), B(20, 10), C(size(A, 1), size(B, 1))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(A, transpose(B))

            call mult(A, B, C, transpB=.true.)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, size(A, 1), raw_B, size(B, 1), raw_C, transpB=.true.)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)
        end block

        block
            real(wp), target :: A(20, 10), B(10, 20), C(size(A, 2), size(B, 1))
            real(wp) :: expected(size(C, 1), size(C, 2))
            real(wp), pointer :: raw_A(:), raw_B(:), raw_C(:)

            call random_number(A)
            call random_number(B)

            expected = matmul(transpose(A), transpose(B))

            call mult(A, B, C, transpA=.true., transpB=.true.)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)

            raw_A(1 : size(A)) => A(:, :)
            raw_B(1 : size(B)) => B(:, :)
            raw_C(1 : size(C)) => C(:, :)

            call mult(raw_A, size(A, 1), raw_B, size(B, 1), raw_C, transpA=.true., transpB=.true.)
            call assert_equals(C, expected, size(expected, 1), size(expected, 2), delta=tolerance)
        end block
    end subroutine

end module test_linalg_mod

program test_linalg

    use fruit
    use test_linalg_mod, only: test_mult

    implicit none
    integer :: failed_count, err

    integer :: n
    block

        call init_fruit()

        call test_linalg_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0)  error stop
    end block

contains

    subroutine test_linalg_driver()
        call run_test_case(test_mult, "test_mult")
    end subroutine
end program test_linalg
