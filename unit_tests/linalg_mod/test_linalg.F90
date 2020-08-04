module test_linalg_mod
    use fruit
    use definitions, only: wp
    use linalg_mod, only: mult
    implicit none
    private
    public :: test_mult

contains

    subroutine test_mult()
        real(wp) :: A(2, 2), B(2, 2), C(2, 2)

        call random_number(A)
        call random_number(B)

        call mult(A, B, C)
        call assert_equals(C, matmul(A, B), 2, 2)

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
