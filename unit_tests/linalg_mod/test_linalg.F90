module test_linalg_mod
    use fruit
    use definitions, only: wp
    implicit none
    private



contains

end module test_gasci_mod

program test_gasci_program

    use mpi
    use fruit

    implicit none
    integer :: failed_count, err

    integer :: n
    block

        call init_fruit()

        call test_linalg_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_linalg_program', 'failed_tests')
    end block

contains

    subroutine test_linalg_driver()
    end subroutine
end program test_linalg_mod
