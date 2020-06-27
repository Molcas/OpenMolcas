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
    use definitions, only: wp, r8
    implicit none
    private
    public :: mult


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
        real(wp), intent(in) :: A(:)
        integer, intent(in) :: rows_A
        real(wp), intent(in) :: B(:)
        integer, intent(in) :: rows_B
        real(wp), intent(out) :: C(:)
        logical, intent(in), optional :: transpA, transpB
        logical :: transpA_, transpB_

        integer :: shapeA(2), shapeB(2), shapeC(2)
        integer :: M, N, K
        integer :: K_1, K_2

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

        K_1 = merge(shapeA(2), shapeA(1), .not. transpA_)
        K_2 = merge(shapeB(1), shapeB(2), .not. transpB_)
        call assert_(K_1 == K_2, 'Shape mismatch.')
        K = K_1

        call assert_(wp == r8, 'Precision mismatch for DGEMM')
        call dgemm_(merge('T', 'N', transpA_), merge('T', 'N',transpB_), &
                    M, N, K, 1._wp, A, shapeA(1), B, shapeB(1), &
                    0._wp, C, shapeC(1))
    end subroutine

    subroutine abort_(message)
        character(*), intent(in) :: message
        call WarningMessage(2, message)
        call QTrace()
        call Abend()
    end subroutine

    subroutine assert_(test_expression, message)
        logical, intent(in) :: test_expression
        character(*), intent(in) :: message
        if (.not. test_expression) then
            call abort_(message)
        end if
    end subroutine

end module
