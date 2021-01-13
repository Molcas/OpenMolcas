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
! Copyright (C) 2019, Oskar Weser                                      *
!***********************************************************************
      module index_symmetry
        implicit none
        private
        public :: one_el_idx, one_el_idx_flatten, two_el_idx,
     &      two_el_idx_flatten
        save

        interface one_el_idx
          module procedure array_one_el_idx, tuple_one_el_idx
        end interface

        interface two_el_idx
          module procedure array_two_el_idx, tuple_two_el_idx
        end interface

        interface one_el_idx_flatten
          module procedure array_1el_idx_flatten, tuple_1el_idx_flatten
        end interface

        interface two_el_idx_flatten
          module procedure array_2el_idx_flatten,
     &        tuple_2el_idx_flatten, tuple_2el_idx_flatten_2
        end interface

      contains

        pure subroutine tuple_one_el_idx(n, i, j)
          integer, intent(in) :: n
          integer, intent(out) :: i, j

          i = ceiling(-0.5d0 + sqrt(2.0d0 * n))
          j = n - (i - 1) * i / 2
        end subroutine

        pure subroutine array_one_el_idx(n, idx)
          integer, intent(in) :: n
          integer, intent(out) :: idx(2)

          idx(1) = ceiling(-0.5d0 + sqrt(2.0d0 * n))
          idx(2) = n - (idx(1) - 1) * idx(1) / 2
        end subroutine

        pure function array_1el_idx_flatten(idx) result(n)
          integer, intent(in) :: idx(2)
          integer :: n
          n = tuple_1el_idx_flatten(idx(1), idx(2))
        end function

        pure function tuple_1el_idx_flatten(i, j) result(n)
          integer, intent(in) :: i, j
          integer :: n
          integer :: p, q
          p = max(i, j)
          q = min(i, j)
          n = q + p * (p - 1) / 2
        end function

        pure subroutine tuple_two_el_idx(n, iorb, jorb, korb, lorb)
          integer, intent(in) :: n
          integer, intent(out) :: iorb, jorb, korb, lorb

          integer :: ijidx, klidx

          ijidx = ceiling(-0.5d0 + sqrt(2.0d0 * n))
          klidx = n - (ijidx - 1) * ijidx / 2

          iorb = ceiling(-0.5d0 + sqrt(2.0d0 * ijidx))
          jorb = ijidx - (iorb - 1) * iorb / 2
          korb = ceiling(-0.5d0 + sqrt(2.0d0 * klidx))
          lorb = klidx - (korb - 1) * korb / 2
        end subroutine

        pure subroutine array_two_el_idx(n, idx)
          integer, intent(in) :: n
          integer, intent(out) :: idx(4)

          integer :: ijidx, klidx

          ijidx = ceiling(-0.5d0 + sqrt(2.0d0 * n))
          klidx = n - (ijidx - 1) * ijidx / 2

          idx(1) = ceiling(-0.5d0 + sqrt(2.0d0 * ijidx))
          idx(2) = ijidx - (idx(1) - 1) * idx(1) / 2
          idx(3) = ceiling(-0.5d0 + sqrt(2.0d0 * klidx))
          idx(4) = klidx - (idx(3) - 1) * idx(3) / 2
        end subroutine

        function array_2el_idx_flatten(idx) result(n)
          integer, intent(in) :: idx(4)
          integer :: n
          n = tuple_2el_idx_flatten(idx(1), idx(2), idx(3), idx(4))
        end function

        function tuple_2el_idx_flatten(p, q, r, s) result(pqrs)
          integer, intent(in) :: p, q, r, s
          integer :: pqrs
          integer :: pq, rs
          pqrs = tuple_2el_idx_flatten_2(p, q, r, s, pq, rs)
        end function

        function tuple_2el_idx_flatten_2(p, q, r, s, pq, rs)result(pqrs)
          integer, intent(in) :: p, q, r, s
          integer, intent(out) :: pq, rs
          integer :: pqrs
          if (p >= q) pq = p * (p - 1) / 2 + q
          if (p < q) pq = q * (q - 1) / 2 + p
          if (r >= s) rs = r * (r - 1) / 2 + s
          if (r < s) rs = s * (s - 1) / 2 + r
          if (pq >= rs) pqrs = pq * (pq - 1) / 2 + rs
          if (pq < rs) pqrs = rs * (rs - 1) / 2 + pq
        end function
      end module index_symmetry
