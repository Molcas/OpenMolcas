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
          module procedure array_1el_idx, tuple_1el_idx
        end interface

        interface one_el_idx_flatten
          module procedure array_1el_idx_flatten, tuple_1el_idx_flatten
        end interface

        interface two_el_idx
          module procedure array_2el_idx, tuple_2el_idx
        end interface

        interface two_el_idx_flatten
          module procedure array_2el_idx_flatten,
     &        tuple_2el_idx_flatten, tuple_2el_idx_flatten_2
        end interface

      contains

        pure function array_1el_idx(n) result(idx)
          integer, intent(in) :: n
          integer :: idx(2)
          idx(1) = ceiling(-0.5d0 + sqrt(2.0d0 * n))
          idx(2) = n - (idx(1) - 1) * idx(1) / 2
        end function

        function tuple_1el_idx(n, i, j) result(idx)
          integer, intent(in) :: n
          integer, intent(out) :: i, j
          integer :: idx(2)
          idx = array_1el_idx(n)
          i = idx(1)
          j = idx(2)
        end function

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

        function tuple_2el_idx(n, iorb, jorb, korb, lorb) result(idx)
          integer, intent(in) :: n
          integer, intent(out) :: iorb, jorb, korb, lorb
          integer :: idx(4)

          integer :: ijidx, klidx

          ijidx = ceiling(-0.5d0 + sqrt(2.0d0 * n))
          klidx = n - (ijidx - 1) * ijidx / 2

          iorb = ceiling(-0.5d0 + sqrt(2.0d0 * ijidx))
          jorb = ijidx - (iorb - 1) * iorb / 2
          korb = ceiling(-0.5d0 + sqrt(2.0d0 * klidx))
          lorb = klidx - (korb - 1) * korb / 2
          idx = [iorb, jorb, korb, lorb]
        end function

        function array_2el_idx(n) result(idx)
          integer, intent(in) :: n
          integer :: idx(4)
          integer :: iorb, jorb, korb, lorb
          idx = tuple_2el_idx(n, iorb, jorb, korb, lorb )
        end function

        pure function array_2el_idx_flatten(idx) result(n)
          integer, intent(in) :: idx(4)
          integer :: n
          n = tuple_2el_idx_flatten(idx(1), idx(2), idx(3), idx(4))
        end function

        pure function tuple_2el_idx_flatten(p, q, r, s) result(pqrs)
          integer, intent(in) :: p, q, r, s
          integer :: pqrs
          integer :: pq, rs
          if (p >= q) pq = p * (p - 1) / 2 + q
          if (p < q) pq = q * (q - 1) / 2 + p
          if (r >= s) rs = r * (r - 1) / 2 + s
          if (r < s) rs = s * (s - 1) / 2 + r
          if (pq >= rs) pqrs = pq * (pq - 1) / 2 + rs
          if (pq < rs) pqrs = rs * (rs - 1) / 2 + pq
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
