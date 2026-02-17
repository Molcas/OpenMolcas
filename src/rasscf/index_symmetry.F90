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

use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Two, Half
use Definitions, only: iwp

implicit none
private

public :: one_el_idx, one_el_idx_flatten, two_el_idx, two_el_idx_flatten

interface one_el_idx
  module procedure :: array_one_el_idx, tuple_one_el_idx
end interface

interface two_el_idx
  module procedure :: array_two_el_idx, tuple_two_el_idx
end interface

interface one_el_idx_flatten
  module procedure :: array_1el_idx_flatten, tuple_1el_idx_flatten
end interface

interface two_el_idx_flatten
  module procedure :: array_2el_idx_flatten, tuple_2el_idx_flatten, tuple_2el_idx_flatten_2
end interface

contains

pure subroutine tuple_one_el_idx(n,i,j)

  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp), intent(out) :: i, j

  i = ceiling(-Half+sqrt(Two*n))
  j = n-nTri_Elem(i-1)

end subroutine tuple_one_el_idx

pure subroutine array_one_el_idx(n,idx)

  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp), intent(out) :: idx(2)

  idx(1) = ceiling(-Half+sqrt(Two*n))
  idx(2) = n-nTri_Elem(idx(1)-1)

end subroutine array_one_el_idx

pure function array_1el_idx_flatten(idx) result(n)

  integer(kind=iwp), intent(in) :: idx(2)
  integer(kind=iwp) :: n

  n = tuple_1el_idx_flatten(idx(1),idx(2))

end function array_1el_idx_flatten

pure function tuple_1el_idx_flatten(i,j) result(n)

  integer(kind=iwp) :: n
  integer(kind=iwp), intent(in) :: i, j

  n = iTri(i,j)

end function tuple_1el_idx_flatten

pure subroutine tuple_two_el_idx(n,iorb,jorb,korb,lorb)

  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp), intent(out) :: iorb, jorb, korb, lorb
  integer(kind=iwp) :: ijidx, klidx

  ijidx = ceiling(-Half+sqrt(Two*n))
  klidx = n-nTri_Elem(ijidx-1)

  iorb = ceiling(-Half+sqrt(Two*ijidx))
  jorb = ijidx-nTri_Elem(iorb-1)
  korb = ceiling(-Half+sqrt(Two*klidx))
  lorb = klidx-nTri_Elem(korb-1)

end subroutine tuple_two_el_idx

pure subroutine array_two_el_idx(n,idx)

  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp), intent(out) :: idx(4)
  integer(kind=iwp) :: ijidx, klidx

  ijidx = ceiling(-Half+sqrt(Two*n))
  klidx = n-nTri_Elem(ijidx-1)

  idx(1) = ceiling(-Half+sqrt(Two*ijidx))
  idx(2) = ijidx-nTri_Elem(idx(1)-1)
  idx(3) = ceiling(-Half+sqrt(Two*klidx))
  idx(4) = klidx-nTri_Elem(idx(3)-1)

end subroutine array_two_el_idx

function array_2el_idx_flatten(idx) result(n)

  integer(kind=iwp) :: n
  integer(kind=iwp), intent(in) :: idx(4)

  n = tuple_2el_idx_flatten(idx(1),idx(2),idx(3),idx(4))

end function array_2el_idx_flatten

function tuple_2el_idx_flatten(p,q,r,s) result(pqrs)

  integer(kind=iwp) :: pqrs
  integer(kind=iwp), intent(in) :: p, q, r, s
  integer(kind=iwp) :: pq, rs

  pqrs = tuple_2el_idx_flatten_2(p,q,r,s,pq,rs)

end function tuple_2el_idx_flatten

function tuple_2el_idx_flatten_2(p,q,r,s,pq,rs) result(pqrs)

  integer(kind=iwp) :: pqrs
  integer(kind=iwp), intent(in) :: p, q, r, s
  integer(kind=iwp), intent(out) :: pq, rs

  pq = iTri(p,q)
  rs = iTri(r,s)
  pqrs = iTri(pq,rs)

end function tuple_2el_idx_flatten_2

end module index_symmetry
