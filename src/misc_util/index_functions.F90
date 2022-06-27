!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module Index_Functions

use Definitions, only: iwp

implicit none
private

public :: C_Ind, C_Ind3, nTri_Elem, nTri_Elem1, nTri3_Elem, nTri3_Elem1

#include "macros.fh"

contains

! Index for a total order of l=lx+ly+lz, sorted with priority lx > ly > lz
! e.g., for l=2 this gives:
!
!   lx ly lz  C_Ind
!   ---------------
!    2  0  0    1
!    1  1  0    2
!    1  0  1    3
!    0  2  0    4
!    0  1  1    5
!    0  0  2    6
pure function C_Ind(l,lx,lz)
  integer(kind=iwp) :: C_Ind
  integer(kind=iwp), intent(in) :: l, lx, lz
  C_Ind = (l-lx)*(l-lx+1)/2+lz+1
end function C_Ind

! Same as C_Ind, but from lx, ly, lz directly
pure function C_Ind3(lx,ly,lz)
  integer(kind=iwp) :: C_Ind3
  integer(kind=iwp), intent(in) :: lx, ly, lz
  unused_var(lx)
  C_Ind3 = (ly+lz)*(ly+lz+1)/2+lz+1
end function C_Ind3

! Number of elements in a triangular matrix of side n
!
!   n  nTri_Elem
!   ------------
!   0     0
!   1     1
!   2     3
!   3     6
!   4    10
pure function nTri_Elem(n)
  integer(kind=iwp) :: nTri_Elem
  integer(kind=iwp), intent(in) :: n
  nTri_Elem = n*(n+1)/2
end function nTri_Elem

! Same as nTri_Elem, but with n=l+1
pure function nTri_Elem1(l)
  integer(kind=iwp) :: nTri_Elem1
  integer(kind=iwp), intent(in) :: l
  nTri_Elem1 = (l+1)*(l+2)/2
end function nTri_Elem1

! Like nTri_Elem, but for 3-dimensional arrays
!
!   n  nTri3_Elem
!   -------------
!   0     0
!   1     1
!   2     4
!   3    10
!   4    20
pure function nTri3_Elem(n)
  integer(kind=iwp) :: nTri3_Elem
  integer(kind=iwp), intent(in) :: n
  nTri3_Elem = n*(n+1)*(n+2)/6
end function nTri3_Elem

! Same as nTri3_Elem, but with n=l+1
pure function nTri3_Elem1(l)
  integer(kind=iwp) :: nTri3_Elem1
  integer(kind=iwp), intent(in) :: l
  nTri3_Elem1 = (l+1)*(l+2)*(l+3)/6
end function nTri3_Elem1

end module
