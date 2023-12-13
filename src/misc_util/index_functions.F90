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

public :: C3_Ind, C3_Ind3, C_Ind, C_Ind3, C_Ind3_Rev, iTri, iTri_Rev, nTri3_Elem, nTri3_Elem1, nTri_Elem, nTri_Elem1

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

! Cumulative C_Ind, including all previous l values
pure function C3_Ind(l,lx,lz)
  integer(kind=iwp) :: C3_Ind
  integer(kind=iwp), intent(in) :: l, lx, lz
  C3_Ind = l*(l+1)*(l+2)/6+(l-lx)*(l-lx+1)/2+lz+1
end function C3_Ind

! Same as C3_Ind, but from lx, ly, lz directly
pure function C3_Ind3(lx,ly,lz)
  integer(kind=iwp) :: C3_Ind3
  integer(kind=iwp), intent(in) :: lx, ly, lz
  C3_Ind3 = (lx+ly+lz)*(lx+ly+lz+1)*(lx+ly+lz+2)/6+(ly+lz)*(ly+lz+1)/2+lz+1
end function C3_Ind3

! Inverse of C_Ind3: from the index and l, return lx, ly, lz
pure function C_Ind3_Rev(lxyz,l)
  use Constants, only: Seven, Eight
  integer(kind=iwp) :: C_Ind3_Rev(3)
  integer(kind=iwp), intent(in) :: lxyz, l
  integer(kind=iwp) :: lx, ly, lyz, lz
  lyz = (int(sqrt(Eight*lxyz-Seven))-1)/2
  lz = lxyz-lyz*(lyz+1)/2-1
  ly = lyz-lz
  lx = l-lyz
  C_Ind3_Rev(:) = [lx,ly,lz]
end function C_Ind3_Rev

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

! Index of element i,j in a triangularly stored matrix
pure function iTri(i,j)
  integer(kind=iwp) :: iTri
  integer(kind=iwp), intent(in) :: i, j
  if (j > i) then
    iTri = j*(j-1)/2+i
  else
    iTri = i*(i-1)/2+j
  end if
end function iTri

! Inverse of iTri: from the index, return i and j
pure function iTri_Rev(ij)
  use Constants, only: Seven, Eight
  integer(kind=iwp) :: iTri_Rev(2)
  integer(kind=iwp), intent(in) :: ij
  integer(kind=iwp) :: i, j
  i = (int(sqrt(Eight*ij-Seven))+1)/2
  j = ij-i*(i-1)/2
  iTri_Rev(:) = [i,j]
end function iTri_Rev

end module
