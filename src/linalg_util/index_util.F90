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

module Index_util

use Definitions, only: iwp

implicit none
private
public :: iTri, iTri0, nTriElem, nTri0Elem

! iTri, nTriElem:   Index and number of elements in a triangular storage
! iTri0, nTri0Elem: idem for 0-based indices (e.g. angular momentum)

contains

pure function iTri(i,j)
  integer(kind=iwp) :: iTri
  integer(kind=iwp), intent(in) :: i,j
  iTri = max(i,j)*(max(i,j)-1)/2+min(i,j)
end function iTri

pure function iTri0(i,j)
  integer(kind=iwp) :: iTri0
  integer(kind=iwp), intent(in) :: i,j
  iTri0 = max(i,j)*(max(i,j)+1)/2+min(i,j)+1
end function iTri0

pure function nTriElem(n)
  integer(kind=iwp) :: nTriElem
  integer(kind=iwp), intent(in) :: n
  nTriElem = n*(n+1)/2
end function nTriElem

pure function nTri0Elem(n)
  integer(kind=iwp) :: nTri0Elem
  integer(kind=iwp), intent(in) :: n
  nTri0Elem = (n+1)*(n+2)/2
end function nTri0Elem

end module Index_util
