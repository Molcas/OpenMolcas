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

integer function IELSUM(IVEC,NELMNT)
! Sum elements of integer vector IVEC

implicit none

integer IVEC(*), NELMNT
integer ISUM, IELMNT

ISUM = 0
do IELMNT=1,NELMNT
  ISUM = ISUM+IVEC(IELMNT)
end do

IELSUM = ISUM

end function IELSUM
