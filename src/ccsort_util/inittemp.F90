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

subroutine inittemp(num)
! this routine initializes status matrix
! num - number of files to be used (I)

implicit real*8(a-h,o-z)
#include "reorg.fh"
integer num
! help variables
integer nhelp

do nhelp=1,num
  stattemp(nhelp) = 0
  nrectemp(nhelp) = 0
  lrectemp(nhelp) = 0
end do

return

end subroutine inittemp
