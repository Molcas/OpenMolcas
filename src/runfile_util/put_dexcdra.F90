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

subroutine Put_dExcdRa(dExcdRa,ndExcdRa)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ndExcdRa
real(kind=wp) :: dExcdRa(ndExcdRa)
character(len=24) :: Label

Label = 'dExcdRa'
call Put_dArray(Label,dExcdRa,ndExcdRa)

return

end subroutine Put_dExcdRa
