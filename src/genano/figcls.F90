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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine prints the trailer part to the postscript figure of     *
! occupation numbers.                                                  *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine FigCls(lu)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lu

write(lu,'(a)') '%---'
write(lu,'(a)') 'grestore'
write(lu,'(a)') 'showpage'
write(lu,'(a)') '%.AutoFeed'

return

end subroutine FigCls
