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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine bufio_end_cvb()

use casvb_global, only: file_id, nbuf
use Definitions, only: wp

implicit none
real(kind=wp) :: dnbuf(1)

call bufio_wrbuf_cvb()
dnbuf(1) = real(nbuf,kind=wp)
call wrlow_cvb(dnbuf,1,file_id,0)

return

end subroutine bufio_end_cvb
