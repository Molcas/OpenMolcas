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

subroutine VanishT1(wrk,wrksize)
! this routine does:
! Vanish T1o

implicit none
#include "chcc1.fh"
#include "wrk.fh"
! help variables
integer len_

len_ = no*nv
call mv0zero(len_,len_,wrk(PosT1o))

return

end subroutine VanishT1
