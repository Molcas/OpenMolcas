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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine DIAG_CPF(ICASE,JSY,HDIAG,FC,FIJ,FJI)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), JSY(*)
real(kind=wp), intent(_OUT_) :: HDIAG(*)
real(kind=wp), intent(in) :: FC(*), FIJ(*), FJI(*)

call IIJJ_CPF(ICASE,JSY,HDIAG,FC,FIJ,FJI)
call IJIJ_CPF(JSY,HDIAG,FJI)

return

end subroutine DIAG_CPF
