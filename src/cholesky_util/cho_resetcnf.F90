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

subroutine CHO_RESETCNF()
!
! Purpose: reset configuration of decomposition to that read from
!          restart file. Original configuration will be saved in
!          restart common block.

implicit none
#include "cholesky.fh"

call CHO_DSWAP(THRCOM,XTHRCOM)
call CHO_DSWAP(THRDIAG,XTHRDIAG)
call CHO_DSWAP(DAMP(1),XDAMP(1))
call CHO_DSWAP(DAMP(2),XDAMP(2))
call CHO_DSWAP(SPAN,XSPAN)
call CHO_DSWAP(THRNEG,XTHRNEG)
call CHO_DSWAP(WARNEG,XWARNEG)
call CHO_DSWAP(TOONEG,XTOONEG)
call CHO_LSWAP(SCDIAG,XSCDIAG)

end subroutine CHO_RESETCNF
