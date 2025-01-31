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

subroutine get_diag(diag,ndet)
! Copies CI diagonal from Lucia enviroment to RASSCF envirmonet

use lucia_data, only: LUDIA
use lucia_data, only: IDISK

implicit none
real*8 diag(*)
integer nDet
integer idummy(1), iDet, IMZERO, I_AM_PACKED

ndet = 0
IDISK(LUDIA) = 0
100 continue
call IDAFILE(LUDIA,2,IDUMMY,1,IDISK(LUDIA))
IDET = IDUMMY(1)
call IDAFILE(LUDIA,2,IDUMMY,1,IDISK(LUDIA))
if (idet == -1) return
call frmdsc(diag(ndet+1),idet,-1,ludia,imzero,i_am_packed)
ndet = ndet+idet
goto 100

end subroutine get_diag
