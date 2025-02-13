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

use lucia_data, only: IDISK, LUDIA
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: diag(*)
integer(kind=iwp), intent(out) :: nDet
integer(kind=iwp) :: I_AM_PACKED, iDet, idummy(1), IMZERO

ndet = 0
IDISK(LUDIA) = 0
do
  call IDAFILE(LUDIA,2,IDUMMY,1,IDISK(LUDIA))
  IDET = IDUMMY(1)
  call IDAFILE(LUDIA,2,IDUMMY,1,IDISK(LUDIA))
  if (idet == -1) exit
  call frmdsc(diag(ndet+1),idet,-1,ludia,imzero,i_am_packed)
  ndet = ndet+idet
end do

end subroutine get_diag
