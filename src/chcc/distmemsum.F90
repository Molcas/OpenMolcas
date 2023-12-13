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

subroutine DistMemSum(maxdim,PosV1,PosV2,PosV3,PosH1,PosH2,PosT)
! This routine does:
! define initial positions of V and H
! used in summary routine

use chcc_global, only: no, printkey
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: maxdim
integer(kind=iwp), intent(out) :: PosV1, PosV2, PosV3, PosH1, PosH2
integer(kind=iwp), intent(inout) :: PosT
integer(kind=iwp) :: length

!1 V files

length = no*no*maxdim*maxdim
PosV1 = PosT
PosT = PosT+length
PosV2 = PosT
PosT = PosT+length
PosV3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,99) 'DM V ',PosV1,PosV2,PosV3,length

!2 H files

length = no*maxdim
PosH1 = PosT
PosT = PosT+length
PosH2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,99) 'DM H ',PosH1,PosH2,PosV3,length

if (printkey >= 10) write(u6,99) 'PosT ',PosT

return

99 format(a7,10(i10,1x))

end subroutine DistMemSum
