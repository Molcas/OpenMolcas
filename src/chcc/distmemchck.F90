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

subroutine DistMemChck(PosV1,PosV2,PosV3,PosT)
! Dimensions for DistMemCheck
! dimensions: V1 - no*nv*nv2, nv2*nv2
!             V2 - nc*nv2
!             V3 - nc*nv*no

use Index_Functions, only: nTri_Elem
use chcc_global, only: nc, no, nv, PosFree
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: PosV1, PosV2, PosV3, PosT
integer(kind=iwp) :: len_

PosT = PosFree

PosV1 = PosT
len_ = no*nv*nTri_Elem(nv)
if (nTri_Elem(nv)**2 > len_) len_ = nTri_Elem(nv)**2
PosT = PosT+len_

PosV2 = PosT
len_ = nc*nTri_Elem(nv)
PosT = PosT+len_

PosV3 = PosT
len_ = nc*no*nv

PosT = PosT+len_

write(u6,*) ' Pos ChCk',PosV1,PosV2,PosV3,PosT

return

end subroutine DistMemChck
