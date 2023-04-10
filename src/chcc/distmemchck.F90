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

implicit none
#include "chcc1.fh"
integer PosV1, PosV2, PosV3, PosT
! help variables
integer len_

PosT = PosFree

PosV1 = PosT
len_ = no*nv*nv*(nv+1)/2
if ((nv*(nv+1)*nv*(nv+1)/4) > len_) len_ = nv*(nv+1)*nv*(nv+1)/4
PosT = PosT+len_

PosV2 = PosT
len_ = nc*nv*(nv+1)/2
PosT = PosT+len_

PosV3 = PosT
len_ = nc*no*nv

PosT = PosT+len_

write(6,*) ' Pos ChCk',PosV1,PosV2,PosV3,PosT

return

end subroutine DistMemChck
