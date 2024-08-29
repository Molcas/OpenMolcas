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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

subroutine GET_EFP()

use, intrinsic :: iso_c_binding, only: c_int
use EFP_Module, only: ABC, Coor_Type, EFP_Coors, Frag_Type, lEFP, nEFP_Coor, nEFP_Fragments
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CoorType

call Get_lScalar('EFP',lEFP)
if (lEFP) then
  call Get_iScalar('nEFP_fragments',nEFP_fragments)
  call Get_iScalar('nEFP_Coor',nEFP_Coor)
  call Get_iScalar('Coor_Type',CoorType)
  Coor_Type = int(CoorType,c_int)

  allocate(FRAG_type(nEFP_fragments))
  call Get_cArray('FRAG_Type',FRAG_Type,180*nEFP_fragments)

  allocate(ABC(3,nEFP_fragments))
  call Get_cArray('ABC',ABC,180*3*nEFP_fragments)

  allocate(EFP_Coors(nEFP_Coor,nEFP_fragments))
  call Get_dArray('EFP_COORS',EFP_COORS,nEFP_Coor*nEFP_fragments)
end if

return

end subroutine GET_EFP
