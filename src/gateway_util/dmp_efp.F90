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

subroutine DMP_EFP()

use EFP_Module, only: ABC, Coor_Type, EFP_COORS, FRAG_Type, lEFP, nEFP_Coor, nEFP_fragments
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CoorType

call Put_lScalar('EFP',lEFP)
if (lEFP) then
  call Put_iScalar('nEFP_fragments',nEFP_fragments)
  CoorType = Coor_Type
  call Put_iScalar('Coor_Type',CoorType)
  call Put_cArray('FRAG_Type',FRAG_Type(1),180*nEFP_fragments)
  call Put_cArray('ABC',ABC(1,1),3*180*nEFP_fragments)
  call Put_iScalar('nEFP_Coor',nEFP_Coor)
  call Put_dArray('EFP_COORS',EFP_COORS,nEFP_Coor*nEFP_fragments)
end if

return

end subroutine DMP_EFP
