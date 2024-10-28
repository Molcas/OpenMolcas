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

module EFP_Module

use, intrinsic :: iso_c_binding, only: c_int, c_ptr
use Definitions, only: wp, iwp

implicit none
private

integer(c_int), parameter :: XYZABC_Type = 0, Points_Type = 1, RotMat_Type = 2
integer(kind=iwp) :: nEFP_Coor, nEFP_fragments
integer(c_int) :: Coor_Type
logical(kind=iwp) :: lEFP
type(c_ptr) :: efp_instance
real(kind=wp), allocatable, target :: EFP_COORS(:,:)
character(len=180), allocatable :: ABC(:,:), FRAG_Type(:)

public :: ABC, Coor_Type, EFP_COORS, efp_instance, FRAG_Type, lEFP, nEFP_Coor, nEFP_fragments, Points_Type, RotMat_Type, XYZABC_Type

end module EFP_Module
