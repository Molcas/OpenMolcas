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

integer(c_int) XYZABC_Type, Points_Type, RotMat_Type
parameter(XYZABC_Type=0,Points_Type=1,RotMat_Type=2)
integer(c_int) Coor_Type
integer nEFP_fragments, nEFP_Coor
logical lEFP
character*180, dimension(:), allocatable :: FRAG_Type
character*180, dimension(:,:), allocatable :: ABC
real*8, dimension(:,:), allocatable, target :: EFP_COORS
type(c_ptr) :: efp_instance

end module EFP_Module
