************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Module EFP_Module
      use iso_c_binding
      Integer(c_int) XYZABC_Type, Points_Type, RotMat_Type
      Parameter (XYZABC_Type=0, Points_Type=1, RotMat_Type=2)
      Integer(c_int) Coor_Type
*
      Integer nEFP_fragments, nEFP_Coor
      Logical lEFP
      Character*180, Dimension(:), Allocatable:: FRAG_Type
      Character*180, Dimension(:,:), Allocatable:: ABC
      Real*8, Dimension(:,:), Allocatable, Target:: EFP_COORS
*
      type(c_ptr) :: efp_instance
      End Module EFP_Module
