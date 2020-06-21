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
      Module External_Centers
      Integer :: nEF=0
      Real*8, Allocatable:: EF_Centers(:,:)
      Real*8, Allocatable:: OAM_Center(:)
      Real*8, Allocatable:: OMQ_Center(:)
      Integer :: nDMS=0
      Real*8, Allocatable:: DMS_Centers(:,:)
      Integer :: nWel=0
      Real*8, Allocatable:: Wel_Info(:,:)
      Real*8, Allocatable:: AMP_Center(:)
      Integer :: nRP=0
      Real*8, Target, Allocatable:: RP_Centers(:,:,:)
      End Module External_Centers
