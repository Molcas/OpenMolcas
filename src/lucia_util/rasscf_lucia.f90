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
Module RASSCF_LUCIA
Private
INTEGER, Public:: kvec3_length=0, ini_h0,Memory_Needed_Lucia=0
Logical, Public:: Sigma_on_disk=.FALSE.
Real*8, Allocatable, Public:: CIVec(:)
Real*8, Allocatable, Public:: PAtmp(:)
Real*8, Allocatable, Public:: Pscr(:)
Real*8, Allocatable, Public:: Ptmp(:)
Real*8, Allocatable, Public:: DStmp(:)
Real*8, Allocatable, Public:: Dtmp(:)
Real*8, Allocatable, Public:: RF1(:)
Real*8, Allocatable, Public:: RF2(:)
End Module RASSCF_LUCIA
