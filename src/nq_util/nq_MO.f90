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
Module nq_MO
  Integer:: nMOs=0
  Real*8, Allocatable:: CMO(:)
  Real*8, Allocatable:: D1MO(:), P2MO(:)
  Real*8, Allocatable:: P2_ontop(:,:)
End Module nq_MO
