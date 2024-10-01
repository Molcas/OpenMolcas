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
Module Cntrl_Data
Private
#include "Molcas.fh"
#include "cntrl.fh"
!    SONTOSTATES      Number of state pairs to calculate
Integer, Public:: SONTOSTATES=0
!    SONTO            Array of SO state pairs
     Integer, Allocatable, Public:: SONTO(:,:)
End Module Cntrl_Data
