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
Module Temperatures
use Definitions, only: wp, iwp
Private
! Default temperatures for thermochemistry (MCLR, SLAPAF)
Integer(kind=iwp), Parameter :: NDefTemp=7
Real(kind=wp) :: DefTemp(NDefTemp)=[0.0d0,100.0d0,273.15d0,298.15d0,323.15d0,373.15d0,473.15d0]
Public :: DefTemp
End Module Temperatures
