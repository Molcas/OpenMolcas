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
Module Intgrl
use Definitions, only: iwp
Implicit None
Private iwp
INTEGER(kind=iwp), public:: IF_ERI,IAD2M(3,36*36),NSYMZ,NORBZ(8),NOSHZ(8),LUINTMZ,IL_ERI
ENd Module Intgrl
