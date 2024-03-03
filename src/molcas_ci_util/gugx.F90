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
Module GUGX
use struct, only: SGStruct, CIStruct, ExStruct
Private

Public SGStruct
Type (SGStruct), Public, Target:: SGS
Public CIStruct
Type (CIStruct), Public, Target:: CIS
Public EXStruct
Type (EXStruct), Public, Target:: EXS


Integer, Parameter, Public :: MXLEV=100
Integer, Public:: L2ACT(MXLEV), LEVEL(MXLEV)

INTEGER, Public:: IFRAS

End Module GUGX
