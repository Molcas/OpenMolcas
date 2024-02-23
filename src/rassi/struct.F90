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
module Struct
! Split-Graph descriptor, sizes, addresses...
Type SGStruct
     Integer nSym
     Integer nLev
     Integer, Allocatable:: ISm(:)
     Integer nVert
     Integer, Allocatable:: DRT(:)
     Integer, Allocatable:: Down(:)
     Integer, Allocatable:: Up(:)
     Integer MidLev
     Integer MVSta
     Integer MVEnd
     Integer, Allocatable:: MAW(:)
     Integer, Allocatable:: LTV(:)
End Type SGStruct
! CI Structures, addresses,..
Type CIStruct
     Integer nMidV
     Integer nIpWlk
     Integer, Allocatable:: NOW(:)
     Integer, Allocatable:: IOW(:)
     Integer, Allocatable:: NCSF(:)
     Integer, Allocatable:: NOCSF(:)
     Integer, Allocatable:: IOCSF(:)
     Integer nWalk
     Integer, Allocatable:: ICase(:)
End Type CIStruct
! Excitation operators, coupling coefficients,...
Type EXStruct
     Integer MxEO
     Integer, Allocatable:: NOCP(:)
     Integer, Allocatable:: IOCP(:)
     Integer nICoup
     Integer, Allocatable:: ICoup(:)
     Integer nVTab
     Integer lVTab
     Integer lMVL
     Integer lMVR
     Integer NT1MX
     Integer NT2MX
     Integer NT3MX
     Integer NT4MX
     Integer NT5MX
End Type EXStruct
Public SGStruct, CIStruct
integer, parameter ::  mxlev=100
integer :: LEVEL(mxlev)

public :: mxlev, LEVEL

end module Struct
