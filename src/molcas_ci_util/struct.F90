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
     Integer IA0
     Integer IB0
     Integer IC0
     Integer nLev
     Integer, Allocatable:: ISm(:)
     Integer nVert
     Integer, Allocatable :: DRT(:)
     Integer, Allocatable :: Down(:)
     Integer, Allocatable:: Up(:)
     Integer MidLev
     Integer MVSta
     Integer MVEnd
     Integer, Allocatable:: MAW(:)
     Integer, Allocatable:: LTV(:)
     Integer, Allocatable:: DAW(:)
     Integer, Allocatable:: RAW(:)
     Integer MXUP
     Integer MXDWN
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
     Integer, Allocatable:: ICoup(:)
     Real*8,  Allocatable:: VTab(:)
     Integer, Allocatable:: MVL(:)
     Integer, Allocatable:: MVR(:)
     Integer NT1MX
     Integer NT2MX
     Integer NT3MX
     Integer NT4MX
     Integer NT5MX
End Type EXStruct
Public SGStruct, CIStruct, EXStruct
integer, parameter ::  mxlev=100
integer :: LEVEL(mxlev)

public :: mxlev, LEVEL

end module Struct
