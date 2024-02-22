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
Type CIStruct
     Integer nMidV
     Integer nIpWlk
     Integer lNOW
     Integer lIOW
     Integer lNCSF
     Integer lNOCSF
     Integer LIOCSF
     Integer nWalk
     Integer lICase
End Type CIStruct
! Excitation operators, coupling coefficients,...
!     iXStruct(1)=MxEO
!     iXStruct(2)=lNOCP
!     iXStruct(3)=lIOCP
!     iXStruct(4)=nICoup
!     iXStruct(5)=lICoup
!     iXStruct(6)=nVTab
!     iXStruct(7)=lVTab
!     iXStruct(8 )=lMVL
!     iXStruct(9 )=lMVR
!     iXStruct(10)=NT1MX
!     iXStruct(11)=NT2MX
!     iXStruct(12)=NT3MX
!     iXStruct(13)=NT4MX
!     iXStruct(14)=NT5MX
! Dimensions of structures.
Public SGStruct, CIStruct
integer, parameter :: nCISize=9, nXSize=14, mxlev=100
integer :: LEVEL(mxlev)

public :: nCISize, nXSize, mxlev, LEVEL

end module Struct
