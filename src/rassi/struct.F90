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
!     iSGStruct(1) =nSym
!     iSGStruct(2) =nLev
!     iSGStruct(3) =lISm
!     iSGStruct(4) =nVert
!     iSGStruct(5) =lDRT
!     iSGStruct(6) =lDown
!     iSGStruct(7) =lUp
!     iSGStruct(8) =MidLev
!     iSGStruct(9) =MVSta
!     iSGStruct(10)=MVEnd
!     iSGStruct(11)=lMAW
!     iSGStruct(12)=lLTV
! CI structure, sizes, addresses...
!     iCIStruct(1)=nMidV
!     iCIStruct(2)=nIpWlk
!     iCIStruct(3)=lNOW
!     iCIStruct(4)=lIOW
!     iCIStruct(5)=lNCSF
!     iCIStruct(6)=lNOCSF
!     iCIStruct(7)=lIOCSF
!     iCIStruct(8)=nWalk
!     iCIStruct(9)=lICase
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
integer, parameter :: nSGSize=12, nCISize=9, nXSize=14, mxlev=100
integer :: LEVEL(mxlev)

public :: nSGSize, nCISize, nXSize, mxlev, LEVEL

end module Struct
