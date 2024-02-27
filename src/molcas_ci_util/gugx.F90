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
! Split-Graph descriptor, sizes, addresses...
!Type SGStruct
!     Integer IA0
!     Integer IB0
!     Integer IC0
!     Integer nLev
!     Integer, Allocatable:: ISm(:)
!     Integer nVert
!     Integer, Allocatable:: DRT(:)
!     Integer, Allocatable:: Down(:)
!     Integer, Allocatable:: Up(:)
!     Integer MidLev
!     Integer MVSta
!     Integer MVEnd
!     Integer, Allocatable:: MAW(:)
!     Integer, Allocatable:: LTV(:)
!     Integer, Allocatable:: RAW(:)
!     Integer, Allocatable:: DAW(:)
!     Integer MXUp
!     Integer MXDwn
!End Type SGStruct
Public CIStruct
Type (CIStruct), Public, Target:: CIS
! CI Structures, addresses,..
!Type CIStruct
!     Integer nMidV
!     Integer nIpWlk
!     Integer, Allocatable:: NOW(:)
!     Integer, Allocatable:: IOW(:)
!     Integer, Allocatable:: NCSF(:)
!     Integer, Allocatable:: NOCSF(:)
!     Integer, Allocatable:: IOCSF(:)
!     Integer nWalk
!     Integer, Allocatable:: ICase(:)
!End Type CIStruct
Public EXStruct
Type (EXStruct), Public, Target:: EXS
! Excitation operators, coupling coefficients,...
!Type EXStruct
!     Integer MxEO
!     Integer, Allocatable:: NOCP(:)
!     Integer, Allocatable:: IOCP(:)
!     Integer, Allocatable:: ICoup(:)
!     Real*8,  Allocatable:: VTab(:)
!     Integer, Allocatable:: MVL(:)
!     Integer, Allocatable:: MVR(:)
!     Integer NT1MX
!     Integer NT2MX
!     Integer NT3MX
!     Integer NT4MX
!     Integer NT5MX
!End Type EXStruct



Integer, Parameter, Public :: MXLEV=100
Integer, Public:: L2ACT(MXLEV), LEVEL(MXLEV)

Integer, Allocatable, Public:: USGN(:), LSGN(:)

Real*8,  Allocatable, Public:: SGTMP(:)
Integer,              Public:: NSGTMP

INTEGER, Public:: NVERT0,IFCAS,LV1RAS,LM1RAS,LV3RAS,LM3RAS
End Module GUGX
