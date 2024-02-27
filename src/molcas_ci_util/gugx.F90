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
use struct, only: SGStruct, CIStruct
Private

Public SGStruct
Type (SGStruct), Public, Target:: SGS
! Split-Graph descriptor, sizes, addresses...
!Type SGStruct
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


Integer, Parameter, Public :: MXLEV=100
Integer, Public:: L2ACT(MXLEV), LEVEL(MXLEV)

Integer, Allocatable, Public:: USGN(:), LSGN(:)

Integer, Allocatable, Public:: ICOUP(:)

Integer, Allocatable, Public:: DAW(:), RAW(:)
Integer,              Public::NDAW,    NRAW

Real*8,  Allocatable, Public:: VTAB(:), SGTMP(:)
Integer,              Public::NVTAB,   NSGTMP

Integer, Allocatable, Public:: MVR(:), MVL(:)
Integer,              Public::NMVR,   NMVL

Integer, Allocatable, Public:: NOCP(:), IOCP(:)
Integer,              Public::NNOCP,   NIOCP

INTEGER, Public:: IA0,IB0,IC0,NVERT0,                        &
                  IFCAS,LV1RAS,LM1RAS,LV3RAS,LM3RAS,         &
                  MXUP,MXDWN,MXEO
End Module GUGX
