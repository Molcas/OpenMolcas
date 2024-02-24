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
use struct, only: SGStruct
Private
Integer, Parameter, Public :: MXLEV=100

Public SGStruct
Type (SGStruct), Public:: SGS
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

Integer, Public:: ISM(MXLEV), L2ACT(MXLEV), LEVEL(MXLEV)

Integer, Allocatable, Public:: NOCSF(:), IOCSF(:), USGN(:), LSGN(:), LTV(:)
Integer,              Public::NNOCSF,   NIOCSF,                     NLTV

Integer, Allocatable, Public::  NOW1(:), IOW1(:), ICASE(:), ICOUP(:)
Integer,              Public:: NNOW,     NIOW,   NICASE,   NICOUP

Integer, Allocatable, Public, Target:: DRT(:), DOWN(:)
Integer,              Public::        NDRT,   NDOWN

Integer, Allocatable, Public:: DAW(:), UP(:), RAW(:), MAW(:)
Integer,              Public::NDAW,   NUP,   NRAW,   NMAW

Real*8,  Allocatable, Public:: VTAB(:), SGTMP(:)
Integer,              Public::NVTAB,   NSGTMP

Integer, Allocatable, Public:: MVR(:), MVL(:)
Integer,              Public::NMVR,   NMVL

Integer, Allocatable, Public:: NOCP(:), IOCP(:)
Integer,              Public::NNOCP,   NIOCP

INTEGER, Public:: NLEV,IA0,IB0,IC0,NVERT0,NVERT,                    &
                  IFCAS,LV1RAS,LM1RAS,LV3RAS,LM3RAS,                &
                  MIDLEV,NMIDV,MVSta,MVEnd,        MXUP,MXDWN,      &
                  NWALK,NIPWLK,NCSF(8),MXEO
End Module GUGX
