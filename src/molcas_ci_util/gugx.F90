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
Private
INTEGER, Public:: NLEV,IA0,IB0,IC0,NVERT0,                          &
                  IFCAS,LV1RAS,LM1RAS,LV3RAS,LM3RAS,                &
                  NVERT,                                            &
                  MIDLEV,NMIDV,MIDV1,MIDV2,        MXUP,MXDWN,      &
                  NWALK,NIPWLK,NCSF(8),NNOCSF,NIOCSF, MXEO

Integer, Parameter :: MXLEV=100
Integer, Public:: ISM(MXLEV), L2ACT(MXLEV), LEVEL(MXLEV)

Integer, Allocatable, Public:: NOCSF(:), IOCSF(:), USGN(:), LSGN(:), LTV(:)
Integer, Allocatable, Public::  NOW1(:), IOW1(:), ICASE(:), ICOUP(:)
Integer,              Public:: NNOW,     NIOW,   NICASE,   NICOUP
Integer, Allocatable, Public, Target:: DRT(:), DOWN(:)
Integer,              Public::        NDRT,   NDOWN
Integer, Allocatable, Public:: DAW(:), UP(:), RAW(:)
Integer,              Public::NDAW,   NUP,   NRAW
Real*8,  Allocatable, Public:: VTAB(:), SGTMP(:)
Integer,              Public::NVTAB,   NSGTMP
End Module GUGX
