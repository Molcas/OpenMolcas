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
     Integer :: NSym=0
     Integer IA0
     Integer IB0
     Integer IC0
     Integer nLev
     Integer, Allocatable:: ISm(:)
     Integer nVert
     Integer nVert0
     Integer, Allocatable :: DRT(:,:)
     Integer, Allocatable :: DRT0(:,:)
     Integer, Pointer::      DRTP(:,:)
     Integer, Allocatable :: Down(:,:)
     Integer, Allocatable :: Down0(:,:)
     Integer, Pointer::      DOWNP(:,:)
     Integer, Allocatable:: Up(:,:)
     Integer, Allocatable:: Ver(:)
     Integer MidLev
     Integer MVSta
     Integer MVEnd
     Integer, Allocatable:: MAW(:,:)
     Integer, Allocatable:: LTV(:)
     Integer, Allocatable:: DAW(:,:)
     Integer, Allocatable:: RAW(:,:)
     Integer MXUP
     Integer MXDWN
     Integer LV1RAS
     Integer LM1RAS
     Integer LV3RAS
     Integer LM3RAS
     Integer, Allocatable:: SCR(:,:)
End Type SGStruct
! CI Structures, addresses,..
Type CIStruct
     Integer nMidV
     Integer nIpWlk
     Integer, Allocatable:: NOW(:,:,:)
     Integer, Allocatable:: IOW(:,:,:)
     Integer, Allocatable:: NCSF(:)
     Integer, Allocatable:: NOCSF(:,:,:)
     Integer, Allocatable:: IOCSF(:,:,:)
     Integer nWalk
     Integer, Allocatable:: ICase(:)
End Type CIStruct
! Excitation operators, coupling coefficients,...
Type EXStruct
     Integer MxEO
     Integer, Allocatable:: NOCP(:,:,:)
     Integer, Allocatable:: IOCP(:,:,:)
     Integer, Allocatable:: ICoup(:,:)
     Real*8,  Allocatable:: VTab(:)
     Integer, Allocatable:: MVL(:,:)
     Integer, Allocatable:: MVR(:,:)
     Real*8,  Allocatable:: SGTMP(:)
     Integer, Allocatable:: USGN(:,:)
     Integer, Allocatable:: LSGN(:,:)
End Type EXStruct
Public SGStruct, CIStruct, EXStruct

end module Struct
