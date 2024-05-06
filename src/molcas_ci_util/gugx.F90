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

module GUGX

use Definitions, only: wp, iwp

implicit none
private

! Split-Graph descriptor, sizes, addresses...
type SGStruct
  integer(kind=iwp) :: NSym = 0, nActEl = 0, IFRAS = 0
  integer(kind=iwp) :: IA0, IB0, IC0, iSpin, nLev, nVert, nVert0, MidLev, MVSta, MVEnd, MXUP, MXDWN, LV1RAS, LM1RAS, LV3RAS, LM3RAS
  integer(kind=iwp), allocatable :: ISm(:), DRT(:,:), DRT0(:,:), Down(:,:), Down0(:,:), Up(:,:), Ver(:), MAW(:,:), LTV(:), &
                                    DAW(:,:), RAW(:,:), SCR(:,:)
  integer(kind=iwp), pointer :: DRTP(:,:), DOWNP(:,:)
end type SGStruct

! CI Structures, addresses,..
type CIStruct
  integer(kind=iwp) :: nMidV, nIpWlk, nWalk
  integer(kind=iwp), allocatable :: NOW(:,:,:), IOW(:,:,:), NCSF(:), NOCSF(:,:,:), IOCSF(:,:,:), ICase(:), IVR(:,:), ISGM(:,:)
  real(kind=wp), allocatable :: VSGM(:,:)
end type CIStruct

! Excitation operators, coupling coefficients,...
type EXStruct
  integer(kind=iwp) :: MxEO, nICoup
  integer(kind=iwp), allocatable :: NOCP(:,:,:), IOCP(:,:,:), ICoup(:,:), MVL(:,:), MVR(:,:), USGN(:,:), LSGN(:,:)
  real(kind=wp), allocatable :: VTab(:), SGTMP(:)
end type EXStruct

type(SGStruct), target :: SGS
type(CIStruct), target :: CIS
type(EXStruct), target :: EXS

integer(kind=iwp), parameter :: MXLEV = 100
integer(kind=iwp) :: L2ACT(MXLEV), LEVEL(MXLEV)

public :: CIS, CIStruct, EXS, EXStruct, L2ACT, LEVEL, MXLEV, SGS, SGStruct

end module GUGX
