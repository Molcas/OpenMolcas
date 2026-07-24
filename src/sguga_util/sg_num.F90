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

function SG_NUM(SGS,EXS,IWALK)
! PURPOSE: FOR ANY GIVEN WALK (STEP VECTOR) COMPUTE THE
!          LEXICAL NUMBER IN THE SPLIT GUGA REPRESENTATION

use sguga, only: SGStruct, EXStruct
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: SG_NUM
type (SGStruct), intent(in) :: SGS
type (EXStruct), intent(in) :: EXS
integer(kind=iwp), intent(in) ::  IWALK(SGS%NLEV)

integer(kind=iwp) :: NLEV, NVERT, MIDLEV, MVSta
integer(kind=iwp) :: IC, ICASE, ICONF, IDAWSUM, IRAWSUM, IUW, LEV, LV, MIDV

NLEV=SGS%nLev
NVERT=SGS%nVert
MIDLEV=SGS%MidLev
MVSta=SGS%MVSta

! FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY

MIDV = 1
do LEV=NLEV,(MIDLEV+1),-1
  ICASE = IWALK(LEV)
  MIDV = SGS%DOWN(MIDV,ICASE)
end do
MIDV = MIDV-MVSta+1

! FIND REVERSE ARC WEIGHT FOR THE UPPER WALK

IRAWSUM = 1
LV = 1
do LEV=NLEV,(MIDLEV+1),-1
  IC = IWALK(LEV)
  LV = SGS%DOWN(LV,IC)
  IRAWSUM = IRAWSUM+SGS%RAW(LV,IC)
end do
IUW = EXS%USGN(IRAWSUM,MIDV)

! FIND DIRECT ARC WEIGHT FOR THE LOWER WALK

IDAWSUM = 1
LV = NVERT
do LEV=1,MIDLEV
  IC = IWALK(LEV)
  LV = SGS%UP(LV,IC)
  IDAWSUM = IDAWSUM+SGS%DAW(LV,IC)
end do
ICONF = EXS%LSGN(IDAWSUM,MIDV)

! COMPUTE LEXICAL ORDERING NUMBER

SG_NUM = IUW+ICONF

end function SG_NUM
