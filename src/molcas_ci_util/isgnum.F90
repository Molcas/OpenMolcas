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

function ISGNUM(NLEV,NVERT,MIDLEV,MVSta,NMIDV,MXUP,MXDWN,IDOWN,IUP,IDAW,IRAW,IUSGNUM,ILSGNUM,IWALK)
! PURPOSE: FOR ANY GIVEN WALK (STEP VECTOR) COMPUTE THE
!          LEXICAL NUMBER IN THE SPLIT GUGA REPRESENTATION

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISGNUM
integer(kind=iwp), intent(in) :: NLEV, NVERT, MIDLEV, MVSta, NMIDV, MXUP, MXDWN, IDOWN(NVERT,0:3), IUP(NVERT,0:3), &
                                 IDAW(NVERT,0:4), IRAW(NVERT,0:4), IUSGNUM(MXUP,NMIDV), ILSGNUM(MXDWN,NMIDV), IWALK(NLEV)
integer(kind=iwp) :: IC, ICASE, ICONF, IDAWSUM, IRAWSUM, IUW, LEV, LV, MIDV

! FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY

MIDV = 1
do LEV=NLEV,(MIDLEV+1),-1
  ICASE = IWALK(LEV)
  MIDV = IDOWN(MIDV,ICASE)
end do
MIDV = MIDV-MVSta+1

! FIND REVERSE ARC WEIGHT FOR THE UPPER WALK

IRAWSUM = 1
LV = 1
do LEV=NLEV,(MIDLEV+1),-1
  IC = IWALK(LEV)
  LV = IDOWN(LV,IC)
  IRAWSUM = IRAWSUM+IRAW(LV,IC)
end do
IUW = IUSGNUM(IRAWSUM,MIDV)

! FIND DIRECT ARC WEIGHT FOR THE LOWER WALK

IDAWSUM = 1
LV = NVERT
do LEV=1,MIDLEV
  IC = IWALK(LEV)
  LV = IUP(LV,IC)
  IDAWSUM = IDAWSUM+IDAW(LV,IC)
end do
ICONF = ILSGNUM(IDAWSUM,MIDV)

! COMPUTE LEXICAL ORDERING NUMBER

ISGNUM = IUW+ICONF

! EXIT

return

end function ISGNUM
