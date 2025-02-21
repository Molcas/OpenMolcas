!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
function ICHECK_OCC_IN_ACCSPC(IOCC,IMINMAX,NGAS,MXPNGAS)
! Check if Occupation of GAS Spaces defined by IOCC are
! within the constraints of IMINMAX chosen by the user
!
! Giovanni Li Manni 7 Nov 2011, for BK implementation

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: ICHECK_OCC_IN_ACCSPC
integer(kind=iwp), intent(in) :: NGAS, IOCC(NGAS), MXPNGAS, IMINMAX(MXPNGAS,2)
integer(kind=iwp) :: I_AM_IN

I_AM_IN = 1
if (any(IOCC(:) < IMINMAX(1:NGAS,1))) I_AM_IN = 0
if (any(IOCC(:) > IMINMAX(1:NGAS,2))) I_AM_IN = 0
ICHECK_OCC_IN_ACCSPC = I_AM_IN

#ifdef _DEBUGPRINT_
write(u6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IMINMAX'
call IWRTMA(IMINMAX,NGAS,2,MXPNGAS,2)
write(u6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IOCC'
call IWRTMA(IOCC,1,NGAS,1,NGAS)
write(u6,*) ' And the verdict is ',I_AM_IN
#endif

end function ICHECK_OCC_IN_ACCSPC
