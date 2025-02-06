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

function ICHECK_OCC_IN_ACCSPC(IOCC,IMINMAX,NGAS,MXPNGAS)
! Check if Occupation of GAS Spaces defined by IOCC are
! within the constraints of IMINMAX chosen by the user
!
! Giovanni Li Manni 7 Nov 2011, for BK implementation

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ICHECK_OCC_IN_ACCSPC
integer(kind=iwp) :: NGAS, IOCC(NGAS), MXPNGAS, IMINMAX(MXPNGAS,2)
integer(kind=iwp) :: I_AM_IN, IGAS, NTEST

I_AM_IN = 1
do IGAS=1,NGAS
  if ((IOCC(IGAS) < IMINMAX(IGAS,1)) .or. (IOCC(IGAS) > IMINMAX(IGAS,2))) I_AM_IN = 0
end do
ICHECK_OCC_IN_ACCSPC = I_AM_IN

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IMINMAX'
  call IWRTMA(IMINMAX,NGAS,2,MXPNGAS,2)
end if
if (NTEST >= 10) then
  write(u6,*) ' Input to ICHECK_OCC_IN_ACCSPC, IOCC'
  call IWRTMA(IOCC,1,NGAS,1,NGAS)
  write(u6,*) ' And the verdict is ',I_AM_IN
end if

end function ICHECK_OCC_IN_ACCSPC
