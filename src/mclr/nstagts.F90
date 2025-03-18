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

function NSTAGTS(IGRP,ITP,ISM)
! Number of strings of group IGRP
!                      type  ITP
!                      sym   ISM

use Str_Info, only: NOCTYP, STR

implicit none
integer IGRP, ITP, ISM
integer IADDRESS, NSTAGTS, NTEST

! element (ITP,ISM) corresponds to address
IADDRESS = (ISM-1)*NOCTYP(IGRP)+ITP

NSTAGTS = Str(IGRP)%NSTSO(IADDRESS)

NTEST = 0
if (NTEST /= 0) write(6,*) ' Number of strings of group,type,sym',IGRP,ITP,ISM,' is ',NSTAGTS

end function NSTAGTS
