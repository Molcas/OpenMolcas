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
! Copyright (C) 1990,1991,1993,1996, Roland Lindh                      *
!               1990, IBM                                              *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine Close_SemiDSCF()

use IOBUF, only: iBuf, iPos, iStatIO, lStRec, Mode_None, OnDisk

implicit none

!write(u6,*) 'Enter: Close_SemiDSCF'

! If data was transfered to the I/O buffer write buffer on disc.

! If buffer empty force the write :
if (iPos == 1) iPos = 2
if (OnDisk) call WLBuf()

iPos = lStRec+1
iStatIO = Mode_None
iBuf = -99
!write(u6,*) 'Exit: Close_SemiDSCF'

end subroutine Close_SemiDSCF
