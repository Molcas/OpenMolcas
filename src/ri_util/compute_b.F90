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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

function Compute_B(irc,kSOk,lSOl,jAOj,nBasFnc,iOpt)
!***********************************************************************
!                                                                      *
!     Author Jonas Bostrom, June 2010                                  *
!                                                                      *
!     Purpose: To do part of MP2 gradient.                             *
!                                                                      *
!***********************************************************************

use RI_glob, only: BMP2
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Compute_B
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: kSOk, lSOl, jAOj, nBasFnc, iOpt
integer(kind=iwp) :: iOff1, iOff2

iOff1 = jAOj*nBasFnc*nBasFnc+(kSOk-1)*nBasFnc+lSOl
iOff2 = jAOj*nBasFnc*nBasFnc+(lSOl-1)*nBasFnc+kSOk
Compute_B = (Bmp2(iOff1,iOpt)+Bmp2(iOff2,iOpt))*Half
irc = 0

end function Compute_B
