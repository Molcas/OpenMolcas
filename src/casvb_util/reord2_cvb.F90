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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine reord2_cvb(cfrom,cto,imode)
! Front-end routine for molcas reord2, transforms
! from SGA CSFs to split-graph-GUGA CSFs.

use sguga, only: sg_reord
use general_data, only: nConf, STSYM
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: cfrom(nConf)
real(kind=wp), intent(out) :: cto(nConf)
integer(kind=iwp), intent(in) :: imode
integer(kind=iwp), parameter :: iState = 1

call sg_reord(iState,stsym,imode,nConf,cfrom,cto)

end subroutine reord2_cvb
