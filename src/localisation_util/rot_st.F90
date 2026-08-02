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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************

subroutine Rot_st(cMO_s,cMO_t,nBasis,Gamma_rot,Debug)
! Author: Y. Carissan.

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBasis
real(kind=wp), intent(inout) :: cMO_s(nBasis), cMO_t(nBasis)
real(kind=wp), intent(in) :: Gamma_rot
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: iBas
real(kind=wp) :: cosGamma_rot, cs, ct, sinGamma_rot

if (Gamma_rot == Zero) return

cosGamma_rot = cos(Gamma_rot)
sinGamma_rot = sin(Gamma_rot)
if (Debug) then
  write(u6,*) 'cos(Gamma)=',cosGamma_rot
  write(u6,*) 'sin(Gamma)=',sinGamma_rot
end if

do iBas=1,nBasis
  cs = cMO_s(iBas)
  ct = cMO_t(iBas)
  cMO_s(iBas) = cosGamma_rot*cs+sinGamma_rot*ct
  cMO_t(iBas) = -sinGamma_rot*cs+cosGamma_rot*ct
end do

return

end subroutine Rot_st
