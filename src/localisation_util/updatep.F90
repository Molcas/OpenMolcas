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
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine UpdateP(PACol,BName,nBas_Start,nOrb2Loc,nAtoms,PA,gamma_rot,iMO_s,iMO_t,Debug)
! Author: Yannick Carissan.
!
! Modifications:
!    - October 6, 2005 (Thomas Bondo Pedersen):
!      Reduce operation count and use BLAS.

use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nAtoms, nBas_Start(nAtoms), nOrb2Loc, iMO_s, iMO_t
real(kind=wp), intent(out) :: PACol(nOrb2Loc,2)
character(len=LenIn8), intent(in) :: BName(*)
real(kind=wp), intent(inout) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(in) :: gamma_rot
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: iAt
#ifdef _DEBUGPRINT_
real(kind=wp) :: PA_ts, Tst
#endif
real(kind=wp) :: cos2g, cosg, cosing, PA_ss, PA_st, PA_tt, sin2g, sing
character(len=LenIn8) :: PALbl

cosg = cos(gamma_rot)
sing = sin(gamma_rot)
cos2g = cosg*cosg
sin2g = sing*sing
cosing = cosg*sing

do iAt=1,nAtoms
  !call RecPrt('PA(1,1,iAt)',' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)

  ! Copy out the PA_ss, PA_tt, and PA_st elements.

  PA_ss = PA(iMO_s,iMO_s,iAt)
  PA_st = PA(iMO_s,iMO_t,iAt)
  PA_tt = PA(iMO_t,iMO_t,iAt)
  !write(u6,*) 'updateP:',PA_ss,PA_st,PA_tt
# ifdef _DEBUGPRINT_
  PA_ts = PA(iMO_t,iMO_s,iAt)
  Tst = PA_st-PA_ts
  if (abs(Tst) > 1.0e-14_wp) then
    write(u6,*) 'Broken symmetry in UpdateP!!'
    write(u6,*) 'MOs s and t: ',iMO_s,iMO_t
    write(u6,*) 'PA_st = ',PA_st
    write(u6,*) 'PA_ts = ',PA_ts
    write(u6,*) 'Diff = ',Tst
    call SysAbendMsg('UpdateP','Broken symmetry!',' ')
  end if
# endif

  ! Copy out columns s and t of PA.

  PACol(:,1) = PA(:,iMO_s,iAt)
  PACol(:,2) = PA(:,iMO_t,iAt)

  ! Compute transformed columns.

  PA(:,iMO_s,iAt) = cosg*PACol(:,1)+sing*PACol(:,2)
  PA(:,iMO_t,iAt) = cosg*PACol(:,2)-sing*PACol(:,1)

  ! Compute PA_ss, PA_tt, PA_st, and PA_ts (= PA_st).

  PA(iMO_s,iMO_s,iAt) = PA_ss*cos2g+PA_tt*sin2g+Two*PA_st*cosing
  PA(iMO_t,iMO_s,iAt) = (PA_tt-PA_ss)*cosing+PA_st*(cos2g-sin2g)
  PA(iMO_s,iMO_t,iAt) = PA(iMO_t,iMO_s,iAt)
  PA(iMO_t,iMO_t,iAt) = PA_tt*cos2g+PA_ss*sin2g-Two*PA_st*cosing

  ! Copy columns to rows.

  PA(iMO_s,:,iAt) = PA(:,iMO_s,iAt)
  PA(iMO_t,:,iAt) = PA(:,iMO_t,iAt)

end do

if (Debug) then
  write(u6,*) 'In UpdateP'
  write(u6,*) '----------'
  do iAt=1,nAtoms
    PALbl = 'PA__'//BName(nBas_Start(iAt))(1:LenIn)
    call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
  end do
end if

return

end subroutine UpdateP
