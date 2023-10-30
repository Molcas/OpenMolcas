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

!***********************************************************************
!*                                                                     *
!*  CIWEIGHT  := Chirgwin-Couson/Lowdin/inverse-overlap weights of     *
!*            := full CASSCF vector and residual.                      *
!*                                                                     *
!***********************************************************************
subroutine ciweight_cvb(civec,civbs,civb,citmp,vec5,orbs,sorbs,orbinv,owrk)

use casvb_global, only: nalf, nbet, ndet, nel, norb
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: civec(0:ndet), civbs(0:ndet), civb(0:ndet), citmp(0:ndet), vec5(0:ndet)
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), intent(out) :: sorbs(norb,norb), orbinv(norb,norb), owrk(norb,norb)
integer(kind=iwp) :: ion, ionmax, ionmin, iretval1, iretval2, mxasg, mxdetcas, mxrem, mxsng, ncnfcas

ionmin = max(nel-norb,0)
ionmax = nbet
mxrem = norb-ionmin
mxsng = nel-2*ionmin
mxasg = nalf-ionmin
call icomb_cvb(mxsng,mxasg,mxdetcas)
! Work out number of configurations in CASSCF vector:
ncnfcas = 0
do ion=ionmin,ionmax
  call icomb_cvb(norb,ion,iretval1)
  call icomb_cvb(norb-ion,nel-2*ion,iretval2)
  ncnfcas = ncnfcas+iretval1*iretval2
end do

call ciweight2_cvb(civec,civbs,civb,citmp,vec5,orbs,sorbs,orbinv,owrk,ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)

return

end subroutine ciweight_cvb
