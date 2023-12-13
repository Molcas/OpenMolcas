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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine ChoLoc_xp(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc,iD)
! Same as ChoLoc_p but handles differently the irc=102

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nBas, nOcc
real(kind=wp), intent(inout) :: Dens(nBas,nBas)
real(kind=wp), intent(out) :: CMO(nBas,nOcc), xNrm
real(kind=wp), intent(in) :: Thrs
integer(kind=iwp), intent(inout) :: iD(nBas)
integer(kind=iwp) :: nVec
character(len=*), parameter :: SecNam = 'ChoLoc_xp'
real(kind=wp), external :: ddot_

irc = 0
xNrm = -huge(xNrm)

nVec = 0
call CD_InCore_p(Dens,nBas,CMO,nOcc,iD,nVec,Thrs,irc)
if (irc == 102) then
  irc = 0  ! reset because it is most likely a numerical noise
else if (irc /= 0) then
  write(u6,*) SecNam,': CD_InCore_p returned ',irc
  return
else if (nVec /= nOcc) then
  write(u6,*) SecNam,': nVec /= nOcc'
  write(u6,*) '   nVec,nOcc = ',nVec,nOcc
  irc = 1
  return
end if

xNrm = sqrt(dDot_(nBas*nOcc,CMO,1,CMO,1))

end subroutine ChoLoc_xp
