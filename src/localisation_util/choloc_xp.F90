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

implicit none
integer irc, nBas, nOcc, iD(nBas)
real*8 Thrs, xNrm
real*8 Dens(nBas,nBas), CMO(nBas,nOcc)
character*9 SecNam
parameter(SecNam='ChoLoc_xp')
integer nVec
real*8 ddot_
external ddot_

irc = 0
xNrm = -9.9d9

nVec = 0
call CD_InCore_p(Dens,nBas,CMO,nOcc,iD,nVec,Thrs,irc)
if ((irc /= 0) .and. (irc /= 102)) then
  write(6,*) SecNam,': CD_InCore_p returned ',irc
  return
else if (irc == 102) then
  irc = 0  ! reset because it is most likely a numerical noise
else if (nVec /= nOcc) then
  write(6,*) SecNam,': nVec /= nOcc'
  write(6,*) '   nVec,nOcc = ',nVec,nOcc
  irc = 1
  return
end if

xNrm = sqrt(dDot_(nBas*nOcc,CMO,1,CMO,1))

end subroutine ChoLoc_xp
