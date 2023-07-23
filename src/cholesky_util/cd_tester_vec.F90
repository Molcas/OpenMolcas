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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine CD_Tester_Vec(iVec1,nVec,Buf,lBuf,nDim,iOpt)

use CDTHLP, only: Vec
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iVec1, nVec, lBuf, nDim, iOpt
real(kind=wp) :: Buf(lBuf)
integer(kind=iwp) :: kOff, lTot
character(len=*), parameter :: SecNam = 'CD_Tester_Vec'

if (iOpt == 1) then
  kOff = 1+nDim*(iVec1-1)
  lTot = nDim*nVec
  call dCopy_(lTot,Buf,1,Vec(kOff),1)
else if (iOpt == 2) then
  kOff = 1+nDim*(iVec1-1)
  lTot = nDim*nVec
  call dCopy_(lTot,Vec(kOff),1,Buf,1)
else
  write(u6,*)
  write(u6,*) 'WARNING! WARNING! WARNING! WARNING!'
  write(u6,*) SecNam,': illegal option: iOpt = ',iOpt
  write(u6,*) 'WARNING! WARNING! WARNING! WARNING!'
  write(u6,*)
  call xFlush(u6)
end if

end subroutine CD_Tester_Vec
