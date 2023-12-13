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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_OpenB(iOpt,iSym,iBatch)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: open (iOpt=1), close and keep (iOpt=2), or close and
!          delete (iOpt=3) Cholesky vector files for MP2 program
!          (batch vectors).
!          For iOpt=0, the units are initialized (to -1).

use ChoMP2, only: LnT1am, lUnit
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iOpt, iSym, iBatch
integer(kind=iwp) :: lU
character(len=6) :: BtchNm
character(len=*), parameter :: BaseNm = '_I', SecNam = 'ChoMP2_OpenB'

! Initialize units and return for iOpt=0.
! ---------------------------------------

if (iOpt == 0) then
  lUnit(iSym,iBatch) = -1
  return
end if

! Open or close files.
! --------------------

if (iOpt == 1) then
  if (LnT1am(iSym,iBatch) > 0) then
    if (iBatch < 10) then
      write(BtchNm,'(A2,I1,A2,I1)') BaseNm,iSym,'__',iBatch
    else if (iBatch < 100) then
      write(BtchNm,'(A2,I1,A1,I2)') BaseNm,iSym,'_',iBatch
    else if (iBatch < 1000) then
      write(BtchNm,'(A2,I1,I3)') BaseNm,iSym,iBatch
    else ! note: due to restriction in filename length...
      call SysAbendMsg(SecNam,'Too many batches','(Current max. is 999)')
      BtchNm = '?!?!?!' ! too avoid compiler warnings...
    end if
    lU = 7
    call daName_MF_WA(lU,BtchNm)
  else
    lU = -1
  end if
  lUnit(iSym,iBatch) = lU
else if (iOpt == 2) then
  lU = lUnit(iSym,iBatch)
  if (lU > 0) then
    call daClos(lU)
    lUnit(iSym,iBatch) = -1
  end if
else if (iOpt == 3) then
  lU = lUnit(iSym,iBatch)
  if (lU > 0) then
    call daEras(lU)
    lUnit(iSym,iBatch) = -1
  end if
else
  call SysAbendMsg(SecNam,'iOpt out of bounds',' ')
end if

end subroutine ChoMP2_OpenB
