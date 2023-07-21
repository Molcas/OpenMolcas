!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_X_Init_Par_DF(irc)
!
! Purpose: setup for parallel DF.

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs, Is_Real_Par
#endif

implicit none
integer irc
character(len=17), parameter :: SecNam = 'Cho_X_Init_Par_DF'
#ifdef _DEBUGPRINT_
logical, parameter :: LocDbg = .true.
#else
logical, parameter :: LocDbg = .false.
#endif
#ifdef _MOLCAS_MPP_
#include "cholesky.fh"
integer nV(8)
integer iSym
logical isSerial

irc = 0

! Return if serial.
! -----------------

isSerial = (nProcs == 1) .or. (.not. Is_Real_Par())
if (isSerial) then
  if (LocDbg) then
    write(6,*) SecNam,': serial run, nothing to do...'
    write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
  end if
  return
else
  if (LocDbg) then
    write(6,*) SecNam,': parallel run...'
    write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
  end if
end if

! Reset number of vectors to the number on this node as stored on
! the runfile.
! ---------------------------------------------------------------

call iCopy(nSym,NumCho,1,nV,1)
call Get_iArray('nVec_RI',NumCho,nSym)
NumChT = NumCho(1)
do iSym=2,nSym
  NumChT = NumChT+NumCho(iSym)
end do

! Debug print.
! ------------

if (LocDbg) then
  write(6,*)
  write(6,*) 'Output from ',SecNam,':'
  write(6,*) 'NumCho before: ',(nV(iSym),iSym=1,nSym)
  write(6,*) 'NumCho after : ',(NumCho(iSym),iSym=1,nSym)
end if

#else

irc = 0
if (LocDbg) write(6,*) SecNam,': serial run, nothing to do...'

#endif

end subroutine Cho_X_Init_Par_DF
