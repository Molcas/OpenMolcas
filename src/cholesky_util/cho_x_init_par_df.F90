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
use Para_Info, only: Is_Real_Par, nProcs
#ifdef _DEBUGPRINT_
use Para_Info, only: MyRank
#endif
use Cholesky, only: nSym, NumCho, NumChT
#endif
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: irc
character(len=*), parameter :: SecNam = 'Cho_X_Init_Par_DF'
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: nV(8)
logical(kind=iwp) :: isSerial
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iSym
#endif

irc = 0

! Return if serial.
! -----------------

isSerial = (nProcs == 1) .or. (.not. Is_Real_Par())
if (isSerial) then
# ifdef _DEBUGPRINT_
  write(u6,*) SecNam,': serial run, nothing to do...'
  write(u6,*) '#nodes: ',nProcs,'  myRank: ',myRank
# endif
  return
# ifdef _DEBUGPRINT_
else
  write(u6,*) SecNam,': parallel run...'
  write(u6,*) '#nodes: ',nProcs,'  myRank: ',myRank
# endif
end if

! Reset number of vectors to the number on this node as stored on
! the runfile.
! ---------------------------------------------------------------

nV(1:nSym) = NumCho(1:nSym)
call Get_iArray('nVec_RI',NumCho,nSym)
NumChT = sum(NumCho(1:nSym))

#ifdef _DEBUGPRINT_
! Debug print.
! ------------

write(u6,*)
write(u6,*) 'Output from ',SecNam,':'
write(u6,*) 'NumCho before: ',(nV(iSym),iSym=1,nSym)
write(u6,*) 'NumCho after : ',(NumCho(iSym),iSym=1,nSym)
#endif

#else

irc = 0
#ifdef _DEBUGPRINT_
write(u6,*) SecNam,': serial run, nothing to do...'
#endif

#endif

end subroutine Cho_X_Init_Par_DF
