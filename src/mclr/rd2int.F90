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
! Copyright (C) 1993, Johan Lorentzon                                  *
!               1993, Jeppe Olsen                                      *
!               1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Rd2Int(iPL)
!***********************************************************************
!                                                                      *
!     Read header of the two-electron integral file                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     J. Lorentzon, J. Olsen and M.P. Fuelscher                        *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use input_mclr, only: CasInt, nBas, nSkip, nSym, TimeDep
use Molcas, only: MxSym
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPL
integer(kind=iwp) :: iRC, iSym, nBasX(mxSym), nSymX, ntSkip
logical(kind=iwp) :: SqSym

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
iRc = -1
call GetOrd(iRc,SqSym,nSymX,nBasX,nSkip)
if (iRc /= 0) then
  write(u6,*) 'Rd2Int: Error reading ORDINT'
  call Abend()
end if
if (iPL >= 2) then
  if (SqSym) write(u6,*) 'OrdInt status: squared'
  if (.not. SqSym) write(u6,*) 'OrdInt status: non-squared'
end if
if (nSymX /= nSym) then
  write(u6,*) 'Rd2Int: nSymX /= nSym'
  write(u6,*) 'nSymX,nSym=',nSymX,nSym
  call Abend()
end if
do iSym=1,nSym
  if (nBas(iSym) /= nBasX(iSym)) then
    write(u6,*) 'Rd2Int: nBas(iSym) /= nBasX(iSym)'
    write(u6,*) 'nBas(iSym),nBasX(iSym)=',nBas(iSym),nBasX(iSym)
    call Abend()
  end if
end do
ntSkip = sum(nSkip(1:nSym))
if (ntSkip /= 0) then
  write(u6,*) 'Rd2Int: ntSkip /= 0'
  write(u6,*) 'ntSkip=',ntSkip
  call Abend()
end if
if ((.not. SqSym) .and. (.not. TimeDep)) then
  CASINT = .true.
else
  CASINT = .false.
end if
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine Rd2Int
