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

subroutine CHO_PUTRED(IPASS,IRED)
!
! Purpose: write reduced set indices to disk and set address for
!          next write.

use Cholesky, only: IndRed, IndRSh, InfRed, iSP2F, LuPri, MaxRed, mmBstRT, nnBstRSh, nnBstRT, nnShl, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IPASS, IRED
character(len=*), parameter :: SECNAM = 'CHO_PUTRED'

if (IPASS > MAXRED) then
  write(LUPRI,*) SECNAM,': integral pass ',IPASS
  write(LUPRI,*) SECNAM,': max. allowed is ',MAXRED
  write(LUPRI,*) SECNAM,': please increase max. allowed!'
  call CHO_QUIT('Too many integral passes in '//SECNAM,104)
else if (IPASS == 1) then
  call CHO_PUTRED1(INFRED,nnBstRSh(:,:,1),IndRed(:,1),INDRSH,iSP2F,MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,1)
  if (MAXRED > 1) INFRED(IPASS+1) = INFRED(IPASS)+NSYM*NNSHL+2*NNBSTRT(1)+NNSHL
else if (IPASS == MAXRED) then
  call CHO_PUTRED1(INFRED,nnBstRSh(:,:,IRED),IndRed(:,IRED),INDRSH,iSP2F,MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,IRED)
else
  call CHO_PUTRED1(INFRED,nnBstRSh(:,:,IRED),IndRed(:,IRED),INDRSH,iSP2F,MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,IRED)
  INFRED(IPASS+1) = INFRED(IPASS)+NSYM*NNSHL+NNBSTRT(IRED)
end if

end subroutine CHO_PUTRED
