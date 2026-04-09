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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine Rd2Int_SCF()
!***********************************************************************
!                                                                      *
!     purpose: Read basis set informations from two-electron file      *
!              and compare them with those read from 2-el. file        *
!                                                                      *
!***********************************************************************

use InfSCF, only: nBas, nSkip, nSym
use Molcas, only: MxSym
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iRC, iSym, nBasX(MxSym), nSymX
logical(kind=iwp) :: SqI2

iRc = -1
call GetOrd(iRc,SqI2,nSymX,nBasX,nSkip)
if (iRc /= 0) then
  write(u6,*) 'The program failed to read the header of ORDINT.'
  call Abend()
end if

if (nSymX /= nSym) then
  write(u6,*) 'nSymX /= nSym, nSymX, nSym=',nSymX,nSym
  call Abend()
end if
do iSym=1,nSym
  if (nBas(iSym) /= nBasX(iSym)) then
    write(u6,*) 'nBas(iSym) /= nBasX(iSym)'
    write(u6,*) 'nBas=',nBas
    write(u6,*) 'nBasX=',nBasX
    call Abend()
  end if
end do

end subroutine Rd2Int_SCF
