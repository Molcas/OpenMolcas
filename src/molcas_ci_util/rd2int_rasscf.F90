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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Rd2Int_RASSCF()
!***********************************************************************
!                                                                      *
!     Read header of the two-electron integral file                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use UnixInfo, only: ProgName
use rasscf_global, only: lSquare
use Definitions, only: iwp, u6

implicit none
#include "rasdim.fh"
#include "general.fh"
integer(kind=iwp) :: i, IERR, iRc, iSym, nBasX(8), nSymX

!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
iRc = -1
call GetOrd(iRc,lSquare,nSymX,nBasX,nSkipX)
if (iRc /= 0) then
  write(u6,*) 'RD2INT Error: Failed to read from ORDINT file.'
  write(u6,*) Progname,' tried to read two-electron integrals from'
  write(u6,*) 'the ORDINT file, but failed. Something is wrong'
  write(u6,*) 'with the file. Perhaps it is missing?'
  call Quit_OnUserError()
end if
if (nSymX /= nSym) then
  write(u6,*) 'RD2INT Error: Wrong size of symmetry group.'
  write(u6,*) ProgName,' tried to use two-electron integrals from'
  write(u6,*) 'a file that was evidently created for some other'
  write(u6,*) 'program run.'
  write(u6,'(1x,a,2i8)') 'nSymX,nSym:',nSymX,nSym
  call Quit_OnUserError()
end if
IERR = 0
do iSym=1,nSym
  if (nBas(iSym) /= nBasX(iSym)) then
    IERR = 1
    exit
  end if
end do
if (IERR == 1) then
  write(u6,*) 'RD2INT Error: Wrong nr of basis functions.'
  write(u6,*) 'RASSCF tried to use two-electron integrals from'
  write(u6,*) 'a file that was evidently created for some other'
  write(u6,*) 'program run.'
  write(u6,'(1x,a,8i8)') 'nBas :',(nBas(i),i=1,nSym)
  write(u6,'(1x,a,8i8)') 'nBasX:',(nBasX(i),i=1,nSym)
  call Quit_OnUserError()
end if
!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*

end subroutine Rd2Int_RASSCF
