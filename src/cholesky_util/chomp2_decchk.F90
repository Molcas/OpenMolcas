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
! Copyright (C) 2008, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_DecChk(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
!
! Thomas Bondo Pedersen, Jan. 2008.
!
! Purpose: check Cholesky decomposition of the (ai|bj) integrals
!          or MP2 amplitudes (sym. block iSym).
!          The columns of the matrix are
!          compared nCol columns at a time. This
!          implies that the memory requirement of this routine
!          should be limited to approximately that of the
!          decomposition itself. Note, however, that since all
!          integrals are computed, this routine will consume
!          significantly more CPU time.
!          Files are assumed open.
!          On exit,
!          ErrStat(1) = min error
!          ErrStat(2) = max error
!          ErrStat(3) = rms error

use ChoMP2, only: iOption_MP2CD
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iSym, nDim, nCol, lWrk
real(kind=wp), intent(inout) :: Col(nDim,nCol)
real(kind=wp), intent(out) :: Wrk(lWrk), ErrStat(3)
character(len=*), parameter :: SecNam = 'ChoMP2_DecChk'

if (iOption_MP2CD == 1) then ! (ai|bj) int
  call ChoMP2_DecChk_1(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
else if (iOption_MP2CD == 2) then ! MP2 amp
  call ChoMP2_DecChk_2(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
else
  write(u6,*) SecNam,': WARNING! Unknown option, iOption_MP2CD = ',iOption_MP2CD
  irc = -123456
end if

end subroutine ChoMP2_DecChk
