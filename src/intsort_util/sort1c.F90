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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine SORT1C(nUt,vInt,nSqNum,nSyBlk)
!***********************************************************************
!                                                                      *
!     Purpose: Provided SEWARD was allowed to use a virtual disk,      *
!              scatter the 2el integrals into the approriate place.    *
!                                                                      *
!     Called from: SORT1A                                              *
!                                                                      *
!     Calling Parameters:                                              *
!     nUt    : number of 2el integrals in the buffers vInt,nSqNum      *
!              and nSyBlk                                              *
!     vInt   : Buffer of 2el integral values                           *
!     nSqNum : sequence number of the integral relative to             *
!              the first adress of the symmetry block                  *
!     nSyBlk : symmetry block number of an integral                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use TwoDat, only: nBatch, RAMD
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nUt
real(kind=wp), intent(in) :: vInt(nUt), nSqNum(nUt), nSyBlk(nUt)
integer(kind=iwp) :: iBatch, iOff, iSyBlk, iUt

!----------------------------------------------------------------------*
!     scatter the 2el integrals to the appropriate position            *
!     of the virtual disk                                              *
!----------------------------------------------------------------------*

!write(u6,'(2X,5(I4,I8,F12.8))') (nSyBlk(iUt),RAMD%adr(nBatch(nSyBlk(iUt)))+nSqNum(iUt)-1,vInt(iUt),iUt=1,nUt)
do iUt=1,nUt
  iSyBlk = int(nSyBlk(iUt),kind=iwp)
  iBatch = nBatch(iSyBlk)
  iOff = RAMD%adr(iBatch)+int(nSqNum(iUt),kind=iwp)
  RAMD%ints(iOff) = vInt(iUt)
end do

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SORT1C
