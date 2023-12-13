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
! Copyright (C) 1993,1996, Markus P. Fuelscher                         *
!               1993, Per Ake Malmqvist                                *
!***********************************************************************

subroutine SORT1A(nUt,vInt,nSqNum,nSyBlk)
!***********************************************************************
!                                                                      *
!     Purpose: Phase 1 of the bin sorting algorithm                    *
!             Distribute the integrals into the bins where the indices *
!             and the integral values are stored in different buffers. *
!             When a bin is completed append the buffer containing the *
!             integral values to the file LuTwo and the buffer         *
!             containing the indice to LuTmp. Before draining the      *
!             buffers they are sorted in ascending order using the     *
!             ESSL routine ISORTX.                                     *
!                                                                      *
!     Called from: PLF2,INDSFT2                                        *
!                                                                      *
!     Calls to : PKI4,PKR8,SetVec,ISORTX,I4Len,R8Len                   *
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
!     M. P. Fuelscher and P.-AA. Malmqvist                             *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     - modified to use a virtual disk                                 *
!       M. P. Fuelscher, University of Lund, Sweden, 1996              *
!                                                                      *
!***********************************************************************

use TwoDat, only: RAMD
use sort_data, only: lBin, lwIBin, lwVBin, mInt, n_Int
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nUt
real(kind=wp), intent(in) :: vInt(nUt), nSqNum(nUt), nSyBlk(nUt)
#include "print.fh"
integer(kind=iwp) :: iBin, iOpt, iPrint, iRout, iUt, next

!----------------------------------------------------------------------*
!     pick up print level                                              *
!----------------------------------------------------------------------*

iRout = 81
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  write(u6,*) ' >>> Enter SORT1A <<<'
  call dVcPrt('nSqNum',' ',nSqNum,nUt)
  call dVcPrt('nSyBlk',' ',nSyBlk,nUt)
  call dVcPrt('vInt',' ',vInt,nUt)
end if

if (RAMD%act) then
  call Untested('Sort1a (RAMD)')
  call SORT1C(nUt,vInt,nSqNum,nSyBlk)
  return
end if

iOpt = 0 ! Always tight!

!----------------------------------------------------------------------*
!     put the 2el integrals into the bins                              *
!----------------------------------------------------------------------*

do iUt=1,nUt
  iBin = int(nSyBlk(iUt),kind=iwp)
  next = n_Int(iBin)+1
  lwVBin(next,iBin) = vInt(iUt)
  lwIBin(next,iBin) = int(nSqNum(iUt),kind=iwp)
  n_Int(iBin) = next
  mInt(1,iBin) = mInt(1,iBin)+1

  !--------------------------------------------------------------------*
  !   If a bin is completed pack the buffer and put it to disc         *
  !   First, sort the the integrals such that the indices are given in *
  !   strictly ascending order.                                        *
  !--------------------------------------------------------------------*

  if (next >= (lBin-1)) then
    call SaveBin(iBin,iOpt)
  end if
end do

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SORT1A
