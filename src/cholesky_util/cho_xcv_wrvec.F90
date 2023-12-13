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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_XCV_WrVec(irc,Vec,l_Vec,NVT,l_NVT,myRankSP,l_myRankSP,SP)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: write partial Cholesky vectors to disk.
!          (Parallel two-step algorithm)

use Cholesky, only: Cho_Real_Par
#ifdef _DEBUGPRINT_
use Cholesky, only: nnBstRSh, nnShl, nSym, NumCho
#endif
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_Vec, l_NVT, NVT(l_NVT), l_myRankSP, myRankSP(l_myRankSP), SP
real(kind=wp), intent(_IN_) :: Vec(l_Vec)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iSP, iSym, n

if ((l_NVT < nSym) .or. (l_myRankSP < nnShl)) then
  irc = -1
  return
end if
if ((SP < 1) .or. (SP > nnShl)) then
  irc = -2
  return
end if
iSP = myRankSP(SP)
n = nnBstRSh(1,iSP,2)*NVT(1)
do iSym=2,nSym
  n = n+nnBstRSh(iSym,iSP,2)*NVT(iSym)
end do
if (l_Vec < n) then
  irc = -3
  return
end if

if (.not. Cho_Real_Par) then
  n = 0
  do iSym=1,nSym
    if (NVT(iSym) /= NumCho(iSym)) n = n+1
  end do
  if (n /= 0) then
    irc = -4
    return
  end if
  if (SP /= myRankSP(SP)) then
    irc = -5
    return
  end if
end if
#endif

if (Cho_Real_Par) then
  ! Parallel: block write to temp files
  call Cho_XCV_WrVec_Par(irc,Vec,NVT,myRankSP,SP)
else
  ! Serial: write directly to vector files
  call Cho_XCV_WrVec_Ser(irc,Vec,SP)
end if

end subroutine Cho_XCV_WrVec
