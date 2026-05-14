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

subroutine Cho_XCV_RdVec(irc,Vec,l_Vec,NVT,myRankSP,n_myRankSP,J1,J2,iSym)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: Read partial Cholesky vectors J1 to J2 on disk.
!          (Parallel two-step algorithm)

#ifdef _DEBUGPRINT_
use Cholesky, only: Cho_Real_Par, nnBstRSh, nnShl, nSym, NumCho
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_Vec, NVT, n_myRankSP, myRankSP(n_myRankSP), J1, J2, iSym
real(kind=wp), intent(out) :: Vec(l_Vec)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, n
#endif

irc = 0
if ((n_myRankSP == 0) .or. (J2-J1+1 == 0)) return ! nothing to do

#ifdef _DEBUGPRINT_
if ((n_myRankSP < 1) .or. (n_myRankSP > nnShl) .or. (iSym < 1) .or. (iSym > nSym) .or. (NVT < 1)) then
  irc = -1
  return
end if
if ((J1 < 1) .or. (J1 > NVT) .or. (J2 < 1) .or. (J2 > NVT) .or. (J1 > J2) .or. (J2-J1+1 > NVT)) then
  irc = -2
  return
end if
n = nnBstRSh(iSym,myRankSP(1),2)*(J2-J1+1)
do i=2,n_myRankSP
  n = n+nnBstRSh(iSym,myRankSP(i),2)*(J2-J1+1)
end do
if (l_Vec < n) then
  irc = -3
  return
end if

if (.not. Cho_Real_Par) then
  if (NVT /= NumCho(iSym)) then
    irc = -4
    return
  end if
  n = 0
  do i=1,n_myRankSP
    if (myRankSP(i) /= i) n = n+1
  end do
  if (n /= 0) then
    irc = -5
    return
  end if
end if
#endif

! Block read on temp files
call Cho_XCV_RdVec_(irc,Vec,myRankSP,n_myRankSP,NVT,J1,J2,iSym)

end subroutine Cho_XCV_RdVec
