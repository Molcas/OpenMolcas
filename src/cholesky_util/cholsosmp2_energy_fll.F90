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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoLSOSMP2_Energy_Fll(N,w,t,EOcc,EVir,Delete,EMP2,irc)
!
! Thomas Bondo Pedersen, December 2012.
!
! Compute Laplace-SOS-MP2 energy correction from full Cholesky
! vectors (i.e., not batched).

use Cholesky, only: NumCho, nSym
use ChoMP2, only: DecoMP2, Laplace_BlockSize, nMP2Vec, nT1am
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: w(N), t(N), EOcc(*), EVir(*)
logical(kind=iwp), intent(in) :: Delete
real(kind=wp), intent(out) :: EMP2
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: iSym, l, l_V, l_X, Nai, need, nEnrVec(8)
character(len=*), parameter :: SecNam = 'ChoLSOSMP2_Energy_Fll'

! check if there is enough memory to read through vector files
! only once
if (DecoMP2) then
  nEnrVec(1:nSym) = nMP2Vec(1:nSym)
else
  nEnrVec(1:nSym) = NumCho(1:nSym)
end if
l_X = 0
l_V = 0
do iSym=1,nSym
  Nai = nT1am(iSym)
  if ((Nai > 0) .and. (nEnrVec(iSym) > 0)) then
    l_X = max(l_X,min(Laplace_BlockSize,nEnrVec(iSym)))
    l_V = max(l_V,Nai*nEnrVec(iSym))
  end if
end do
need = l_X+2*l_V
call mma_maxDBLE(l)
l = l-100
if ((l < 1) .or. (need >= l)) then ! not enough memory for one read
  call ChoLSOSMP2_Energy_Fll2(N,w,t,EOcc,EVir,Delete,EMP2,irc)
  if (irc /= 0) write(u6,'(A,A,I10)') SecNam,': Cho_LSOSMP2_Energy_Fll2 returned',irc
else ! enough memory for one read through vector files
  call ChoLSOSMP2_Energy_Fll1(N,w,t,EOcc,EVir,Delete,EMP2,irc)
  if (irc /= 0) write(u6,'(A,A,I10)') SecNam,': Cho_LSOSMP2_Energy_Fll1 returned',irc
end if

end subroutine ChoLSOSMP2_Energy_Fll
