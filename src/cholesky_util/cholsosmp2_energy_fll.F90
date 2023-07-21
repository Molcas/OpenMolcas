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

use stdalloc

implicit none
integer N
real*8 w(N)
real*8 t(N)
real*8 EOcc(*)
real*8 EVir(*)
logical Delete
real*8 EMP2
integer irc
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "cholesky.fh"
character(len=21), parameter :: SecNam = 'ChoLSOSMP2_Energy_Fll'
integer nEnrVec(8)
integer l_X
integer l_V
integer iSym
integer Nai
integer need
integer l

! check if there is enough memory to read through vector files
! only once
if (DecoMP2) then
  call iCopy(nSym,nMP2Vec,1,nEnrVec,1)
else
  call iCopy(nSym,NumCho,1,nEnrVec,1)
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
  if (irc /= 0) write(6,'(A,A,I10)') SecNam,': Cho_LSOSMP2_Energy_Fll2 returned',irc
else ! enough memory for one read through vector files
  call ChoLSOSMP2_Energy_Fll1(N,w,t,EOcc,EVir,Delete,EMP2,irc)
  if (irc /= 0) write(6,'(A,A,I10)') SecNam,': Cho_LSOSMP2_Energy_Fll1 returned',irc
end if

end subroutine ChoLSOSMP2_Energy_Fll
