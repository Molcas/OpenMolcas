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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_DiaSP()
!
! Thomas Bondo Pedersen, March 2006.
!
! Purpose: prescreening of diagonal.

use Index_Functions, only: iTri
use Cholesky, only: Cho_PreScreen, iSP2F, nnShl, nnShl_tot, nShell, Thr_PreScreen
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, ij, j
real(kind=wp) :: Tau, Tmax_All
real(kind=wp), allocatable :: TMax(:,:)

if (Cho_PreScreen) then ! prescreening with approx. diagonal

  call mma_allocate(Tmax,nShell,nShell,Label='nShell')

  call Shell_MxSchwz(nShell,Tmax)
  Tmax_All = Tmax(1,1)
  do i=2,nShell
    do j=1,i
      Tmax_All = max(Tmax_All,Tmax(i,j))
    end do
  end do

  Tau = Thr_PreScreen
  nnShl = 0
  do i=1,nShell
    do j=1,i
      if (Tmax_All*Tmax(i,j) > Tau) nnShl = nnShl+1
    end do
  end do
  call mma_allocate(iSP2F,nnShl,Label='iSP2F')

  ij = 0
  do i=1,nShell
    do j=1,i
      if (Tmax_All*Tmax(i,j) > Tau) then
        ij = ij+1
        iSP2F(ij) = iTri(i,j)
      end if
    end do
  end do

  call mma_deallocate(TMax)

else ! no prescreening, include all shell pairs.

  nnShl = nnShl_Tot
  call mma_allocate(iSP2F,nnShl,Label='iSP2F')

  do ij=1,nnShl
    iSP2F(ij) = ij
  end do

end if

#ifdef _DEBUGPRINT_
if (.not. Cho_PreScreen) Tau = Zero
write(LuPri,*) '>>> Exit from Cho_DiaSP:'
write(LuPri,*) '    Screening threshold               : ',Tau
write(LuPri,*) '    Total number of shell pairs       : ',nnShl_Tot
write(LuPri,*) '    Contributing number of shell pairs: ',nnShl
if (nnShl_Tot /= 0) write(LuPri,*) '    Screening-%: ',1.0e2_wp*real(nnShl_Tot-nnShl,kind=wp)/real(nnShl_Tot,kind=wp)
#endif

end subroutine Cho_DiaSP
