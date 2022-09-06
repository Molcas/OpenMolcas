!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_Reorder_RI(Vec,lVec,nVec,iSym)

use ChoArr, only: iRS2F

implicit real*8(a-h,o-z)
real*8 Vec(lVec,nVec)
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"
real*8, allocatable :: Scr(:)
integer, allocatable :: iF2RS(:)
! Statement functions
MulD2h(i,j) = ieor(i-1,j-1)+1
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

if (nVec < 1) return
if (lVec < 1) return
if ((lVec /= nnBstR(iSym,1)) .or. (nVec > NumCho(iSym))) then
  call SysAbendMsg('Cho_Reorder_RI','Input argument error!',' ')
end if
if (nnShl /= nnShl_Tot) then
  call SysAbendMsg('Cho_Reorder_RI','Screening is not allowed!','(nnShl /= nnShl_Tot)')
end if

! Set mapping from global address to reduced set.
! -----------------------------------------------

liF2RS = nBasT*(nBasT+1)/2
call mma_allocate(iF2RS,liF2RS,Label='iF2RS')
iF2RS(:) = 0
do iRS=1,nnBstR(iSym,1)
  iRS_tot = iiBstR(iSym,1)+iRS
  na = iRS2F(1,iRS_tot)
  nb = iRS2F(2,iRS_tot)
  nab = iTri(na,nb)
  iF2RS(nab) = iRS
end do

! Reorder.
! --------

lScr = lVec
call mma_allocate(Scr,lScr,Label='Scr')
do iVec=1,nVec

  Scr(:) = Vec(:,iVec)
  kFrom = 0
  do iSymb=1,nSym

    iSyma = MulD2h(iSymb,iSym)

    if (iSyma > iSymb) then
      do ib=1,nBas(iSymb)
        nb = iBas(iSymb)+ib
        do ia=1,nBas(iSyma)
          na = iBas(iSyma)+ia
          nab = iTri(na,nb)
          iRS = iF2RS(nab)
#         ifdef _DEBUGPRINT_
          if ((iRS < 1) .or. (iRS > nnBstR(iSym,1))) then
            call SysAbendMsg('Cho_Reorder_RI','Index out of bounds',' ')
          end if
#         endif
          kFrom = kFrom+1
          Vec(iRS,iVec) = Scr(kFrom)
        end do
      end do
    else if (iSyma == iSymb) then
      do ia=1,nBas(iSyma)
        na = iBas(iSyma)+ia
        do ib=1,ia
          nb = iBas(iSymb)+ib
          nab = iTri(na,nb)
          iRS = iF2RS(nab)
#         ifdef _DEBUGPRINT_
          if ((iRS < 1) .or. (iRS > nnBstR(iSym,1))) then
            call SysAbendMsg('Cho_Reorder_RI','Index out of bounds',' ')
          end if
#         endif
          kFrom = kFrom+1
          Vec(iRS,iVec) = Scr(kFrom)
        end do
      end do
    end if

  end do

end do
call mma_deallocate(Scr)
call mma_deallocate(iF2RS)

end subroutine Cho_Reorder_RI
