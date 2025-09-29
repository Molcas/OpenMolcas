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

subroutine RdDens(d1,nd1,d2,nd2)

use Index_Functions, only: iTri
use MCLR_Data, only: LuJob
use input_mclr, only: iRoot, iTOC, lRoots, nRoots, ntAsh, ntAsh, Weight
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nd1, nd2
real(kind=wp), intent(out) :: D1(nd1), d2(nd2)
real(kind=wp) :: Fact, rdum(1), W
integer(kind=iwp) :: i, iB, iDij, iDkl, iIJKL, j, jB, jDisk, kB, lB
real(kind=wp), allocatable :: D1t(:), D2t(:), G2tt(:)

!                                                                      *
!***********************************************************************
!                                                                      *
d1(:) = Zero
call mma_allocate(G2tt,nd2,Label='G2tt')
call mma_allocate(D2t,nd2,Label='D2t')
call mma_allocate(D1t,nd1,Label='D1t')
G2tt(:) = Zero
jDisk = ITOC(3)
do i=1,lroots
  W = Zero
  do j=1,nroots
    if (iroot(j) == i) W = Weight(j)
  end do
  call dDaFile(LUJOB,2,D1t,nd1,jDisk)
  call dDaFile(LUJOB,0,rdum,nd1,jDisk)
  call dDaFile(LUJOB,2,D2t,ND2,jDisk)
  call dDaFile(LUJOB,0,rdum,ND2,jDisk)
  if (W /= Zero) then
    G2tt(:) = G2tt(:)+w*D2t(:)
    D1(:) = D1(:)+w*D1t(:)
  end if
  !write(u6,*) i,w,'LUJOB',LUJOB
  !call Triprt('D',' ',D1,ntash)
end do
call Put_dArray('D2av',G2tt,nd2)
call Put_dArray('D1av',D1,nd1)

do iB=1,ntash
  do jB=1,iB
    iDij = iTri(ib,jB)
    do kB=1,ib
      do lB=1,kB
        iDkl = iTri(kB,lB)
        fact = One
        if ((iDij >= iDkl) .and. (kB == lB)) fact = Two
        if ((iDij < iDkl) .and. (iB == jB)) fact = Two
        iijkl = iTri(iDij,iDkl)
        D2(iijkl) = Fact*G2tt(iijkl)
      end do
    end do
  end do
end do

call mma_deallocate(G2tt)
call mma_deallocate(D2t)
call mma_deallocate(D1t)

end subroutine RdDens
