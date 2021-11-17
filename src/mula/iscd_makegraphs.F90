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

subroutine ISCD_MakeGraphs(m_max,maxOrd,Graph1,Graph2,nOsc)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: m_max, nOsc
integer(kind=iwp), intent(out) :: maxOrd, Graph1(m_max+1,nOsc+1), Graph2(m_max+1,m_max+1,nOsc)
integer(kind=iwp) :: i, iOsc, iQ1, iQ2, m, n, nQuanta, nTabDim
integer(kind=iwp), allocatable :: Num(:)

! Initialize.
if (m_max == 0) return
call TabDim(m_max,nOsc,nTabDim)
maxOrd = nTabDim-1

! Set up the vertex table
Graph1(:,:) = 0
Graph1(:,2) = 1
Graph1(1,:) = 1
if (nOsc > 1) then
  do iOsc=2,nOsc
    n = 0
    do nQuanta=0,m_max
      n = n+Graph1(nQuanta+1,iOsc)
      Graph1(nQuanta+1,iOsc+1) = n
    end do
  end do
end if

! set up the arc table
call mma_allocate(Num,[0,m_max],label='Number')
Num(0) = 0
N = 0
do m=1,m_max
  N = N+Graph1(m,nosc+1)
  Num(m) = N
end do
Graph2(:,:,:) = 0
do iOsc=1,nosc
  do iQ1=0,m_max      ! Where we are going
    do iQ2=0,iQ1-1    ! Where we came from
      do i=iQ2+1,iq1  ! Sum over preceding paths
        Graph2(iQ1+1,iQ2+1,iOsc) = Graph1(i+1,iOsc)+Graph2(iQ1+1,iQ2+1,iOsc)
      end do
    end do
  end do
end do

do iQ1=0,m_max  ! Where we are going
  do iQ2=0,iq1  ! Where we came from
    Graph2(iQ1+1,iQ2+1,nOsc) = Graph2(iQ1+1,iQ2+1,nOsc)+Num(iQ1)
  end do
end do

call mma_deallocate(Num)

return

end subroutine ISCD_MakeGraphs
