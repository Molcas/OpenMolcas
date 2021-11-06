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

subroutine ISCD_MakeGraphs(m_max,maxOrd,maxIncOrd,Graph1,Graph2,nOsc)

implicit real*8(a-h,o-z)
integer Graph1(m_max+1,nOsc+1)
integer Graph2(m_max+1,m_max+1,nOsc)
integer nTabDim
#include "WrkSpc.fh"

! Initialize.
if (m_max == 0) return
call TabDim_drv(m_max,nOsc,nTabDim)
maxOrd = nTabDim-1

! Set up the vertex table
do iv=1,m_max+1
  do jv=1,nOsc+1
    Graph1(iv,jv) = 0
  end do
end do
do iv=1,m_max+1
  Graph1(iv,2) = 1
end do
do jv=1,nOsc+1
  Graph1(1,jv) = 1
end do
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
call GetMem('Number','Allo','INTE',ipNumber,m_max+1)
do iv=0,m_max
  iWork(ipNumber+iv) = 0
end do
N = 0
do m=1,m_max
  N = N+Graph1(m,nosc+1)
  iWork(ipNumber+m) = n
end do
do iv=1,m_max+1
  do jv=1,m_max+1
    do kv=1,nOsc
      Graph2(iv,jv,kv) = 0
    end do
  end do
end do
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
    Graph2(iQ1+1,iQ2+1,nOsc) = Graph2(iQ1+1,iQ2+1,nOsc)+iWork(ipNumber+iQ1)
  end do
end do

call GetMem('Number','Free','INTE',ipNumber,m_max+1)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(maxIncOrd)

end subroutine ISCD_MakeGraphs
