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
! Copyright (C) 1995,1998, Niclas Forsberg                             *
!               1995,1998, Anders Bernhardsson                         *
!***********************************************************************

!module TabMod

!  Contains:
!    MakeTab   (m_max,maxOrd,maxIncOrd,mMat,mInc,mDec)
!    TabDim    (nDim,nOsc) Result(nTabDim)
!    iDetNr    (iocc,graph,nosc,m_max)  Result(iDetNr)
!
!  Written by:
!    Niclas Forsberg & Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1995&1998.

!vv private

!contains

subroutine MakeTab2(m_max,maxOrd,maxIncOrd,msiz,mMat,mInc,mDec,nOsc)
!  Purpose:
!    Create tables used in FCval.
!
!  Input:
!    nOsc      : Integer - the the number of oscillators.
!    m_max     : Integer - the maximum sum of the quantum numbers.
!    maxOrd    : Integer - number of rows in mMat.
!    nTabDim   : Integer
!    Osc_Shift : Integer array
!
!  Output:
!    mMat      : Two dimensional integer array
!    mInc      : Two dimensional integer array
!    mDec      : Two dimensional integer array
!
!  Calls:
!    none
!
!  Written by:
!    Niclas Forsberg Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1995&1998

implicit real*8(a-h,o-z)
integer mInc(0:msiz,nosc), mDec(0:msiz,nosc), mmat(0:msiz,nosc)
!integer iocc(10)  ! test
!integer nTabDim, nvTabDim
#include "WrkSpc.fh"

call GetMem('Graph1','Allo','Inte',ipGraph1,(m_max+1)*(nOsc+1))
call GetMem('Graph2','Allo','Inte',ipGraph2,(m_max+1)*(m_max+1)*nOsc)

call MakeTab2_a(m_max,maxOrd,maxIncOrd,msiz,mMat,mInc,mDec,nOsc,iWork(ipgraph1),iWork(ipgraph2))

call GetMem('Graph1','Free','Inte',ipGraph1,(m_max+1)*(nOsc+1))
call GetMem('Graph2','Free','Inte',ipGraph2,(m_max+1)*(m_max+1)*nOsc)

end subroutine MakeTab2
!####
subroutine MakeTab2_a(m_max,maxOrd,maxIncOrd,msiz,mMat,mInc,mDec,nOsc,graph1,graph2)
!  Purpose:
!    Create tables used in FCval.
!
!  Input:
!    nOsc      : Integer - the the number of oscillators.
!    m_max     : Integer - the maximum sum of the quantum numbers.
!    maxOrd    : Integer - number of rows in mMat.
!    nTabDim   : Integer
!    Osc_Shift : Integer array
!
!  Output:
!    mMat      : Two dimensional integer array
!    mInc      : Two dimensional integer array
!    mDec      : Two dimensional integer array
!
!  Calls:
!    none
!
!  Written by:
!    Niclas Forsberg Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1995&1998

implicit real*8(a-h,o-z)
integer mInc(0:msiz,nosc), mDec(0:msiz,nosc), mmat(0:msiz,nosc)
integer Graph1(m_max+1,nOsc+1)
integer Graph2(m_max+1,m_max+1,nOsc)
!integer iocc(10)  ! test
integer nTabDim, nvTabDim
#include "WrkSpc.fh"

! Initialize.

do iv=0,msiz
  do jv=1,nosc
    mInc(iv,jv) = 0
    mDec(iv,jv) = 0
    mMat(iv,jv) = 0
  end do
end do
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

call GetMem('ivec','Allo','INTE',ipivec,nOsc)
do iQuanta=1,m_max
  do iv=1,nOsc
    iWork(ipiVec+iv-1) = 0
  end do
  iQ = -1
  iWork(ipiVec) = -1

  call TabDim2_drv(iQuanta,nOsc,nd)
  call TabDim2_drv(iQuanta-1,nOsc,nvTabDim)

  nd = nd-nvTabDim

  do iDet=1,nD
    iWork(ipiVec) = iWork(ipiVec)+1
    iQ = iQ+1
    if (iQ > iQuanta) then
      do i=1,nOsc-1
        if (iQ <= iQuanta) exit
        iQ = iQ-iWork(ipiVec+i-1)+1
        iWork(ipiVec+i-1) = 0
        iWork(ipiVec+i) = iWork(ipiVec+i)+1
      end do
    end if
    iWork(ipiVec+nOsc-1) = iQuanta-iq
    iDNR = iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
    do iv=1,nOsc
      mMat(iDnr,iv) = iWork(ipiVec+iv-1)
    end do
  end do
end do

! Create mInc.
!minc = -1
do iv=0,msiz
  do jv=1,nosc
    mInc(iv,jv) = -1
  end do
end do

call TabDim2_drv(m_max-1,nosc,nvTabDim)
maxIncOrd = nvTabDim-1
do i=0,maxIncOrd
  do iv=1,nOsc
    iWork(ipivec+iv-1) = mMat(i,iv)
  end do
  do j=1,nOsc
    iWork(ipivec+j-1) = iWork(ipivec+j-1)+1
    mInc(i,j) = iDetnr(iWork(ipivec),Graph2,nosc,m_max)
    iWork(ipivec+j-1) = iWork(ipivec+j-1)-1
  end do
end do

! Create mDec.

do iv=1,nOsc
  mdec(0,iv) = -1
end do
do i=1,maxOrd
  do j=1,nOsc
    if (mmat(i,j) /= 0) then
      do iv=1,nOsc
        iWork(ipivec+iv-1) = mMat(i,iv)
      end do
      iWork(ipivec+j-1) = iWork(ipivec+j-1)-1
      mDec(i,j) = iDetnr(iWork(ipivec),Graph2,nosc,m_max)
      do iv=1,nOsc
        iWork(ipivec+iv-1) = iWork(ipivec+j-1)+1
      end do
    else
      mdec(i,j) = -1
    end if
  end do
end do

call GetMem('ivec','Free','INTE',ipivec,nOsc)

end subroutine MakeTab2_a
