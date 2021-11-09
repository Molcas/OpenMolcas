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
!               2009, Giovanni Ghigo                                   *
!***********************************************************************

!  Written by:
!    Niclas Forsberg & Anders Bernhardsson,
!    Dept. of Theoretical Chemistry, Lund University, 1995&1998.
!  Modified by:
!    Giovanni Ghigo for InterSystem Crossing
!    Dept. of Chimica Generale e Chimica Organica, University of Torino, 2009 June.
!
!  Purpose:
!    Create tables used in FCval.
!
!  Contains:
!    MakeTab   (m_max,maxOrd,maxIncOrd,mMat,mInc,mDec)
!    TabDim    (nDim,nOsc) Result(nTabDim)
!    iDetNr    (iocc,graph,nosc,m_max)  Result(iDetNr)
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
subroutine ISC_MakeTab2(m_max,maxOrd,maxIncOrd,mSiz,mMat,Graph1,Graph2,nOsc)

implicit real*8(a-h,o-z)
integer mMat(0:msiz,nosc)
integer Graph1(m_max+1,nOsc+1)
integer Graph2(m_max+1,m_max+1,nOsc)
integer nTabDim, nvTabDim
#include "WrkSpc.fh"

!write(u6,*)
!write(u6,*) 'CGGt[ISC_MakeTab2_a] Infos:                   '
!write(u6,*) 'Graph1(',m_max+1,',',nOsc+1,')'
!write(u6,*) 'Graph2(',m_max+1,',',m_max+1,',',nOsc,')'
!write(u6,*) 'mInc,mDec,mMat(0:',msiz,',',nosc,')'
!write(u6,*) '-------------------------------------------   '

! Initialize.

!write(u6,*)
!write(u6,*) '                 msiz,nosc==',msiz,nosc
!call GetMem('Test_1','LIST','INTE',iDum,iDum)
do iv=0,msiz
  do jv=1,nosc
    !mInc(iv,jv) = 0
    !mDec(iv,jv) = 0
    mMat(iv,jv) = 0
  end do
end do
if (m_max == 0) return
!write(u6,*) '                 Calling TabDim_drv'
call TabDim_drv(m_max,nOsc,nTabDim)
maxOrd = nTabDim-1

! Set up the vertex table
!write(u6,*) '                 Graph1,dim=',(m_max+1)*(nOsc+1)
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
!write(u6,*) '                 ipNumber,dim=',(m_max)
call GetMem('Number','Allo','INTE',ipNumber,m_max+1)
do iv=0,m_max
  iWork(ipNumber+iv) = 0
end do
N = 0
do m=1,m_max
  N = N+Graph1(m,nosc+1)
  iWork(ipNumber+m) = n
end do
!write(u6,*) '                 Graph2,dim=',nOsc*((m_max+1)**2)
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

!write(u6,*) 'CGGt[MakeTab2_a] Vec,dim=',nOsc
call GetMem('ivec','Allo','INTE',ipiVec,nOsc)
do iQuanta=1,m_max
  !write(u6,*) 'CGGt[MakeTab2_a] iQuanta=',iQuanta
  do iv=1,nOsc
    iWork(ipiVec+iv-1) = 0
  end do
  iQ = -1
  iWork(ipiVec) = -1

  !write(u6,*) 'CGGt[MakeTab2_a] Call TabDim2_drv - 1 '
  !write(u6,*) '    iQuanta,nOsc,nd==',iQuanta,nOsc,nd
  call TabDim2_drv(iQuanta,nOsc,nd)
  !write(u6,*) 'CGGt[MakeTab2_a] Call TabDim2_drv - 2 '
  !write(u6,*) '     iQuanta-1,nOsc,nvTabDim==',iQuanta-1,nOsc,nvTabDim
  call TabDim2_drv(iQuanta-1,nOsc,nvTabDim)

  nd = nd-nvTabDim

  do iDet=1,nD
    iWork(ipiVec) = iWork(ipiVec)+1
    iQ = iQ+1
    if (iQ > iQuanta) then
      do i=1,nOsc-1
        if (iQ <= iQuanta) goto 99
        iQ = iQ-iWork(ipiVec+i-1)+1
        iWork(ipiVec+i-1) = 0
        iWork(ipiVec+i) = iWork(ipiVec+i)+1
      end do
    end if
99  continue
    iWork(ipiVec+nOsc-1) = iQuanta-iq
    iDNR = iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
    do iv=1,nOsc
      mMat(iDNR,iv) = iWork(ipiVec+iv-1)
    end do
  end do
end do

call GetMem('iVec','Free','INTE',ipiVec,nOsc)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(MaxIncOrd)

end subroutine ISC_MakeTab2
!####
subroutine Mk_nIncDec(m_max,nOrd,msiz,mInc,mDec,mMat,Graph2,nOsc)

implicit real*8(a-h,o-z)
integer mInc(0:msiz,nosc), mDec(0:msiz,nosc)
integer mMat(0:nOrd,nosc)
integer Graph2(m_max+1,m_max+1,nOsc)
integer m_max, nOrd, msiz
#include "WrkSpc.fh"

call GetMem('iVec','Allo','INTE',ipiVec,nOsc)

! Create mInc.

do iv=0,msiz
  do jv=1,nosc
    mInc(iv,jv) = -1
  end do
end do
do i=0,msiz
  do iv=1,nOsc
    iWork(ipiVec+iv-1) = mMat(i,iv)
  end do
  do j=1,nOsc
    iWork(ipiVec+j-1) = iWork(ipiVec+j-1)+1
    mInc(i,j) = iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
    iWork(ipiVec+j-1) = iWork(ipiVec+j-1)-1
  end do
end do

! Create mDec.

do iv=1,nOsc
  mDec(0,iv) = -1
end do
do i=1,msiz
  do j=1,nOsc
    if (mMat(i,j) /= 0) then
      do iv=1,nOsc
        iWork(ipiVec+iv-1) = mMat(i,iv)
      end do
      iWork(ipiVec+j-1) = iWork(ipiVec+j-1)-1
      mDec(i,j) = iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
      do iv=1,nOsc
        iWork(ipiVec+iv-1) = iWork(ipiVec+j-1)+1
      end do
    else
      mDec(i,j) = -1
    end if
  end do
end do

call GetMem('ivec','Free','INTE',ipiVec,nOsc)

return

end subroutine Mk_nIncDec
