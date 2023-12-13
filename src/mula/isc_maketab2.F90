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
subroutine ISC_MakeTab2(m_max,maxOrd,mSiz,mMat,Graph1,Graph2,nOsc)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: m_max, mSiz, nOsc
integer(kind=iwp), intent(out) :: maxOrd, mMat(0:mSiz,nOsc), Graph1(m_max+1,nOsc+1), Graph2(m_max+1,m_max+1,nOsc)
integer(kind=iwp) :: i, iDet, iDNR, iOsc, iQ, iQ1, iQ2, iQuanta, m, n, nd, nQuanta, nTabDim, nvTabDim
integer(kind=iwp), allocatable :: iVec(:), Num(:)
integer(kind=iwp), external :: iDetnr

!write(u6,*)
!write(u6,*) 'CGGt[ISC_MakeTab2] Infos:'
!write(u6,*) 'Graph1(',m_max+1,',',nOsc+1,')'
!write(u6,*) 'Graph2(',m_max+1,',',m_max+1,',',nOsc,')'
!write(u6,*) 'mInc,mDec,mMat(0:',msiz,',',nosc,')'
!write(u6,*) '-------------------------------------------'

! Initialize.

!write(u6,*)
!write(u6,*) '                 msiz,nosc==',msiz,nosc
!mInc(:,:) = 0
!mDec(:,:) = 0
mMat(:,:) = 0
if (m_max == 0) return
!write(u6,*) '                 Calling TabDim'
call TabDim(m_max,nOsc,nTabDim)
maxOrd = nTabDim-1

! Set up the vertex table
!write(u6,*) '                 Graph1,dim=',(m_max+1)*(nOsc+1)
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
!write(u6,*) '                 Num,dim=',m_max
call mma_allocate(Num,[0,m_max],label='Number')
Num(0) = 0
N = 0
do m=1,m_max
  N = N+Graph1(m,nosc+1)
  Num(m) = N
end do
!write(u6,*) '                 Graph2,dim=',nOsc*((m_max+1)**2)
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

!write(u6,*) 'CGGt[MakeTab2] Vec,dim=',nOsc
call mma_allocate(iVec,nOsc,label='iVec')
do iQuanta=1,m_max
  !write(u6,*) 'CGGt[MakeTab2] iQuanta=',iQuanta
  iVec(:) = 0
  iQ = -1
  iVec(1) = -1

  !write(u6,*) 'CGGt[MakeTab2] Call TabDim - 1'
  !write(u6,*) '    iQuanta,nOsc,nd==',iQuanta,nOsc,nd
  call TabDim(iQuanta,nOsc,nd)
  !write(u6,*) 'CGGt[MakeTab2] Call TabDim - 2'
  !write(u6,*) '     iQuanta-1,nOsc,nvTabDim==',iQuanta-1,nOsc,nvTabDim
  call TabDim(iQuanta-1,nOsc,nvTabDim)

  nd = nd-nvTabDim

  do iDet=1,nD
    iVec(1) = iVec(1)+1
    iQ = iQ+1
    if (iQ > iQuanta) then
      do i=1,nOsc-1
        if (iQ <= iQuanta) exit
        iQ = iQ-iVec(i)+1
        iVec(i) = 0
        iVec(i+1) = iVec(i+1)+1
      end do
    end if
    iVec(nOsc) = iQuanta-iq
    iDNR = iDetnr(iVec,Graph2,nosc,m_max)
    mMat(iDNR,:) = iVec(:)
  end do
end do

call mma_deallocate(iVec)

return

end subroutine ISC_MakeTab2
