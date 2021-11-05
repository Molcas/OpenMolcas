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
!!-----------------------------------------------------------------------!
!!
!       Module TabMod
!!
!!  Contains:
!!    MakeTab   (m_max,maxOrd,maxIncOrd,mMat,mInc,mDec)
!!    TabDim    (nDim,nOsc) Result(nTabDim)
!!    iDetNr    (iocc,graph,nosc,m_max)  Result(iDetNr)
!!
!!  Written by:
!!    Niclas Forsberg & Anders Bernhardsson,
!!    Dept. of Theoretical Chemistry, Lund University, 1995&1998.
!!
!!-----------------------------------------------------------------------!
!!
!vv       Private

!       Contains

!!-----------------------------------------------------------------------!

!!-----------------------------------------------------------------------!
!!
      Subroutine MakeTab2(                                              &
     &      m_max,maxOrd,maxIncOrd,msiz,mMat,mInc,mDec,nOsc)
!!
!!  Purpose:
!!    Create tables used in FCval.
!!
!!  Input:
!!    nOsc      : Integer - the the number of oscillators.
!!    m_max     : Integer - the maximum sum of the quantum numbers.
!!    maxOrd    : Integer - number of rows in mMat.
!!    nTabDim   : Integer
!!    Osc_Shift : Integer array
!!
!!  Output:
!!    mMat      : Two dimensional integer array
!!    mInc      : Two dimensional integer array
!!    mDec      : Two dimensional integer array
!!
!!  Calls:
!!    none
!!
!!  Written by:
!!    Niclas Forsberg Anders Bernhardsson,
!!    Dept. of Theoretical Chemistry, Lund University, 1995&1998
!!
      Implicit Real*8 ( a-h,o-z )
      Integer mInc (0:msiz,nosc),mDec (0:msiz,nosc),                    &
     &  mmat (0:msiz,nosc) !  (0:msiz,nosc)
!       Integer iocc(10)  ! test
!       Integer  nTabDim,nvTabDim
#include "WrkSpc.fh"

      Call GetMem('Graph1','Allo','Inte',                               &
     &  ipGraph1,(m_max+1)*(nOsc+1))
      Call GetMem('Graph2','Allo','Inte',                               &
     &  ipGraph2,(m_max+1)*(m_max+1)*nOsc)

      Call MakeTab2_a(                                                  &
     &      m_max,maxOrd,maxIncOrd,msiz,mMat,mInc,mDec,                 &
     & nOsc,iWork(ipgraph1),iWork(ipgraph2))

      Call GetMem('Graph1','Free','Inte',                               &
     &  ipGraph1,(m_max+1)*(nOsc+1))
      Call GetMem('Graph2','Free','Inte',                               &
     &  ipGraph2,(m_max+1)*(m_max+1)*nOsc)

      End


!!-----------------------------------------------------------------------!
!!
      Subroutine MakeTab2_a(                                            &
     &      m_max,maxOrd,maxIncOrd,msiz,mMat,mInc,mDec,                 &
     &  nOsc,graph1,graph2)
!!
!!  Purpose:
!!    Create tables used in FCval.
!!
!!  Input:
!!    nOsc      : Integer - the the number of oscillators.
!!    m_max     : Integer - the maximum sum of the quantum numbers.
!!    maxOrd    : Integer - number of rows in mMat.
!!    nTabDim   : Integer
!!    Osc_Shift : Integer array
!!
!!  Output:
!!    mMat      : Two dimensional integer array
!!    mInc      : Two dimensional integer array
!!    mDec      : Two dimensional integer array
!!
!!  Calls:
!!    none
!!
!!  Written by:
!!    Niclas Forsberg Anders Bernhardsson,
!!    Dept. of Theoretical Chemistry, Lund University, 1995&1998
!!
      Implicit Real*8 ( a-h,o-z )
      Integer mInc (0:msiz,nosc),mDec (0:msiz,nosc),                    &
     &  mmat (0:msiz,nosc) !  (0:msiz,nosc)
      Integer Graph1 (m_max+1,nOsc+1)
      Integer Graph2 (m_max+1,m_max+1,nOsc)
!       Integer iocc(10)  ! test
      Integer  nTabDim,nvTabDim
#include "WrkSpc.fh"
!!
!!---- Initialize.

      do iv=0,msiz
      do jv=1,nosc
      mInc(iv,jv) = 0
      mDec(iv,jv) = 0
      mMat(iv,jv) = 0
      enddo
      enddo
      If ( m_max.eq.0 ) Return
      Call TabDim_drv(m_max,nOsc,nTabDim)
      maxOrd = nTabDim-1
!!
!!---- Set up the vertex table
      do iv=1,m_max+1
      do jv=1,nOsc+1
      Graph1(iv,jv) = 0
      enddo
      enddo
      do iv=1,m_max+1
      Graph1(iv,2) = 1
      enddo
      do jv=1,nOsc+1
      Graph1(1,jv) = 1
      enddo
      If ( nOsc.gt.1 ) Then
      Do iOsc = 2,nOsc
      n = 0
      Do nQuanta = 0,m_max
      n = n+Graph1(nQuanta+1,iOsc)
      Graph1(nQuanta+1,iOsc+1) = n
      End Do
      End Do
      End If
!!
!!---- set up the arc table
      Call GetMem('Number','Allo','INTE',ipNumber,m_max+1)
      do iv=0,m_max
      iWork(ipNumber+iv) = 0
      enddo
      N = 0
      Do m = 1,m_max
      N = N+Graph1(m,nosc+1)
      iWork(ipNumber+m) = n
      End Do
      do iv=1,m_max+1
      do jv=1,m_max+1
      do kv=1,nOsc
      Graph2(iv,jv,kv) = 0
      enddo
      enddo
      enddo
      Do iOsc = 1,nosc
      Do iQ1 = 0,m_max         ! Where we are going
      Do iQ2 = 0,iQ1-1      ! Where we came from
      Do i = iQ2+1,iq1   ! Sum over preceding paths
      Graph2(iQ1+1,iQ2+1,iOsc) = Graph1(i+1,iOsc)+                      &
     &             Graph2(iQ1+1,iQ2+1,iOsc)
      End Do
      End Do
      End Do
      End Do
!!
      Do iQ1 = 0,m_max            ! Where we are going
      Do iQ2 = 0,iq1           ! Where we came from
      Graph2(iQ1+1,iQ2+1,nOsc) = Graph2(iQ1+1,iQ2+1,nOsc)+              &
     &       iWork(ipNumber+iQ1)
      End Do
      End Do
!!
      Call GetMem('Number','Free','INTE',ipNumber,m_max+1)
!!
!!
      Call GetMem('ivec','Allo','INTE',ipivec,nOsc)
      Do iQuanta=1,m_max
        do iv=1,nOsc
          iWork(ipiVec+iv-1)=0
        enddo
      iQ=-1
      iWork(ipiVec)=-1

      Call TabDim2_drv(iQuanta,nOsc,nd)
      Call TabDim2_drv(iQuanta-1,nOsc,nvTabDim)

      nd=nd-nvTabDim

      Do iDet=1,nD
      iWork(ipiVec)=iWork(ipiVec)+1
      iQ=iQ+1
      If (iQ.gt.iQuanta) Then
       Do i=1,nOsc-1
        if(iQ.le.iQuanta) goto 99
        iQ=iQ-iWork(ipiVec+i-1)+1
        iWork(ipiVec+i-1)=0
        iWork(ipiVec+i)=iWork(ipiVec+i)+1
       End Do
      End If
  99  Continue
      iWork(ipiVec+nOsc-1)=iQuanta-iq
      iDNR=iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
      do iv=1,nOsc
      mMat(iDnr,iv)=iWork(ipiVec+iv-1)
      enddo
      End Do
      End Do
!!
!!---- Create mInc.
!       minc=-1
      do iv=0,msiz
      do jv=1,nosc
      mInc(iv,jv) = -1
      enddo
      enddo

      Call TabDim2_drv(m_max-1,nosc,nvTabDim)
      maxIncOrd = nvTabDim-1
      Do i = 0,maxIncOrd
      do iv=1,nOsc
      iWork(ipivec+iv-1) = mMat(i,iv)
      enddo
      Do j = 1,nOsc
      iWork(ipivec+j-1) = iWork(ipivec+j-1)+1
      mInc(i,j)=iDetnr(iWork(ipivec),Graph2,nosc,m_max)
      iWork(ipivec+j-1) = iWork(ipivec+j-1)-1
      End Do
      End Do
!!
!!---- Create mDec.

      do iv=1,nOsc
      mdec(0,iv)=-1
      enddo
      Do i = 1,maxOrd
      Do j = 1,nOsc
      If (mmat(i,j).ne.0)Then
      do iv=1,nOsc
      iWork(ipivec+iv-1) = mMat(i,iv)
      enddo
      iWork(ipivec+j-1)=iWork(ipivec+j-1)-1
      mDec(i,j)=iDetnr(iWork(ipivec),Graph2,nosc,m_max)
      do iv=1,nOsc
      iWork(ipivec+iv-1) = iWork(ipivec+j-1)+1
      enddo
      Else
      mdec(i,j)=-1
      End IF
      End Do
      End Do
!!

      Call GetMem('ivec','Free','INTE',ipivec,nOsc)
      End
