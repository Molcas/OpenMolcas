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
!!-----------------------------------------------------------------------!
!!
!!  Written by:
!!    Niclas Forsberg & Anders Bernhardsson,
!!    Dept. of Theoretical Chemistry, Lund University, 1995&1998.
!!  Modified by:
!!    Giovanni Ghigo for InterSystem Crossing
!!    Dept. of Chimica Generale e Chimica Organica, University of Torino, 2009 June.
!!
!!  Purpose:
!!    Create tables used in FCval.
!!
!!  Contains:
!!    MakeTab   (m_max,maxOrd,maxIncOrd,mMat,mInc,mDec)
!!    TabDim    (nDim,nOsc) Result(nTabDim)
!!    iDetNr    (iocc,graph,nosc,m_max)  Result(iDetNr)
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
!!-----------------------------------------------------------------------!
!!
      Subroutine ISC_MakeTab2(m_max,maxOrd,maxIncOrd,mSiz,mMat,         &
     &                        Graph1,Graph2,nOsc)
!!
      Implicit Real*8 ( a-h,o-z )
      Integer mMat (0:msiz,nosc) !  (0:msiz,nosc)
      Integer Graph1 (m_max+1,nOsc+1)
      Integer Graph2 (m_max+1,m_max+1,nOsc)
      Integer  nTabDim,nvTabDim
#include "WrkSpc.fh"
!      Write(6,*)                                                  ! CGGt
!      Write(6,*)'CGGt[ISC_MakeTab2_a] Infos:                   '  ! CGGt
!      Write(6,*)'Graph1(',m_max+1,',',nOsc+1,')'                  ! CGGt
!      Write(6,*)'Graph2(',m_max+1,',',m_max+1,',',nOsc,')'        ! CGGt
!      Write(6,*)'mInc,mDec,mMat(0:',msiz,',',nosc,')'             ! CGGt
!      Write(6,*)'-------------------------------------------   '  ! CGGt
!!
!!---- Initialize.

!      Write(6,*)                                                  ! CGGt
!      Write(6,*)'                 msiz,nosc==',msiz,nosc          ! CGGt
!      Call GetMem('Test_1','LIST','INTE',iDum,iDum)               ! CGGt
      do iv=0,msiz
         do jv=1,nosc
!            mInc(iv,jv) = 0
!            mDec(iv,jv) = 0
            mMat(iv,jv) = 0
         enddo
      enddo
      If ( m_max.eq.0 ) Return
!      Write(6,*)'                 Calling TabDim_drv'             ! CGGt
      Call TabDim_drv(m_max,nOsc,nTabDim)
      maxOrd = nTabDim-1
!!
!!---- Set up the vertex table
!      Write(6,*)'                 Graph1,dim=',(m_max+1)*(nOsc+1) ! CGGt
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
!      Write(6,*)'                 ipNumber,dim=',(m_max)          ! CGGt
      Call GetMem('Number','Allo','INTE',ipNumber,m_max+1)
      do iv=0,m_max
        iWork(ipNumber+iv) = 0
      enddo
      N = 0
      Do m = 1,m_max
        N = N+Graph1(m,nosc+1)
        iWork(ipNumber+m) = n
      End Do
!      Write(6,*)'                 Graph2,dim=',nOsc*((m_max+1)**2) ! CGG
      do iv=1,m_max+1
        do jv=1,m_max+1
          do kv=1,nOsc
            Graph2(iv,jv,kv) = 0
          enddo
        enddo
      enddo
      Do iOsc = 1,nosc
        Do iQ1 = 0,m_max         ! Where we are going
          Do iQ2 = 0,iQ1-1       ! Where we came from
            Do i = iQ2+1,iq1     ! Sum over preceding paths
              Graph2(iQ1+1,iQ2+1,iOsc) = Graph1(i+1,iOsc)+              &
     &             Graph2(iQ1+1,iQ2+1,iOsc)
            End Do
          End Do
        End Do
      End Do
!!
      Do iQ1 = 0,m_max            ! Where we are going
        Do iQ2 = 0,iq1           ! Where we came from
          Graph2(iQ1+1,iQ2+1,nOsc) = Graph2(iQ1+1,iQ2+1,nOsc)+          &
     &       iWork(ipNumber+iQ1)
        End Do
      End Do
!!
      Call GetMem('Number','Free','INTE',ipNumber,m_max+1)
!!
!!
!      Write(6,*)'CGGt[MakeTab2_a] Vec,dim=',nOsc                  ! CGGt
      Call GetMem('ivec','Allo','INTE',ipiVec,nOsc)
      Do iQuanta=1,m_max
!      Write(6,*)'CGGt[MakeTab2_a] iQuanta=',iQuanta               ! CGGt
        do iv=1,nOsc
          iWork(ipiVec+iv-1)=0
        enddo
        iQ=-1
        iWork(ipiVec)=-1

!      Write(6,*)'CGGt[MakeTab2_a] Call TabDim2_drv - 1 '           ! CGG
!      Write(6,*)'    iQuanta,nOsc,nd==',iQuanta,nOsc,nd            ! CGG
        Call TabDim2_drv(iQuanta,nOsc,nd)
!      Write(6,*)'CGGt[MakeTab2_a] Call TabDim2_drv - 2 '           ! CGG
!      Write(6,*)'     iQuanta-1,nOsc,nvTabDim==',iQuanta-1,nOsc,nvTabDim
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
  99      Continue
          iWork(ipiVec+nOsc-1)=iQuanta-iq
          iDNR=iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
          do iv=1,nOsc
            mMat(iDNR,iv)=iWork(ipiVec+iv-1)
          enddo
        End Do
      End Do

      Call GetMem('iVec','Free','INTE',ipiVec,nOsc)
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(MaxIncOrd)
      End


      Subroutine Mk_nIncDec(m_max,nOrd,msiz,mInc,mDec,mMat,             &
     &                      Graph2,nOsc)
!!
      Implicit Real*8 ( a-h,o-z )
      Integer mInc(0:msiz,nosc), mDec(0:msiz,nosc)
      Integer mMat(0:nOrd,nosc)
      Integer Graph2(m_max+1,m_max+1,nOsc)
      Integer m_max,nOrd,msiz
#include "WrkSpc.fh"
!!
      Call GetMem('iVec','Allo','INTE',ipiVec,nOsc)
!!
!!---- Create mInc.
!!
      do iv=0,msiz
        do jv=1,nosc
          mInc(iv,jv) = -1
        enddo
      enddo
      Do i = 0, msiz
        do iv=1, nOsc
          iWork(ipiVec+iv-1) = mMat(i,iv)
        enddo
        Do j = 1,nOsc
          iWork(ipiVec+j-1) = iWork(ipiVec+j-1)+1
          mInc(i,j)=iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
          iWork(ipiVec+j-1) = iWork(ipiVec+j-1)-1
        End Do
      End Do
!!
!!---- Create mDec.
!!
      do iv=1,nOsc
        mDec(0,iv)=-1
      enddo
      Do i = 1, msiz
        Do j = 1,nOsc
          If (mMat(i,j).ne.0)Then
            do iv=1,nOsc
              iWork(ipiVec+iv-1) = mMat(i,iv)
            enddo
            iWork(ipiVec+j-1)=iWork(ipiVec+j-1)-1
            mDec(i,j)=iDetnr(iWork(ipiVec),Graph2,nosc,m_max)
            do iv=1,nOsc
              iWork(ipiVec+iv-1) = iWork(ipiVec+j-1)+1
            enddo
          Else
            mDec(i,j)=-1
          End IF
        End Do
      End Do
!!
      Call GetMem('ivec','Free','INTE',ipiVec,nOsc)
      Return
      End
