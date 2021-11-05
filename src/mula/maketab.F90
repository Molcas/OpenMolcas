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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************
      Subroutine MakeTab(m_max,maxOrd,maxIncOrd,mMat,mInc,mDec,nOsc)
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
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Implicit Real*8 ( a-h,o-z )
      Integer mMat  (0:maxord,nosc)
      Integer mInc(0:maxord,nosc),mDec(0:maxord,nosc)
      Logical    equal
#include "WrkSpc.fh"
!!
!!---- Initialize.
      maxOrd = 0
      maxIncOrd = 0
      do iv=0,maxord
      do jv=1,nosc
      mMat(iv,jv) = 0
      mInc(iv,jv) = 0
      mDec(iv,jv) = 0
      enddo
      enddo
      If ( m_max.eq.0 ) Return
      If ( nOsc.eq.1 ) Then
      mTempDim = m_max
      Else
      mTempDim = (1-nOsc**(m_max+1))/(1-nOsc)
      End If
      lmTemp=mTempDim*nOsc
      Call GetMem('mTemp','Allo','INTE',ipmTemp,lmTemp)

      Call GetMem('row','Allo','INTE',iprow,nOsc)

!       mTemp = 0
      do iv=0,lmTemp-1
      iWork(ipmTemp+iv)=0
      enddo
      num = 1
      istart_row = 0
      irow = istart_row+1
      Call GetMem('unit','Allo','INTE',ipunit,nOsc*nOsc)
      Do i = 1,nOsc
      Do j = 1,nOsc
      iWork(ipunit+i+nOsc*(j-1)-1) = 0
      End Do
      iWork(ipunit+i+nOsc*(i-1)-1) = 1
      End Do
!!
!!---- Create table mMat.
      mMat_row = 1
      Do m = 1,m_max
      irow = 1
      numtemp = 0
!!---- Produce all combinations for a given total sum = m_max.
      Do n = 1,num
      Do jrow = 1,nOsc
      Do jcol = 1,nOsc
      iWork(ipmTemp+irow+mTempDim*(jcol-1)) =                           &
     &                   mMat(istart_row,jcol)+                         &
     &                   iWork(ipunit+jrow+nOsc*(jcol-1)-1)
      End Do
      irow = irow+1
      numtemp = numtemp+1
      End Do
      istart_row = istart_row+1
      End Do
      num = numtemp
!!---- Remove all entries which occur more than once.
      Do i = 1,num
!             row1 => mTemp(i,:)
      j = istart_row
      jmax = mMat_row
      equal = .false.
      Do While (( j.lt.jmax ).and.( .not.equal ))
!                row2 => mMat(j,:)
      equal = .true.
      Do l = 1,nOsc
!                   If ( row1(l).ne.row2(l) ) Then
      If ( iWork(ipmTemp+i+mTempDim*(l-1)).ne.                          &
     &                              mMat(j,l) ) Then
      equal = .false.
      End If
      End Do
      j = j+1
      End Do
      If ( .not.equal ) Then
      Do jcol = 1,nOsc
      mMat(mMat_row,jcol) =                                             &
     &                       iWork(ipmTemp+i+mTempDim*(jcol-1))
      End Do
      mMat_row = mMat_row+1
      End If
      End Do
      num = mMat_row-istart_row
      End Do
      maxOrd = mMat_row-1
      Call GetMem('mTemp','Free','INTE',ipmTemp,lmTemp)
!!
!!---- Create mInc.
      maxIncOrd = maxOrd-num
      Do i = 0,maxIncOrd
      do iv=1,nOsc
      iWork(iprow+iv-1) = mMat(i,iv)
      enddo
      Do j = 1,nOsc
      iWork(iprow+j-1) = iWork(iprow+j-1)+1
      equal = .false.
      k = i+1
      Do While (( .not.equal ).and.( k.le.maxOrd ))
!                row2 => mMat(k,:)
      equal = .true.
      Do l = 1,nOsc
      If ( iWork(iprow+l-1).ne.mMat(k,l) ) Then
      equal = .false.
      End If
      End Do
      If ( equal ) Then
      mInc(i,j) = k
      Else
      k = k+1
      End If
      End Do
      iWork(iprow+j-1) = iWork(iprow+j-1)-1
      End Do
      End Do
!!
!!---- Create mDec.
      Do i = 1,maxOrd
      Do iv=1,nOsc
      iWork(iprow+iv-1) = mMat(i,iv)
      enddo
      Do j = 1,nOsc
      If ( iWork(iprow+j-1).gt.0 ) Then
      iWork(iprow+j-1) = iWork(iprow+j-1)-1
      equal = .false.
      k = 0
      Do While (( .not.equal ).and.( k.le.maxOrd ))
!                   row2 => mMat(k,:)
      equal = .true.
      Do l = 1,nOsc
      If ( iWork(iprow+l-1).ne.mMat(k,l) ) Then
      equal = .false.
      End If
      End Do
      If ( equal ) Then
      mDec(i,j) = k
      Else
      k = k+1
      End If
      End Do
      iWork(iprow+j-1) = iWork(iprow+j-1)+1
      Else
      mDec(i,j) = 0
      End If
      End Do
      End Do
!!
      Call GetMem('row','Free','INTE',iprow,nOsc)
      Call GetMem('unit','Free','INTE',ipunit,nOsc*nOsc)
!!
      End
