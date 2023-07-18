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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************
      SubRoutine ChoMP2_Vec(iVec1,nVec,Buf,lBuf,nDim,iOpt)
!
!     Thomas Bondo Pedersen, Dec. 2004.
!
!     Purpose: write (iOpt=1) or read (iOpt=2) "new" vectors to buffer.
!
      use ChoMP2_dec, only: NowSym
      Implicit None
      Integer iVec1, nVec, lBuf, nDim, iOpt
      Real*8  Buf(lBuf)
#include "chomp2.fh"

      Character(LEN=3), Parameter::  ThisNm = 'Vec'
      Character(LEN=10), Parameter:: SecNam = 'ChoMP2_Vec'
      Integer :: iSym, iJob, lTot, iAdr

      Logical DoClose

      iSym = NowSym
      DoClose = .false.

      If (iOpt .eq. 1) Then

         If (lUnit_F(iSym,2) .lt. 1) Then
            Call ChoMP2_OpenF(1,2,iSym)
            DoClose = .true.
         End If

         iJob = 1
         lTot = nDim*nVec
         iAdr = nDim*(iVec1 - 1) + 1
         Call ddaFile(lUnit_F(iSym,2),iJob,Buf,lTot,iAdr)

      Else If (iOpt .eq. 2) Then

         If (lUnit_F(iSym,2) .lt. 1) Then
            Call ChoMP2_OpenF(1,2,iSym)
            DoClose = .true.
         End If

         iJob = 2
         lTot = nDim*nVec
         iAdr = nDim*(iVec1 - 1) + 1
         Call ddaFile(lUnit_F(iSym,2),iJob,Buf,lTot,iAdr)

      Else

         Write(6,*) SecNam,': illegal option: iOpt = ',iOpt
         Call ChoMP2_Quit(SecNam,'illegal option',' ')

      End If

      If (DoClose) Then
         Call ChoMP2_OpenF(2,2,iSym)
         DoClose = .false.
      End If

      End
