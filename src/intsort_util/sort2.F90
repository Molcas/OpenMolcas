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
! Copyright (C) 1993,1996, Markus P. Fuelscher                         *
!               1993, Per Ake Malmqvist                                *
!               1998, Roland Lindh                                     *
!***********************************************************************

subroutine SORT2()
!***********************************************************************
!                                                                      *
!     Purpose: Control phase 2 of bin sorting algorithm. First,        *
!              reload all integrals belonging to the same slice        *
!              and sort them ( cf. SORT2A ). While reading them        *
!              keep a trace of the records which have been read.       *
!              They are used to generate forward chaining indices.     *
!              Also integrals which are of zero value are placed       *
!              back into the list. In the second step put the sorted   *
!              integrals back onto disk. If new records are needed     *
!              append them to the end of the list.                     *
!                                                                      *
!     Called from: Seward_main                                         *
!                                                                      *
!     Calls to : Sort2A,Sort2B,Sort3                                   *
!                                                                      *
!     Calling parameters: none                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.-AA. Malmqvist                             *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     - modified to possibly use a virtual disk                        *
!       M. P. Fuelscher,University of Lund, Sweden, 1996               *
!     - modified to compute MxOrd                                      *
!       R. Lindh,University of Lund, Sweden, 1998                      *
!                                                                      *
!***********************************************************************

use TwoDat, only: lStRec, nBatch, RAMD
use sort_data, only: IndBin, lSll, LuTwo, MxOrd, mxSyP, nBs, nSkip, nSln, nSyOp, Square, ValBin
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: ib, iBatch, iBin, ibj, iDisk, iErr, iOff, iOpt, iOrd, iPrint, iRout, iSkip, iSlice, iStk, iSyblj, iSyBlk, &
                     iSymi, iSymj, jb, jSkip, jSymj, kb, kbl, kSkip, kSybll, kSymk, kSyml, kSymMx, lb, lSkip, lSlice, lSrtA, &
                     lSrtA_, lStk, lStk_Max, lSyml, lSymMx, mxij, nij, nSlice, nStk, nSym
integer(kind=iwp), allocatable :: IOStk(:)
real(kind=wp), allocatable :: SrtA(:), Scr(:)

!----------------------------------------------------------------------*
!     pick up the print level                                          *
!----------------------------------------------------------------------*

iRout = 84
iPrint = nPrint(iRout)
if (iPrint >= 10) write(u6,*) ' >>> Enter SORT2 <<<'

!----------------------------------------------------------------------*
!     Initialize the IO-stack                                          *
!----------------------------------------------------------------------*

call mma_MaxINT(lStk_Max)
lStk = max(64*1024,lStk_Max/5)
call mma_allocate(IOStk,lStk,Label='IOStk')
IOStk(:) = 0
nStk = 0

!----------------------------------------------------------------------*
!     loop over all submatrices of 2el integrals (slices)              *
!----------------------------------------------------------------------*

iOrd = 0
iBin = 0
nSym = nSyOp
do iSymi=1,nSym
  ib = nBs(iSymi)
  iSkip = nSkip(iSymi)
  do jSymj=1,iSymi
    iSymj = 1+ieor(iSymi-1,jSymj-1)
    iSyblj = jSymj+iSymi*(iSymi-1)/2
    jb = nBs(jSymj)
    ibj = ib*jb
    if (iSymi == jSymj) ibj = ib*(ib+1)/2
    jSkip = nSkip(jSymj)
    kSymMx = iSymi
    if (Square) kSymMx = nSym
    do kSymk=1,kSymMx
      kb = nBs(kSymk)
      kSkip = nSkip(kSymk)
      lSymMx = jSymj
      if ((kSymk /= iSymi) .or. Square) lSymMx = kSymk
      do lSyml=1,lSymMx
        kSyml = 1+ieor(kSymk-1,lSyml-1)
        kSybll = lSyml+kSymk*(kSymk-1)/2
        if (ieor(iSymj-1,kSyml-1) == 0) then
          lb = nBs(lSyml)
          kbl = kb*lb
          if (kSymk == lSyml) kbl = kb*(kb+1)/2
          lSkip = nSkip(lSyml)
          iSyBlk = kSybll+mxSyP*(iSyblj-1)

          if ((iSkip+jSkip+kSkip+lSkip == 0) .and. (ibj*kbl /= 0)) then
            !                                                          *
            !***********************************************************
            !                                                          *
            if (RAMD%act) then

              ! RAMD option

              lSrtA = ibj*kbl
              iBin = iBin+1
              iBatch = nBatch(iSyBlk)
              iOff = RAMD%adr(iBatch)+1
              call SORT2B(iBin,lSrtA,iOrd,lSrtA,RAMD%ints(iOff),IOStk,lStk,nStk)
              !                                                        *
              !*********************************************************
              !                                                        *
            else
              !                                                        *
              !*********************************************************
              !                                                        *
              ! Sorting option
              !
              nSlice = nSln(iSyBlk)
              lSlice = lSll(iSyBlk)
              nij = lSlice/kbl
              mxij = ibj

              lSrtA_ = min(mxij,nij)*kbl
              call mma_allocate(SrtA,lSrtA_,Label='SrtA')

              !--------------------------------------------------------*
              !     Allocate space to keep all integrals of a slice in *
              !     memory fetch them from disk and sort them.         *
              !--------------------------------------------------------*

              do iSlice=1,nSlice
                iBin = iBin+1
                lSrtA = min(mxij,nij)*kbl
                SrtA(1:lSrtA) = Zero
                call SORT2A(iBin,lSrtA,SrtA,IOStk,lStk,nStk)

                !------------------------------------------------------*
                ! Sort the IO-stack in ascending order, i.e., prefer-  *
                ! ence is always given to the lowest available disk    *
                ! adresses.                                            *
                !------------------------------------------------------*

                call ILASRT('D',nStk,IOStk,iErr)

                !------------------------------------------------------*
                ! Restore integrals on LuTwo                           *
                !------------------------------------------------------*

                call SORT2B(iBin,lSrtA,iOrd,lSrtA,SrtA,IOStk,lStk,nStk)
                mxij = mxij-nij
              end do

              call mma_deallocate(SrtA)
              !                                                        *
              !*********************************************************
              !                                                        *
            end if
            !                                                          *
            !***********************************************************
            !                                                          *
          end if
        end if
      end do
    end do
  end do
end do

!----------------------------------------------------------------------*
!     Clean the remaining records in the I/O Stack                     *
!----------------------------------------------------------------------*

call mma_allocate(Scr,lStRec,Label='Scr')
Scr(:) = Zero
do iStk=1,nStk
  iOrd = iOrd+1
  iDisk = IOStk(iStk)
  Scr(2) = real(iOrd,kind=wp)
  iOpt = 1
  call dDAFILE(LuTwo,iOpt,Scr,lStRec,iDisk)
end do
MxOrd = iOrd
call mma_deallocate(Scr)
call mma_deallocate(IOStk)

if (.not. RAMD%act) then
  call mma_deallocate(ValBin)
  call mma_deallocate(IndBin)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SORT2
