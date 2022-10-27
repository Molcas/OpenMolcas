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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine LoadInts(iRc,iOpt)
!***********************************************************************
!                                                                      *
!    Purpose: Read the table of content from the OrdInt file           *
!                                                                      *
!    Calling parameters:                                               *
!    iRc    : return code                                              *
!             (zero upon successful completion)                        *
!    iOpt   : option code                                              *
!             iOpt = 0  ==> Square =.true.                             *
!             iOpt = 1  ==> Square =.false.                            *
!                                                                      *
!    Global data declarations (Include files) :                        *
!    TwoDat : table of contents and auxiliary information              *
!    TowRc  : Table of return code                                     *
!                                                                      *
!    Local data declarations: none                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul

implicit integer(A-Z)
#include "Molcas.fh"
#include "TwoDat.fh"
logical Square

! loop over all symmetry blocks

Square = iOpt == 1
iOff = RAMD_anchor
nInts = 0
nSym = TocTwo(isSym)
mxSyP = nSym*(nSym+1)/2
do iSymi=1,nSym
  ib = TocTwo(isBas+iSymi-1)
  iSkip = TocTwo(isSkip+iSymi-1)
  do jSymj=1,iSymi
    iSymj = Mul(iSymi,jSymj)
    iSyblj = jSymj+iSymi*(iSymi-1)/2
    jb = TocTwo(isBas+jSymj-1)
    ibj = ib*jb
    if (iSymi == jSymj) ibj = ib*(ib+1)/2
    jSkip = TocTwo(isSkip+jSymj-1)
    kSymMx = iSymi
    if (Square) kSymMx = nSym
    do kSymk=1,kSymMx
      kb = TocTwo(isBas+kSymk-1)
      kSkip = TocTwo(isSkip+kSymk-1)
      lSymMx = jSymj
      if ((kSymk /= iSymi) .or. Square) lSymMx = kSymk
      do lSyml=1,lSymMx
        kSyml = Mul(kSymk,lSyml)
        kSybll = lSyml+kSymk*(kSymk-1)/2
        if (Mul(iSymj,kSyml) == 1) then
          lb = TocTwo(isBas+lSyml-1)
          kbl = kb*lb
          if (kSymk == lSyml) kbl = kb*(kb+1)/2
          lSkip = TocTwo(isSkip+lSyml-1)
          iSyBlk = kSybll+mxSyP*(iSyblj-1)
          if ((iSkip+jSkip+kSkip+lSkip == 0) .and. (ibj*kbl /= 0)) then

            ! check if there is enough space available

            lBuf = ibj*kbl
            nInts = nInts+lBuf
            if (nInts >= RAMD_size) then
              iRc = 001
              write(6,*)
              write(6,'(2X,A,I3.3,A)') '*** (W)-level message LOADINTS',iRc,' ***'
              write(6,'(2X,A)') 'There is not enough space on the RAM disk'
              write(6,'(2X,A)') 'The program will resume normal activity'
              write(6,*)
              return
            end if

            ! save start address of this symmetry block

            iBatch = nBatch(iSyBlk)
            RAMD_adr(iBatch) = iOff

            ! load integrals

            iRc = 0
            iOpt = 1
            call RdOrd(iRc,iOpt,iSymi,jSymj,kSymk,lSyml,RAMD_ints(iOff),lBuf+1,npq)

            ! update pointers

            iOff = iOff+lBuf

          end if
        end if
      end do
    end do
  end do
end do

! define initial load point

RAMD_next = RAMD_adr(1)

return

end subroutine LoadInts
