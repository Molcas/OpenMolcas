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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!               1991, Per Ake Malmqvist                                *
!***********************************************************************

subroutine MkOrd(iDisk)
!***********************************************************************
!                                                                      *
!    Purpose: Create the table of content of the OrdInt file           *
!                                                                      *
!    Called from: Sort0 ans Sort3                                      *
!                                                                      *
!    Calls to : DaFile                                                 *
!                                                                      *
!    Calling parameters:                                               *
!    iDisk  : First free entry on disk after table of contents         *
!                                                                      *
!    Global data declarations (Include files) :                        *
!    TwoDat  : table of contents and auxiliary information             *
!              on the ordered 2el file                                 *
!    TowID   : Table of file identifiers                               *
!    TwoDef  : definitions of the record structure                     *
!    Srt0    : common block containing information pertinent to        *
!              the calculation of 2el integral sequence numbers        *
!    Srt1    : common block containing information the number of       *
!              bins and partitioning of symmetry blocks                *
!    Srt2    : common block containing information pertinent to        *
!              the bin sorting algorithm                               *
!    PkCtl   : packing table                                           *
!                                                                      *
!    Local data declarations: none                                     *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use srt2, only: LuTwo, mDaTwo, iDVBin
use Integral_Parameters, only: iPack
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iDisk
#include "FileIDs.fh"
#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "PkCtl.fh"
integer(kind=iwp) :: iAssm, iBatch, iBin, iExp, iFlit, ijPair, Init_do_setup_d, Init_do_setup_e, Init_do_setup_l, iOpt, iSyblj, &
                     iSyBlk, iSymi, iSymj, iToc, jFlit, jSymj, kFlit, klPair, kSybll, kSymk, kSyml, kSymMx, lFlit, lSyml, lSymMx, &
                     nPairs

!----------------------------------------------------------------------*
!     Initialize table of content                                      *
!----------------------------------------------------------------------*

do iToc=1,lTocTwo
  TocTwo(iToc) = iNoNum
end do

!----------------------------------------------------------------------*
!     Write file identifier                                            *
!----------------------------------------------------------------------*

TocTwo(isId) = IDtwo
TocTwo(isVer) = VNtwo
TocTwo(isForm) = 0

!----------------------------------------------------------------------*
!     Write ordring mode                                               *
!----------------------------------------------------------------------*

TocTwo(isOrd) = 0

!----------------------------------------------------------------------*
!     Write symmetry and basis set information                         *
!----------------------------------------------------------------------*

TocTwo(isSym) = nSyOp
do iSymi=0,nSyOp-1
  TocTwo(isSkip+iSymi) = nSkip(iSymi+1)
end do
nPairs = nSyOp*(nSyOp+1)/2
iBatch = 0
do iSymi=1,nSyOp
  do jSymj=1,iSymi
    do kSymk=1,nSyOp
      do lSyml=1,kSymk
        if (ieor(iSymi-1,jSymj-1) == ieor(kSymk-1,lSyml-1)) then
          ijPair = jSymj+iSymi*(iSymi-1)/2
          klPair = lSyml+kSymk*(kSymk-1)/2
          iSyBlk = (ijPair-1)*nPairs+klPair
          iBatch = iBatch+1
          nBatch(iSyBlk) = iBatch
        end if
      end do
    end do
  end do
end do

!----------------------------------------------------------------------*
!     Write the skip parameters                                        *
!----------------------------------------------------------------------*

TocTwo(isSym) = nSyOP
do iSymi=0,nSyOp-1
  TocTwo(isBas+iSymi) = nBs(iSymi+1)
end do

!----------------------------------------------------------------------*
!     Write disk access table                                          *
!----------------------------------------------------------------------*

do iBatch=0,175
  TocTwo(isDAdr+iBatch) = 0
end do

iBin = 1
do iSymi=1,nSyOp
  iFlit = nSkip(iSymi)
  do jSymj=1,iSymi
    jFlit = nSkip(jSymj)
    iSymj = 1+ieor(iSymi-1,jSymj-1)
    iSyblj = jSymj+iSymi*(iSymi-1)/2
    kSymMx = iSymi
    if (Square) kSymMx = nSyOp
    do kSymk=1,kSymMx
      kFlit = nSkip(kSymk)
      lSymMx = jSymj
      if ((kSymk /= iSymi) .or. Square) lSymMx = kSymk
      do lSyml=1,lSymMx
        lFlit = nSkip(lSyml)
        kSyml = 1+ieor(kSymk-1,lSyml-1)

        kSybll = lSyml+kSymk*(kSymk-1)/2
        if (ieor(iSymj-1,kSyml-1) == 0) then
          if ((iFlit+jFlit+kFlit+lFlit) == 0) then
            iSyBlk = kSybll+mxSyP*(iSyblj-1)
            iBatch = nBatch(iSyBlk)
            TocTwo(isDAdr+iBatch-1) = iDVBin(2,iBin)
            iBin = iBin+nSln(iSyBlk)
          end if
        end if
      end do
    end do
  end do
end do

TocTwo(isMxDa) = mDaTwo

!----------------------------------------------------------------------*
!     Write packing information                                        *
!----------------------------------------------------------------------*

call Real2Int(PkThrs,TocTwo(isPkTh))
call Real2Int(PkCutof,TocTwo(isPkCt))
call Real2Int(PkScal,TocTwo(isPkSc))
iPack = 1
if (Pack) iPack = 0
TocTwo(isPkPa) = iPack
iAssm = 1
if (Assm) iAssm = 0
TocTwo(isPkAs) = iAssm
do iExp=0,4095
  TocTwo(isPkTb+iExp) = PkTab(iExp)
end do

!----------------------------------------------------------------------*
!     Transfer table of content to disk                                *
!----------------------------------------------------------------------*

iDisk = 0
iOpt = 1
LuTwo = AuxTwo(isUnit)
call iDAFILE(LuTwo,iOpt,TocTwo,lTocTwo,iDisk)
AuxTwo(isDaDa) = iDisk

return

end subroutine MkOrd
