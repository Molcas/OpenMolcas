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
!    Purpose: Create the table of contents of the OrdInt file          *
!                                                                      *
!    Called from: Sort0 and Sort3                                      *
!                                                                      *
!    Calls to : DaFile                                                 *
!                                                                      *
!    Calling parameters:                                               *
!    iDisk  : First free entry on disk after table of contents         *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use TwoDat, only: AuxTwo, FInfoTwo, iNoNum, isBas, isDAdr, isId, isForm, isMxDa, isOrd, isPkPa, isPkTh, isSkip, isSym, isVer, &
                  lTocTwo, nBatch, TocTwo
use Pack_mod, only: isPack, PkThrs
use sort_data, only: iDVBin, LuTwo, mDaTwo, mxSyP, nBs, nSkip, nSln, nSyOp, Square
use Gateway_global, only: iPack
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iDisk
integer(kind=iwp) :: iBatch, iBin, iFlit, ijPair, iOpt, iSyblj, iSyBlk, iSymi, iSymj, jFlit, jSymj, kFlit, klPair, kSybll, kSymk, &
                     kSyml, kSymMx, lFlit, lSyml, lSymMx, nPairs

!----------------------------------------------------------------------*
!     Initialize table of content                                      *
!----------------------------------------------------------------------*

TocTwo(:) = iNoNum

!----------------------------------------------------------------------*
!     Write file identifier                                            *
!----------------------------------------------------------------------*

TocTwo(isId) = FInfoTwo%ID
TocTwo(isVer) = FInfoTwo%VN
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
iPack = merge(0,1,isPack)
TocTwo(isPkPa) = iPack

!----------------------------------------------------------------------*
!     Transfer table of content to disk                                *
!----------------------------------------------------------------------*

iDisk = 0
iOpt = 1
LuTwo = AuxTwo%Unt
call iDAFILE(LuTwo,iOpt,TocTwo,lTocTwo,iDisk)
AuxTwo%DaDa = iDisk

return

end subroutine MkOrd
