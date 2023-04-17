!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module TwoDat

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

!----------------------------------------------------------------------*
!                                                                      *
! Return codes:                                                        *
!   %good - No error                                                   *
!   %CL01 - file has not been opened                                   *
!   %TC01 - file has not been opened                                   *
!   %TC02 - unknown ordering parameter                                 *
!   %TC03 - illegal number of symmetry operations                      *
!   %RD01 - illegal combination of symmetry labels                     *
!   %RD02 - illegal ordering of symmetry labels                        *
!   %RD03 - inconsistent symmetry batch                                *
!   %RD04 - illegal buffer size                                        *
!   %RD05 - buffer size is to small                                    *
!   %RD06 - unknown option                                             *
!   %RD07 - symmetry block is not available                            *
!   %RD08 - file has not been opened                                   *
!   %RD09 - packing table has not been loaded                          *
!   %RD10 - pq1 out of bounds                                          *
!   %RD11 - bad initialization of genint                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Switches for OpnOrd:                                                 *
!   sNew   - create a new file                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Entry points to auxiliary information on 2el. file:              *
!     (AuxTwo)                                                         *
!                                                                      *
!     %Unt  : logical unit number of 2el. integral file                *
!     %Opn  : status identifier (is it opened?)                        *
!     %Dada : disk address of last access                              *
!     %Upk8 : pointer to next integral to unpack                       *
!     %lBf1 : number of left integrals in TWOBUF                       *
!     %Npq  : number of pq pairs per symmetry batch which has not      *
!             been read currently                                      *
!                                                                      *
!     Entry points to table of contents for ordered 2el. integrals:    *
!                                                                      *
!     lTocTwo : length of table of content                             *
!     isId    : location of file identifier                            *
!     isVer   : location of version number                             *
!     isOrd   : location of ordering key                               *
!     isForm  : format of data on disk                                 *
!     isSym   : location of number of irred. representations           *
!     isBas   : first location of basis function counters              *
!     isSkip  : first location of the skip flags                       *
!     isDAdr  : first location of disk adresses for each sym. batch    *
!     isMxDa  : highest disk adress written                            *
!     isPkTh  : location of threshold for packing                      *
!     isPkCt  : location of cutof for packing (unused)                 *
!     isPkSc  : location of scaling constant for packing (unused)      *
!     isPkPa  : location of packing key                                *
!     isPkAs  : location of assmbler key (unused)                      *
!     isPkTb  : first location of packing table (unused)               *
!     isFree  : first free location for suplementary information       *
!                                                                      *
!     Parameters for the symmetry block number to batch number         *
!     translation table                                                *
!                                                                      *
!     lBatch  : length of the translation table                        *
!                                                                      *
!     Dummy number used to initialize AuxTwo and TocTwo                *
!                                                                      *
!     iNoNum  : dummy constant = -1                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! RAMD: save control structures needed if the two-electron integrals   *
!       if random access memory (RAM) is used as disk to keep them.    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Global Parameters used for the I/O                               *
!                                                                      *
!     lDaRec  : minimal record (in 4Byte words) transferred by         *
!               DAFILE                                                 *
!     nSect   : integer multiple of lDaRec which constitutes the       *
!               standard record length                                 *
!     lStRec  : standard record length                                 *
!     lTop    : each record has its own header of length lTop          *
!                                                                      *
!----------------------------------------------------------------------*

type FInfo_type
  integer(kind=iwp) :: ID = 4098, VN = 1024
end type FInfo_type
type(FInfo_type), parameter :: FInfoTwo = FInfo_type()

type rc_type
  integer(kind=iwp) :: good = 0, &
                       CL01 = 1, &
                       TC01 = 2, &
                       TC02 = 3, &
                       TC03 = 4, &
                       RD01 = 5, &
                       RD02 = 6, &
                       RD03 = 7, &
                       RD04 = 8, &
                       RD05 = 9, &
                       RD06 = 10, &
                       RD07 = 11, &
                       RD08 = 12, &
                       RD09 = 13, &
                       RD10 = 14, &
                       RD11 = 15
end type rc_type
type(rc_type), parameter :: rcTwo = rc_type()

integer(kind=iwp), parameter :: sNew = 0
integer(kind=iwp), parameter :: isId    = 1, &
                                isVer   = isId   +1, &
                                isOrd   = isVer  +1, &
                                isForm  = isOrd  +1, &
                                isSym   = isForm +1, &
                                isBas   = isSym  +1, &
                                isSkip  = isBas  +8, &
                                isDAdr  = isSkip +8, &
                                isMxDa  = isDadr +176, &
                                isPkTh  = isMxDa +1, &
                                isPkCt  = isPkTh +2, &   ! unused
                                isPkSc  = isPkCt +2, &   ! unused
                                isPkPa  = isPkSc +2, &
                                isPkAs  = isPkPa +1, &   ! unused
                                isPkTb  = isPkAs +1, &   ! unused
                                isFree  = isPkTb +4096, &
                                lTocTwo = isFree +10
integer(kind=iwp), parameter :: iNoNum = -1, lBatch = 1296, lDaRec = 1024, nSect = 32, lStRec = nSect*lDaRec, lTop = 4

type AuxTwo_type
  integer(kind=iwp) :: Unt = 0, DaDa = 0, Upk8 = 0, lBf1 = 0, Npq = 0
  logical(kind=iwp) :: Opn = .false.
end type AuxTwo_type

type RAMD_type
  logical(kind=iwp) :: act = .false.
  integer(kind=iwp) :: adr(176) = 0, next = 0
  real(kind=wp) :: ints(8) = Zero
end type RAMD_type

integer(kind=iwp) :: nBatch(lBatch), TocTwo(lTocTwo)
type(AuxTwo_type) :: AuxTwo
type(RAMD_type) :: RAMD

public :: AuxTwo, FInfoTwo, iNoNum, isBas, isDAdr, isForm, isFree, isId, isMxDa, isOrd, isPkPa, isPkTh, isSkip, isSym, isVer, &
          lDaRec, lStRec, lTocTwo, lTop, nBatch, nSect, RAMD, rcTwo, sNew, TocTwo

end module TwoDat
