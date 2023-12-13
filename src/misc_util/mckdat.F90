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

module MckDat

use Definitions, only: iwp, ItoB

implicit none
private

!----------------------------------------------------------------------*
!                                                                      *
! Return codes:                                                        *
!   %good - No error                                                   *
!   %CL01 - file is not opened                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Switches for one electron integral handlers:                         *
!   sOpSiz - Return only size of array                                 *
!   sNoOri - Do not read origin of operator                            *
!   sNoNuc - Do not read nuclear contribution                          *
!   sRdFst - Read first operator                                       *
!   sRdNxt - Read next operator                                        *
!   sNew   - New file, create.                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Other definitions                                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
!   NaN    - Not a Number = Variable undefined                         *
!   NotNaN - A Number = Variable is defined                            *
!   nTitle - length of title.                                          *
!   PhyRec - physical buffer size for DAFILE                           *
!   nBuf   - logical internal buffer size for reading/writing matrices *
!   MxOp   - Max number of operators                                   *
!   LenOp  - Length of TOC field defining operator                     *
!                                                                      *
! Pointers:                                                            *
!   pFID     - File identifier                                         *
!   pVersN   - Version number                                          *
!   pTitle   - Titleof the problem                                     *
!   pOp      - Operator list                                           *
!   pSym     - Number of irred. representations                        *
!   pSymOp   - generator of irred. representation                      *
!   pBas     - Number of basis functions per irred. representation     *
!   pASh     -   ?                                                     *
!   pish     -   ?                                                     *
!   pChdisp  -   ?                                                     *
!   pndisp   -   ?                                                     *
!   pldisp   -   ?                                                     *
!   pnrdisp  -   ?                                                     *
!   pdegdisp -   ?                                                     *
!   ptdisp   -   ?                                                     *
!   pPert    -   ?                                                     *
!   pNext    - Next free record                                        *
!                                                                      *
! Offsets:                                                             *
!   oLabel - Label of operator                                         *
!   oComp  - Component number                                          *
!   oSymLb - Symmetry label of operator                                *
!   oAddr  - Disk address                                              *
!                                                                      *
!----------------------------------------------------------------------*

#include "Molcas.fh"

type FInfo_type
  integer(kind=iwp) :: ID = 4097, VN = 1024
end type FInfo_type
type(FInfo_type), parameter :: FInfoMck = FInfo_type()

type rc_type
  integer(kind=iwp) :: good = 0, &
                       CL01 = 1
end type rc_type
type(rc_type), parameter :: rcMck = rc_type()

integer(kind=iwp), parameter :: sOpSiz  =          0, &
                                sNoOri  = sOpSiz  +1, &
                                sNoNuc  = sNoOri  +1, &
                                sRdFst  = sNoNuc  +1, &
                                sRdNxt  = sRdFst  +1, &
                                sRdCur  = sRdNxt  +1, &
                                sLength = sRdCur  +1, &
                                sDbg    = sLength +1
integer(kind=iwp), parameter :: sNew = 0, &
                                sDmp = sNew +1
integer(kind=iwp), parameter :: NaN = -1, NotNaN = 0, &
                                nTitle = (72*2)/ItoB, &
                                PhyRec = 1024, nBuf = 4*PhyRec, MxOp = 2048, LenOp = 5
integer(kind=iwp), parameter :: pFID     = 1, &
                                pVersN   = pFID     +1, &
                                pTitle   = pVersN   +1, &
                                pOp      = pTitle   +nTitle+1, &
                                pSym     = pOp      +MxOp*LenOp, &
                                pSymOp   = pSym     +1, &
                                pBas     = pSymOp   +int(real(3*MxSym+ItoB-1)/real(ItoB)), &
                                pASh     = pBas     +MxSym, &
                                pish     = pAsh     +MxSym, &
                                pChdisp  = pish     +MxSym, &
                                pndisp   = pchdisp  +5*MxOp, &
                                pldisp   = pndisp   +1, &
                                pnrdisp  = pldisp   +MxSym, &
                                pdegdisp = pnrdisp  +MxOp, &
                                ptdisp   = pdegdisp +MxOp, &
                                pPert    = Ptdisp   +MxOp, &
                                pNext    = ppert    +5, &
                                pEnd     = pNext    +1
integer(kind=iwp), parameter :: oLabel = 0, &
                                oComp  = oLabel +2, &
                                oSymLb = oComp  +1, &
                                oAddr  = oSymLb +1

type AuxMck_type
  integer(kind=iwp) :: Lu = 0
  logical(kind=iwp) :: Opn = .false.
end type AuxMck_type

integer(kind=iwp), parameter :: lTocMck = 1024*int(real(pEnd+1023)/1024.0)
integer(kind=iwp) :: TmpBuf(nBuf), TocMck(lTocMck)
type(AuxMck_type) :: AuxMck

public :: AuxMck, FInfoMck, LenOp, lTocMck, MxOp, NaN, nBuf, NotNaN, nTitle, oAddr, oComp, oLabel, oSymLb, pASh, pBas, pChdisp, &
          pdegdisp, pEnd, pFID, pish, pldisp, pndisp, pNext, pnrdisp, pOp, pPert, pSym, pSymOp, ptdisp, pTitle, pVersN, rcMck, &
          sLength, sDbg, sDmp, sNew, sNoNuc, sNoOri, sOpSiz, sRdCur, sRdFst, sRdNxt, TmpBuf, TocMck

end module MckDat
