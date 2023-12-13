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

module OneDat

use Definitions, only: iwp, ItoB

implicit none
private

!----------------------------------------------------------------------*
!                                                                      *
! Return codes:                                                        *
!   %good - No error                                                   *
!   %CL01 - file is not opened                                         *
!   %RD03 - information not available                                  *
!   %WR11 - to many operators                                          *
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
!   nAuxDt - extra space for origin and nuclear contribution.          *
!   nTitle - length of title.                                          *
!   MxOp   - Max number of operators                                   *
!   LenOp  - Length of TOC field defining operator                     *
!                                                                      *
! Pointers:                                                            *
!   pFID   - File identifier                                           *
!   pVersN - Version number                                            *
!   pTitle - Titleof the problem                                       *
!   pOp    - Operator list                                             *
!   pSym   - Number of irred. representations                          *
!   pSymOp - generator of irred. representation                        *
!   pBas   - Number of basis functions per irred. representation       *
!   pAtom  - Atoms in system                                           *
!   pCoord - Coordinates of atoms in system                            *
!   pPot   - Nuclear-Nuclear repulsion                                 *
!   pCoM   - Coordinates of center of mass                             *
!   pCoC   - Coordinates of center of charge                           *
!   pALbl  - Atom labels                                               *
!   pType  - Basis function symmetry labels                            *
!   pChrge - Charge of atoms in system                                 *
!   pOption- Various options - direct, expert mode, ...                *
!   pNext  - Next free record                                          *
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
  integer(kind=iwp) :: ID = 4101, VN = 1024
end type FInfo_type
type(FInfo_type), parameter :: FInfoOne = FInfo_type()

type rc_type
  integer(kind=iwp) :: good = 0, &
                       CL01 = 1, &
                       RD03 = 2, &
                       WR11 = 3
end type rc_type
type(rc_type), parameter :: rcOne = rc_type()

integer(kind=iwp), parameter :: sOpSiz = 0, &
                                sNoOri = sOpSiz +1, &
                                sNoNuc = sNoOri +1, &
                                sRdFst = sNoNuc +1, &
                                sRdNxt = sRdFst +1, &
                                sRdCur = sRdNxt +1, &
                                sDbg   = sRdCur +1
integer(kind=iwp), parameter :: sNew = 0, &
                                sDmp = sNew +1
integer(kind=iwp), parameter :: NaN = -1, &
                                nAuxDt = 4, nTitle = (72*2)/ItoB, &
                                MxOp = 16384, LenOp = 5
integer(kind=iwp), parameter :: pFID    = 1, &
                                pVersN  = pFID    +1, &
                                pTitle  = pVersN  +1, &
                                pOp     = pTitle  +nTitle+1, &
                                pSym    = pOp     +MxOp*LenOp, &
                                pSymOp  = pSym    +1, &
                                pBas    = pSymOp  +MxSym, &
                                pAtom   = pBas    +MxSym, &
                                pCoord  = pAtom   +1, &
                                pPot    = pCoord  +6*MxAtom+1, &
                                pCoM    = pPot    +2+1, &
                                pCoC    = pCoM    +6+1, &
                                pALbl   = pCoC    +6+1, &
                                pType   = pALbl   +MxAtom+1, &
                                pChrge  = pType   +4*MxBas+1, &
                                pIndex  = pChrge  +2*MxAtom+1, &
                                pNext   = pIndex  +2*MxAtom+1, &
                                pOption = pNext   +1, &
                                pEnd    = pOption +1
integer(kind=iwp), parameter :: oLabel = 0, &
                                oComp  = oLabel +2, &
                                oSymLb = oComp  +1, &
                                oAddr  = oSymLb +1

type AuxOne_type
  integer(kind=iwp) :: Lu = 0
  logical(kind=iwp) :: Opn = .false.
end type AuxOne_type

integer(kind=iwp), parameter :: lTocOne = 1024*int(real(pEnd+1023)/1024.0)
integer(kind=iwp) :: nBas(8), nSym
type(AuxOne_type) :: AuxOne
integer(kind=iwp), allocatable :: TocOne(:)

public :: AuxOne, FInfoOne, LenOp, lTocOne, MxOp, NaN, nAuxDt, nBas, nSym, oAddr, oComp, oLabel, oSymLb, pALbl, pAtom, pBas, &
          pChrge, pCoC, pCoM, pCoord, pEnd, pFID, pIndex, pNext, pOp, pOption, pPot, pSym, pSymOp, pTitle, pType, pVersN, rcOne, &
          sDbg, sDmp, sNew, sNoNuc, sNoOri, sOpSiz, sRdCur, sRdFst, sRdNxt, TocOne

end module OneDat
