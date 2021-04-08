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

module motra_global

use Definitions, only: wp, iwp

implicit none
private

#include "mxdm.fh"

!----------------------------------------------------------------------*
! Allocate space to store the system description                       *
!----------------------------------------------------------------------*
integer(kind=iwp) :: iOper(mxSym), nAtoms, nBas(mxSym), nDel(mxSym), nFro(mxSym), nOrb(mxSym), nSym
real(kind=wp) :: Coor(3,MxAtms), PotNuc
character(len=LenIn8) :: BsLbl(mxOrb)

!----------------------------------------------------------------------*
! Allocate space to store the MO-coefficients and occupations          *
!----------------------------------------------------------------------*
real(kind=wp) :: Occ(mxOrb)

!----------------------------------------------------------------------*
! Allocate space to store the one-electron inntegral file header and   *
! the header of the input source of MO coefficients                    *
!----------------------------------------------------------------------*
character(len=80) :: VecTit
character(len=144) :: Header

!----------------------------------------------------------------------*
! Allocate space to store the title                                    *
!----------------------------------------------------------------------*
integer(kind=iwp) :: nTit
character(len=72) :: Title(mxTit)

!----------------------------------------------------------------------*
! Allocate space to store logical switches for routing and printing    *
!----------------------------------------------------------------------*
integer(kind=iwp) :: Debug, iAutoCut, ihdf5, iOneOnly, iPrint, iRFpert, iVecTyp

!----------------------------------------------------------------------*
! Save cutting threshold for AUTO cut option                           *
!----------------------------------------------------------------------*
real(kind=wp) :: CutThrs(mxSym)

!----------------------------------------------------------------------*
! Store file names and unit numbers.                                   *
!----------------------------------------------------------------------*
! InpOrb, JobIph: Input files for MO coefficients
! OneAO, TwoAO:   one- and two-electron integrals in AO basis
! OneMO, TwoMO:   one- and two-electron integrals in MO basis
! Half:           half transformed two-electron integrals
! Ext:            EXTRACT file
! Com:            COMFILE file
integer(kind=iwp) :: LuCom, LuExt, LuHalf, LuInpOrb, LuJobIph, LuOneAO, LuOneMO, LuTwoAO, LuTwoMO
character(len=8) :: FnCom, FnExt, FnHalf, FnOneAO, FnOneMO, FnTwoAO, FnTwoMO
character(len=180) :: FnInpOrb, FnJobIph

!----------------------------------------------------------------------*
! Allocate space to store the table of contents of various files       *
!----------------------------------------------------------------------*
integer(kind=iwp) :: TcJobIph(1024), TcOneMO(1024)

!----------------------------------------------------------------------*
! Define TOC for the electron repulsion integrals in MO basis          *
!                                                                      *
! Define the buffer size of the electron repulsion integrals in MO     *
! basis                                                                *
!----------------------------------------------------------------------*
#include "tratoc.fh"
integer(kind=iwp), parameter :: kBuf = nTraBuf

integer(kind=iwp) :: iCTonly, iDoInt
integer(kind=iwp) :: nTot, nTot1, nTot2, n2max, nOrbt, nOrbtt, ISP, ISQ, ISR, ISS, NBP, NBQ, NBR, NBS, NOP, NOQ, NOR, NOS, LMOP, &
                     LMOQ, LMOR, LMOS, NBPQ, NBRS, NOVX, INCORE, MEMX, LTUVX, IAD13

public :: &
BsLbl, &
Coor, &
CutThrs, &
Debug, &
FnCom, &
FnExt, &
FnHalf, &
FnInpOrb, &
FnJobIph, &
FnOneAO, &
FnOneMO, &
FnTwoAO, &
FnTwoMO, &
Header, &
IAD13, &
iAutoCut, &
iCTonly, &
iDoInt, &
ihdf5, &
INCORE, &
iOneOnly, &
iOper, &
iPrint, &
iRFpert, &
ISP, &
ISQ, &
ISR, &
ISS, &
iTraToc, &
iVecTyp, &
kBuf, &
LenIn8, &
LMOP, &
LMOQ, &
LMOR, &
LMOS, &
LTUVX, &
LuCom, &
LuExt, &
LuHalf, &
LuInpOrb, &
LuJobIph, &
LuOneAO, &
LuOneMO, &
LuTwoAO, &
LuTwoMO, &
MEMX, &
MxOrb, &
MxRoot, &
MxTit, &
MxSym, &
n2max, &
nAtoms, &
nBas, &
NBP, &
NBPQ, &
NBQ, &
NBR, &
NBRS, &
NBS, &
nDel, &
nFro, &
NOP, &
NOQ, &
NOR, &
nOrb, &
nOrbt, &
nOrbtt, &
NOS, &
NOVX, &
nSym, &
nTit, &
nTot, &
nTot1, &
nTot2, &
nTraToc, &
Occ, &
PotNuc, &
TcJobIph, &
TcOneMO, &
Title, &
VecTit

end module motra_global
