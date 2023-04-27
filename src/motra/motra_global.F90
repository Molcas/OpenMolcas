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

#include "Molcas.fh"
integer(kind=iwp), parameter :: MxTit = 1
!----------------------------------------------------------------------*
! Allocate space to store the system description                       *
!----------------------------------------------------------------------*
integer(kind=iwp) :: nBas(mxSym), nDel(mxSym), nFro(mxSym), nOrb(mxSym), nSym
real(kind=wp) :: PotNuc
character(len=LenIn8), allocatable :: BsLbl(:)

!----------------------------------------------------------------------*
! Allocate space to store the MO-coefficients and occupations          *
!----------------------------------------------------------------------*
real(kind=wp), allocatable :: CMO(:), HOne(:), Kine(:), Occ(:), Ovlp(:)
! Occ is never allocated?

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
! TwoAO:          two-electron integrals in AO basis
! OneMO, TwoMO:   one- and two-electron integrals in MO basis
! Half:           half transformed two-electron integrals
integer(kind=iwp) :: LuHalf, LuInpOrb, LuJobIph, LuOneMO, LuTwoAO, LuTwoMO
character(len=8) :: FnHalf, FnOneMO, FnTwoAO, FnTwoMO
character(len=180) :: FnInpOrb, FnJobIph

integer(kind=iwp) :: iCTonly, iDoInt
integer(kind=iwp) :: nTot1, nTot2, n2max, nOrbt, nOrbtt, ISP, ISQ, ISR, ISS, NBP, NBQ, NBR, NBS, NOP, NOQ, NOR, NOS, LMOP, LMOQ, &
                     LMOR, LMOS, NBPQ, NBRS, NOVX, LTUVX, IAD13

!----------------------------------------------------------------------*
! Option to skip orthogonalization in motra, required for GronOR runs  *
!----------------------------------------------------------------------*
integer(kind=iwp) :: iortho

public :: BsLbl, CMO, CutThrs, Debug, FnHalf, FnInpOrb, FnJobIph, FnOneMO, FnTwoAO, FnTwoMO, Header, HOne, IAD13, iAutoCut, &
          iCTonly, iDoInt, ihdf5, iOneOnly, iortho, iPrint, iRFpert, ISP, ISQ, ISR, ISS, iVecTyp, LMOP, Kine, LMOQ, LMOR, LMOS, &
          LTUVX, LuHalf, LuInpOrb, LuJobIph, LuOneMO, LuTwoAO, LuTwoMO, MxTit, N2MAX, nBas, NBP, NBPQ, NBQ, NBR, NBRS, NBS, nDel, &
          nFro, NOP, NOQ, NOR, nOrb, nOrbt, nOrbtt, NOS, NOVX, nSym, nTit, nTot1, nTot2, Occ, Ovlp, PotNuc, Title, VecTit

end module motra_global
