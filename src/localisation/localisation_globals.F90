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

module Localisation_globals

use Definitions, only: wp, iwp

implicit none
private

! CMO     : the MO coefficients
! nCMO    : size of CMO
! Occ     : occupancy vector
! Eor     : orbital energy vector
! ipInd   : orbital type info vector
! nSym    : number of Irrep
! nBas    : number of basis function per Irrep
! nOrb    : number of occ orb per irrep
! LuSpool : Unit number of the input
! nOrb2Loc: Number of orbital to localise
!
! BName   : Basis function names
!
! LC_FileOrb: orbital file for Seward to read

#include "Molcas.fh"

integer(kind=iwp) :: nBas(MxSym), nOrb(MxSym), nSym, nCMO, LuSpool, nOrb2Loc(MxSym), nFro(MxSym), nAtoms, nMxIter, nOccInp(MxSym), &
                     nVirInp(MxSym), LocModel, nActa, iWave, MxConstr, nConstr(MxSym)
real(kind=wp) :: Thrs, ThrRot, ThrGrad, ThrDomain(2), ThrPairDomain(3), ThrSel
logical(kind=iwp) :: LocNatOrb, LocCanOrb, Wave, Maximisation, ChoStart, DoCNOs, Silent, Test_Localisation, Analysis, PrintMOs, &
                     Timing, AnaAtom, EvalER, Order, LocPAO, AnaPAO, AnaPAO_Save, DoDomain, AnaDomain, Skip !, LocVir
character(len=512) :: LC_FileOrb
character(len=3) :: AnaNrm
integer(kind=iwp), allocatable :: Ind(:)
character(len=LenIn8), allocatable :: BName(:)
character(len=LenIn), allocatable :: NamAct(:)
real(kind=wp), allocatable :: CMO(:), EOrb(:), MOrig(:), Occ(:)

public :: AnaAtom, AnaDomain, Analysis, AnaNrm, AnaPAO, AnaPAO_Save, BName, ChoStart, CMO, DoCNOs, DoDomain, EOrb, EvalER, Ind, &
          iWave, LC_FileOrb, LocCanOrb, LocModel, LocNatOrb, LocPAO, LuSpool, Maximisation, MOrig, MxConstr, nActa, NamAct, &
          nAtoms, nBas, nCMO, nConstr, nFro, nMxIter, nOccInp, nOrb, nOrb2Loc, nSym, nVirInp, Occ, Order, PrintMOs, Silent, Skip, &
          Test_Localisation, ThrDomain, ThrGrad, ThrPairDomain, ThrRot, Thrs, ThrSel, Timing, Wave

end module Localisation_globals

