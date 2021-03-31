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

! ipCMO   : pointer to the MO coefficients
! nCMO    : size of Work(ipCMO)
! ipOcc   : pointer to occupancy vector
! ipEor   : pointer to orbital energy vector
! ipInd   : pointer to orbital type info vector
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

integer(kind=iwp) :: nBas(MxSym), nOrb(MxSym), nSym, ipCMO, nCMO, ipOcc, ipEor, ipMOrig, ipInd, LuSpool, nOrb2Loc(MxSym), &
                     nFro(MxSym), nAtoms, nMxIter, nOccInp(MxSym), nVirInp(MxSym), LocModel, nActa, iWave, MxConstr, nConstr(8)
real(kind=wp) :: Thrs, ThrRot, ThrGrad, ThrDomain(2), ThrPairDomain(3), ThrSel
logical(kind=iwp) :: LocNatOrb, LocCanOrb, Wave, Maximisation, ChoStart, DoCNOs, Silent, Test_Localisation, Analysis, PrintMOs, &
                     Timing, AnaAtom, EvalER, Order, LocPAO, AnaPAO, AnaPAO_Save, DoDomain, AnaDomain, Skip !, LocVir
character(len=LenIn8) :: BName(MxBas)
character(len=LenIn) :: NamAct(mxAtom)
character(len=512) :: LC_FileOrb
character(len=3) :: AnaNrm

public :: &
AnaAtom, &
AnaDomain, &
Analysis, &
AnaNrm, &
AnaPAO, &
AnaPAO_Save, &
ChoStart, &
DoCNOs, &
DoDomain, &
EvalER, &
ipCMO, &
ipEor, &
ipInd, &
ipMOrig, &
ipOcc, &
iWave, &
LC_FileOrb, &
LocCanOrb, &
LocModel, &
LocNatOrb, &
LocPAO, &
LuSpool, &
Maximisation, &
MxConstr, &
nActa, &
NamAct, &
BName, &
nAtoms, &
nBas, &
nCMO, &
nConstr, &
nFro, &
nMxIter, &
nOccInp, &
nOrb, &
nOrb2Loc, &
nSym, &
nVirInp, &
Order, &
PrintMOs, &
Silent, &
Skip, &
Test_Localisation, &
ThrDomain, &
ThrGrad, &
ThrPairDomain, &
ThrRot, &
Thrs, &
ThrSel, &
Timing, &
Wave

end module Localisation_globals

