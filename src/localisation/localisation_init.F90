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
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!               2026, Lila Zapp                                        *
!***********************************************************************

subroutine localisation_init()

use Localisation_globals, only: AnaAtom, AnaDomain, AnalyseLoc, Analysis, AnaNrm, AnaPAO, AnaPAO_Save, ChargeType, ChoStart, &
                                DoCNOs, DoDomain, EvalER, getIMmldn, inpOptMeth, iWave, LocCanOrb, LocModel, LocNatOrb, LocPAO, &
                                Maximisation, MoldMod, nConstr, nFro, NMxIter, nOrb2Loc, nSym, OptMeth, Order, PrintMOs, Skip, &
                                Test_Localisation, ThrDomain, ThrGrad, ThrPairDomain, ThrRot, Thrs, Timing, useFH, Wave
use Definitions, only: wp
use Constants, only: Ten

implicit none
real(kind=wp), parameter :: ThrsDef = 1.0e-6_wp, & !functional change
                            ThrRotDef = 1.0e-10_wp, & !rotation angle in jacobi sweeps
                            ThrGradDef = 1.0e-5_wp

useFH = .false. !use full Pipek-Mezey Hessian in the SGEK
getIMmldn = .false.
MoldMod = -1 !at every MoldMod-th iteration, generate a molden file
nOrb2Loc(1:nSym) = 0
nFro(1:nSym) = 0
nConstr(1:nSym) = 0
Skip = .false.
Maximisation = .true.
ChoStart = .false.
LocModel = 1 ! Pipek-Mezey localisation
OptMeth = 1 ! PM localisation done with Jacobi Sweeps
inpOptMeth = 1 ! PM localisation done with Jacobi Sweeps
ChargeType = 1 ! PM localisation done within the Mulliken population framework
AnalyseLoc = 0
if (nSym > 1) LocModel = 3  ! Cholesky localisation
Test_Localisation = .false.
NMxIter = 100
Thrs = ThrsDef
ThrRot = ThrRotDef
ThrGrad = ThrGradDef
Analysis = .false.
AnaAtom = nSym == 1
AnaNrm = 'Fro'
PrintMOs = .true.
Timing = .true.
EvalER = .false.
Order = .false.
LocPAO = .false.
AnaPAO = .false.
AnaPAO_Save = AnaPAO
DoDomain = .false.
AnaDomain = .false.
ThrDomain(1) = 0.9_wp
ThrDomain(2) = 2.0e-2_wp
ThrPairDomain(1) = 1.0e-10_wp
ThrPairDomain(2) = Ten
ThrPairDomain(3) = 15.0_wp
LocNatOrb = .false.
LocCanOrb = .false.
Wave = .false.
iWave = 0
DoCNOs = .false.

end subroutine localisation_init
