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
!               2026, Lila Zapp (split subroutine for initialisation)  *
!***********************************************************************

subroutine localisation_init()
use Localisation_globals, only: nSym, nOrb2Loc, nFro, nConstr, Skip, LocOrb, Thrs_UsrDef, nOrb2Loc_UsrDef, nFro_UsrDef, Freeze,&
                                Maximisation, ChoStart, LocModel, OptMeth, ChargeType, LocModel_UsrDef,Test_Localisation, &
                                NMxIter, Thrs, ThrRot, ThrGrad, Analysis, AnaAtom, AnaNrm, PrintMOs, Timing, EvalER, Order,&
                                LocPAO, AnaPAO, AnaPAO_Save, DoDomain, AnaDomain, ThrDomain, ThrPairDomain, LocNatOrb, &
                                LocCanOrb, Wave, iWave, DoCNOs
use definitions, only: iwp, wp
use constants, only: Ten

implicit none
integer(kind=iwp) :: iSym
integer(kind=iwp), parameter :: Occupied = 0
real(kind=wp), parameter :: ThrsDef=1.0e-6_wp, ThrRotDef=1.0e-10_wp, ThrGradDef=1.0e-2_wp

do iSym=1,nSym
  nOrb2Loc(iSym) = 0
  nFro(iSym) = 0
  nConstr(iSym) = 0
end do
Skip = .false.
LocOrb = Occupied
Thrs_UsrDef = .false.
nOrb2Loc_UsrDef = .false.
nFro_UsrDef = .false.
Freeze = .false.
Maximisation = .true.
ChoStart = .false.
LocModel = 1  ! Pipek-Mezey localisation
OptMeth = 1 ! PM localisation done with Jacobi Sweeps
ChargeType = 1 ! PM localisation done within the Mulliken population framework
if (nSym > 1) LocModel = 3  ! Cholesky localisation
LocModel_UsrDef = .false.
Test_Localisation = .false.
NMxIter = 300
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

