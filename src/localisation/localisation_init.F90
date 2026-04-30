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
                                LocCanOrb, Wave, iWave, DoCNOs, Loosen, ThrStep, AnalyseLoc, MoldMod, getIMmldn, useFH,&
                                inpOptMeth
use Definitions, only: iwp, wp
use Constants, only: Ten,Five,Half,One,Deg2Rad

implicit none
integer(kind=iwp) :: iSym
integer(kind=iwp), parameter :: Occupied = 0
real(kind=wp), parameter :: ThrsDef=1.0e-6_wp, & !functional change
                            ThrRotDef=1.0e-10_wp, & !rotation angle in jacobi sweeps
                            ThrGradDef=1.0e-4_wp,& !gradient norm
                            ThrStepDef=1.0e-2_wp !kappa norm

useFH=.true. !use full Pipek-Mezey Hessian in the SGEK
getIMmldn = .false.
MoldMod=-1 !at every MoldMod-th iteration, generate a molden file
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
inpOptMeth = 1 ! PM localisation done with Jacobi Sweeps
ChargeType = 1 ! PM localisation done within the Mulliken population framework
AnalyseLoc = 0
if (nSym > 1) LocModel = 3  ! Cholesky localisation
LocModel_UsrDef = .false.
Test_Localisation = .false.
NMxIter = 100
Thrs = ThrsDef
ThrRot = ThrRotDef
ThrGrad = ThrGradDef
ThrStep = ThrStepDef
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

! Default undershoot avoidance settings (for GEK)
Loosen%Thrs = cos(Five*Deg2Rad)
Loosen%Thrs2 = cos(20.0_wp*Deg2Rad)
Loosen%Step = Half*(One+sqrt(Five))
Loosen%Factor = One

end subroutine localisation_init

