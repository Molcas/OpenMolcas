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

subroutine DataSol(IDSolv)
! Database of optical and physical data for various solvent.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IDSolv
real(kind=wp), external :: DKList, RList
#include "rctfld.fh"

! Here begins the collection of data for each solvent.
if (IDSolv == 1) then

  ! Water

  Eps = 78.39_wp
  EpsInf = 1.776_wp
  DerEps = -0.3562_wp
  RSolv = 1.385_wp
  VMol = 18.07_wp
  TCE = 2.57e-4_wp
  !STen = 71.81_wp
  !DSTen = 0.650_wp
  !CMF = 1.277_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 3.348e-2_wp
  !NATyp = 2
  ! Atom O:
  NTT(1) = 1
  RDiff(1) = 1.5_wp
  KT(1) = DKList(8)
  RWT(1) = RList(8)
  ! Atom H:
  NTT(2) = 2
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 2) then

  ! Acetonitrile (Warning: RSolv and VMol for this solvent need to be checked)

  Eps = 36.64_wp
  EpsInf = 1.806_wp
  DerEps = Zero
  RSolv = 2.155_wp
  VMol = 53.68_wp
  TCE = 1.192e-3_wp
  !STen = 36.47_wp
  !DSTen = 1.373_wp
  !CMF = 0.808_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.153e-2_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 3
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6_wp
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 3) then

  ! Methanol

  Eps = 32.63_wp
  EpsInf = 1.758_wp
  DerEps = -0.1984_wp
  RSolv = 1.855_wp
  VMol = 40.7_wp
  TCE = 1.182e-3_wp
  !STen = 22.12_wp
  !DSTen = 1.154_wp
  !CMF = 1.776_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.495e-2_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 4
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5_wp
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 4) then

  ! Ethanol

  Eps = 24.55_wp
  EpsInf = 1.847_wp
  DerEps = -0.1510_wp
  RSolv = 2.180_wp
  VMol = 58.7_wp
  TCE = 1.103e-3_wp
  !STen = 21.89_wp
  !DSTen = 1.146_wp
  !CMF = 1.543_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.032e-2_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5_wp
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 5) then

  ! IsoQuinoline

  Eps = 10.43_wp
  EpsInf = 1.010_wp
  DerEps = Zero
  RSolv = 3.50_wp
  VMol = 117.27_wp
  TCE = 1.255e-3_wp
  !STen = 26.53_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.135e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 9
  RDiff(1) = 1.5_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 7
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6_wp
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 6) then

  ! Quinoline

  Eps = 9.03_wp
  EpsInf = 1.010_wp
  DerEps = Zero
  RSolv = 3.50_wp
  VMol = 117.27_wp
  TCE = 1.255e-3_wp
  ! STen = 26.53_wp
  ! DSTen = Zero
  ! CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  ! Rho = 5.135e-3_wp
  ! NATyp  = 3
  ! Atom C:
  NTT(1) = 9
  RDiff(1) = 1.5_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 7
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6_wp
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 7) then

  ! Chloroform

  Eps = 4.90_wp
  EpsInf = 2.085_wp
  DerEps = Zero
  RSolv = 2.48_wp
  VMol = 80.7_wp
  TCE = 1.255e-3_wp
  !STen = 26.53_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 7.482e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 2.82_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 1
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 3
  RDiff(3) = 1.79_wp
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 8) then

  ! EthylEther

  Eps = 4.335_wp
  EpsInf = Zero
  DerEps = Zero
  RSolv = 2.785_wp
  VMol = 103.84_wp
  TCE = 1.617e-3_wp
  !STen = Zero
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.799e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 4
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 10
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.50_wp
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 9) then

  ! MethyleneChloride

  Eps = 8.93_wp
  EpsInf = 2.020_wp
  DerEps = Zero
  RSolv = 2.27_wp
  VMol = 64.5_wp
  TCE = 1.367e-3_wp
  !STen = 27.33_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 9.74e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 2
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 2
  RDiff(3) = 1.79_wp
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 10) then

  ! DiChloroEthane

  Eps = 10.36_wp
  EpsInf = 2.085_wp
  DerEps = Zero
  RSolv = 2.505_wp
  VMol = 79.4_wp
  TCE = 1.156e-3_wp
  !STen = 31.54_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 7.576e-3_wp
  !NATyp  = 3
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 4
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 2
  RDiff(3) = 1.79_wp
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 11) then

  ! CarbonTetraChloride

  Eps = 2.228_wp
  EpsInf = 2.129_wp
  DerEps = Zero
  RSolv = 2.685_wp
  VMol = 96.5_wp
  TCE = 1.270e-3_wp
  !STen = 26.15_wp
  !DSTen = 1.436_wp
  !CMF = 0.629_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 6.241e-3_wp
  !NATyp = 2
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 2.82_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom Cl:
  NTT(2) = 4
  RDiff(2) = 1.79_wp
  KT(2) = DKList(17)
  RWT(2) = RList(17)

else if (IDSolv == 12) then

  ! Benzene

  Eps = 2.247_wp
  EpsInf = 2.244_wp
  DerEps = Zero
  RSolv = 2.63_wp
  VMol = 88.91_wp
  TCE = 1.380e-3_wp
  !STen = 28.18_wp
  !DSTen = 1.469_wp
  !CMF = 0.629_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 6.773e-3_wp
  !NATyp = 2
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.5_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 13) then

  ! Toluene

  Eps = 2.379_wp
  EpsInf = 2.232_wp
  DerEps = Zero
  RSolv = 2.82_wp
  VMol = 106.3_wp
  TCE = 1.08e-3_wp
  !STen = 27.92_wp
  !DSTen = 1.391_wp
  !CMF = 0.679_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.665e-3_wp
  !NATyp = 2
  ! Atom C:
  NTT(1) = 7
  RDiff(1) = 1.5_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 8
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 14) then

  ! ChloroBenzene

  Eps = 5.621_wp
  EpsInf = 2.320_wp
  DerEps = Zero
  RSolv = 2.805_wp
  VMol = 101.79_wp
  TCE = 0.981e-3_wp
  !STen = 32.69_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.916e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.5_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 5
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 1
  RDiff(3) = 1.79_wp
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 15) then

  ! NitroMethane

  Eps = 38.20_wp
  EpsInf = 1.904_wp
  DerEps = Zero
  RSolv = 2.155_wp
  VMol = 53.68_wp
  TCE = 1.192e-3_wp
  !STen = 36.47_wp
  !DSTen = 1.373_wp
  !CMF = 0.808_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.122e-2_wp
  !NATyp = 4
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 3
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6_wp
  KT(3) = DKList(7)
  RWT(3) = RList(7)
  ! Atom O:
  NTT(4) = 2
  RDiff(4) = 1.5_wp
  KT(4) = DKList(8)
  RWT(4) = RList(8)

else if (IDSolv == 16) then

  ! Heptane

  Eps = 1.92_wp
  EpsInf = 1.918_wp
  DerEps = Zero
  RSolv = 3.125_wp
  VMol = 146.56_wp
  TCE = 1.25e-3_wp
  !STen = 19.80_wp
  !DSTen = 1.505_wp
  !CMF = 0.687_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 4.109e-3_wp
  !NATyp = 2
  ! Atom C:
  NTT(1) = 7
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 16
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 17) then

  ! CycloHexane

  Eps = 2.023_wp
  EpsInf = 2.028_wp
  DerEps = Zero
  RSolv = 2.815_wp
  VMol = 108.10_wp
  TCE = 1.20e-3_wp
  !STen = 24.38_wp
  !DSTen = 1.467_wp
  !CMF = 0.621_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.571e-3_wp
  !NATyp = 2
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 12
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 18) then

  ! Aniline

  Eps = 6.89_wp
  EpsInf = 2.506_wp
  DerEps = Zero
  RSolv = 2.80_wp
  VMol = 91.15_wp
  TCE = 0.85e-3_wp
  !STen = 42.79_wp
  !DSTen = 0.731_wp
  !CMF = 0.972_wp
  ! Atomic parameters for dispersion and repulsion
  !Rho = 6.607e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.5_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 7
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6_wp
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 19) then

  ! Acetone

  Eps = 20.7_wp
  EpsInf = 1.841_wp
  DerEps = -9.77e-2_wp
  RSolv = 2.38_wp
  VMol = 73.52_wp
  TCE = 1.42e-3_wp
  !STen = 22.67_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 8.190e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 3
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5_wp
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 20) then

  ! TetraHydroFuran

  Eps = 7.58_wp
  EpsInf = 1.971_wp
  DerEps = Zero
  RSolv = 2.56_wp
  VMol = 81.11_wp
  TCE = 1.142e-3_wp
  !STen = 26.40_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 7.425e-3_wp
  !NATyp = 3
  ! Atom C:
  NTT(1) = 4
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 8
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5_wp
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 21) then

  ! DiMethylSulfoxide

  Eps = 46.7_wp
  EpsInf = 2.179_wp
  DerEps = -0.1902_wp
  RSolv = 2.455_wp
  VMol = 70.94_wp
  TCE = 9.82e-2_wp
  !STen = 42.86_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 8.49e-3_wp
  !NATyp = 4
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76_wp
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2_wp
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom S:
  NTT(3) = 1
  RDiff(3) = 1.80_wp
  KT(3) = DKList(16)
  RWT(3) = RList(16)
  ! Atom O:
  NTT(4) = 1
  RDiff(4) = 1.5_wp
  KT(4) = DKList(8)
  RWT(4) = RList(8)

else if (IDSolv == 22) then

  ! Argon

  ! Warning: the following data are referred to the absolute
  !          temperature of 118 K (the highest available). Change
  !          the parameter Tabs to this value in actual calculations.
  !          To use other temperatures, one has to provide also the
  !          corresponding density (in g cm-3) to compute RHO (i.e.
  !          the numeral density) and VMol (i.e. the molar volume)
  !          properly. Of course, the proper dielectric constant
  !          must be provided too.

  !Tabs = 118.0_wp
  Eps = 1.430_wp
  EpsInf = 1.430_wp
  DerEps = Zero
  RSolv = 1.875_wp
  VMol = 34.29_wp
  TCE = 9.82e-2_wp
  !STen = 42.86_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.756e-2_wp
  !NATyp = 1
  ! Atom Ar:
  NTT(1) = 1
  RDiff(1) = 1.875_wp
  KT(1) = DKList(18)
  RWT(1) = RList(18)

else if (IDSolv == 23) then

  ! Krypton

  ! Warning: the following data are referred to the absolute
  !          temperature of 168 K (the highest available). Change
  !          the parameter Tabs to this value in actual calculations.
  !          To use other temperatures, one has to provide also the
  !          corresponding density (in g cm-3) to compute RHO (i.e.
  !          the numeral density) and VMol (i.e. the molar volume)
  !          properly. Of course, the proper dielectric constant
  !          must be provided too.

  !Tabs = 168.0_wp
  Eps = 1.519_wp
  EpsInf = 1.519_wp
  DerEps = Zero
  RSolv = 2.07_wp
  VMol = 42.71_wp
  TCE = 9.82e-2_wp
  !STen = 42.86_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.410e-2_wp
  !NATyp = 1
  ! Atom Kr:
  NTT(1) = 1
  RDiff(1) = 2.07_wp
  KT(1) = DKList(36)
  RWT(1) = RList(36)

else if (IDSolv == 24) then

  ! Xenon

  ! Warning: the following data are referred to the absolute
  !          temperature of 210 K (the highest available). Change
  !          the parameter Tabs to this value in actual calculations.
  !          To use other temperatures, one has to provide also the
  !          corresponding density (in g cm-3) to compute RHO (i.e.
  !          the numeral density) and VMol (i.e. the molar volume)
  !          properly. Of course, the proper dielectric constant
  !          must be provided too.

  !Tabs = 210.0_wp
  Eps = 1.706_wp
  EpsInf = 1.706_wp
  DerEps = Zero
  RSolv = 2.20_wp
  VMol = 50.38_wp
  TCE = 9.82e-2_wp
  !STen = 42.86_wp
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.195e-2_wp
  !NATyp = 1
  ! Atom Xe:
  NTT(1) = 1
  RDiff(1) = 2.20_wp
  KT(1) = DKList(53)
  RWT(1) = RList(53)

end if

! Use user specified value of the dielectric constant

if (Eps_User /= -One) Eps = Eps_User
if (EpsInf_User /= Zero) EpsInf = EpsInf_User

return

end subroutine DataSol
