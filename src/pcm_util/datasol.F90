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

implicit real*8(A-H,O-Z)
#include "rctfld.fh"
data One/1.0d0/,Zero/0.0d0/

! Here begins the collection of data for each solvent.
if (IDSolv == 1) then

  ! Water

  Eps = 78.39d0
  EpsInf = 1.776d0
  DerEps = -0.3562d0
  RSolv = 1.385d0
  VMol = 18.07d0
  TCE = 2.57d-4
  !STen = 71.81d0
  !DSTen = 0.650d0
  !CMF = 1.277d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 3.348d-2
  !NATyp = 2
  ! Atom O:
  NTT(1) = 1
  RDiff(1) = 1.5d0
  KT(1) = DKList(8)
  RWT(1) = RList(8)
  ! Atom H:
  NTT(2) = 2
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 2) then

  ! Acetonitrile (Warning: RSolv and VMol for this solvent need to be checked)

  Eps = 36.64d0
  EpsInf = 1.806d0
  DerEps = Zero
  RSolv = 2.155d0
  VMol = 53.68d0
  TCE = 1.192d-3
  !STen = 36.47d0
  !DSTen = 1.373d0
  !CMF = 0.808d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.153d-2
  !NATyp = 3
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 3
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6d0
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 3) then

  ! Methanol

  Eps = 32.63d0
  EpsInf = 1.758d0
  DerEps = -0.1984d0
  RSolv = 1.855d0
  VMol = 40.7d0
  TCE = 1.182d-3
  !STen = 22.12d0
  !DSTen = 1.154d0
  !CMF = 1.776d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.495d-2
  !NATyp = 3
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 4
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5d0
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 4) then

  ! Ethanol

  Eps = 24.55d0
  EpsInf = 1.847d0
  DerEps = -0.1510d0
  RSolv = 2.180d0
  VMol = 58.7d0
  TCE = 1.103d-3
  !STen = 21.89d0
  !DSTen = 1.146d0
  !CMF = 1.543d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.032d-2
  !NATyp = 3
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5d0
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 5) then

  ! IsoQuinoline

  Eps = 10.43d0
  EpsInf = 1.010d0
  DerEps = Zero
  RSolv = 3.50d0
  VMol = 117.27d0
  TCE = 1.255d-3
  !STen = 26.53d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.135d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 9
  RDiff(1) = 1.5d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 7
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6d0
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 6) then

  ! Quinoline

  Eps = 9.03d0
  EpsInf = 1.010d0
  DerEps = Zero
  RSolv = 3.50d0
  VMol = 117.27d0
  TCE = 1.255d-3
  ! STen = 26.53d0
  ! DSTen = Zero
  ! CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  ! Rho = 5.135d-3
  ! NATyp  = 3
  ! Atom C:
  NTT(1) = 9
  RDiff(1) = 1.5d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 7
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6d0
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 7) then

  ! Chloroform

  Eps = 4.90d0
  EpsInf = 2.085d0
  DerEps = Zero
  RSolv = 2.48d0
  VMol = 80.7d0
  TCE = 1.255d-3
  !STen = 26.53d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 7.482d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 2.82d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 1
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 3
  RDiff(3) = 1.79d0
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 8) then

  ! EthylEther

  Eps = 4.335d0
  EpsInf = Zero
  DerEps = Zero
  RSolv = 2.785d0
  VMol = 103.84d0
  TCE = 1.617d-3
  !STen = Zero
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.799d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 4
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 10
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.50d0
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 9) then

  ! MethyleneChloride

  Eps = 8.93d0
  EpsInf = 2.020d0
  DerEps = Zero
  RSolv = 2.27d0
  VMol = 64.5d0
  TCE = 1.367d-3
  !STen = 27.33d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 9.74d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 2
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 2
  RDiff(3) = 1.79d0
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 10) then

  ! DiChloroEthane

  Eps = 10.36d0
  EpsInf = 2.085d0
  DerEps = Zero
  RSolv = 2.505d0
  VMol = 79.4d0
  TCE = 1.156d-3
  !STen = 31.54d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 7.576d-3
  !NATyp  = 3
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 4
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 2
  RDiff(3) = 1.79d0
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 11) then

  ! CarbonTetraChloride

  Eps = 2.228d0
  EpsInf = 2.129d0
  DerEps = Zero
  RSolv = 2.685d0
  VMol = 96.5d0
  TCE = 1.270d-3
  !STen = 26.15d0
  !DSTen = 1.436d0
  !CMF = 0.629d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 6.241d-3
  !NATyp = 2
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 2.82d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom Cl:
  NTT(2) = 4
  RDiff(2) = 1.79d0
  KT(2) = DKList(17)
  RWT(2) = RList(17)

else if (IDSolv == 12) then

  ! Benzene

  Eps = 2.247d0
  EpsInf = 2.244d0
  DerEps = Zero
  RSolv = 2.63d0
  VMol = 88.91d0
  TCE = 1.380d-3
  !STen = 28.18d0
  !DSTen = 1.469d0
  !CMF = 0.629d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 6.773d-3
  !NATyp = 2
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.5d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 13) then

  ! Toluene

  Eps = 2.379d0
  EpsInf = 2.232d0
  DerEps = Zero
  RSolv = 2.82d0
  VMol = 106.3d0
  TCE = 1.08d-3
  !STen = 27.92d0
  !DSTen = 1.391d0
  !CMF = 0.679d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.665d-3
  !NATyp = 2
  ! Atom C:
  NTT(1) = 7
  RDiff(1) = 1.5d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 8
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 14) then

  ! ChloroBenzene

  Eps = 5.621d0
  EpsInf = 2.320d0
  DerEps = Zero
  RSolv = 2.805d0
  VMol = 101.79d0
  TCE = 0.981d-3
  !STen = 32.69d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.916d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.5d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 5
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom Cl:
  NTT(3) = 1
  RDiff(3) = 1.79d0
  KT(3) = DKList(17)
  RWT(3) = RList(17)

else if (IDSolv == 15) then

  ! NitroMethane

  Eps = 38.20d0
  EpsInf = 1.904d0
  DerEps = Zero
  RSolv = 2.155d0
  VMol = 53.68d0
  TCE = 1.192d-3
  !STen = 36.47d0
  !DSTen = 1.373d0
  !CMF = 0.808d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.122d-2
  !NATyp = 4
  ! Atom C:
  NTT(1) = 1
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 3
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6d0
  KT(3) = DKList(7)
  RWT(3) = RList(7)
  ! Atom O:
  NTT(4) = 2
  RDiff(4) = 1.5d0
  KT(4) = DKList(8)
  RWT(4) = RList(8)

else if (IDSolv == 16) then

  ! Heptane

  Eps = 1.92d0
  EpsInf = 1.918d0
  DerEps = Zero
  RSolv = 3.125d0
  VMol = 146.56d0
  TCE = 1.25d-3
  !STen = 19.80d0
  !DSTen = 1.505d0
  !CMF = 0.687d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 4.109d-3
  !NATyp = 2
  ! Atom C:
  NTT(1) = 7
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 16
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 17) then

  ! CycloHexane

  Eps = 2.023d0
  EpsInf = 2.028d0
  DerEps = Zero
  RSolv = 2.815d0
  VMol = 108.10d0
  TCE = 1.20d-3
  !STen = 24.38d0
  !DSTen = 1.467d0
  !CMF = 0.621d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 5.571d-3
  !NATyp = 2
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 12
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)

else if (IDSolv == 18) then

  ! Aniline

  Eps = 6.89d0
  EpsInf = 2.506d0
  DerEps = Zero
  RSolv = 2.80d0
  VMol = 91.15d0
  TCE = 0.85d-3
  !STen = 42.79d0
  !DSTen = 0.731d0
  !CMF = 0.972d0
  ! Atomic parameters for dispersion and repulsion
  !Rho = 6.607d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 6
  RDiff(1) = 1.5d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 7
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom N:
  NTT(3) = 1
  RDiff(3) = 1.6d0
  KT(3) = DKList(7)
  RWT(3) = RList(7)

else if (IDSolv == 19) then

  ! Acetone

  Eps = 20.7d0
  EpsInf = 1.841d0
  DerEps = -9.77d-2
  RSolv = 2.38d0
  VMol = 73.52d0
  TCE = 1.42d-3
  !STen = 22.67d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 8.190d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 3
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5d0
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 20) then

  ! TetraHydroFuran

  Eps = 7.58d0
  EpsInf = 1.971d0
  DerEps = Zero
  RSolv = 2.56d0
  VMol = 81.11d0
  TCE = 1.142d-3
  !STen = 26.40d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 7.425d-3
  !NATyp = 3
  ! Atom C:
  NTT(1) = 4
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 8
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom O:
  NTT(3) = 1
  RDiff(3) = 1.5d0
  KT(3) = DKList(8)
  RWT(3) = RList(8)

else if (IDSolv == 21) then

  ! DiMethylSulfoxide

  Eps = 46.7d0
  EpsInf = 2.179d0
  DerEps = -0.1902d0
  RSolv = 2.455d0
  VMol = 70.94d0
  TCE = 9.82d-2
  !STen = 42.86d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 8.49d-3
  !NATyp = 4
  ! Atom C:
  NTT(1) = 2
  RDiff(1) = 1.76d0
  KT(1) = DKList(6)
  RWT(1) = RList(6)
  ! Atom H:
  NTT(2) = 6
  RDiff(2) = 1.2d0
  KT(2) = DKList(1)
  RWT(2) = RList(1)
  ! Atom S:
  NTT(3) = 1
  RDiff(3) = 1.80d0
  KT(3) = DKList(16)
  RWT(3) = RList(16)
  ! Atom O:
  NTT(4) = 1
  RDiff(4) = 1.5d0
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

  !Tabs = 118.0d0
  Eps = 1.430d0
  EpsInf = 1.430d0
  DerEps = Zero
  RSolv = 1.875d0
  VMol = 34.29d0
  TCE = 9.82d-2
  !STen = 42.86d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.756d-2
  !NATyp = 1
  ! Atom Ar:
  NTT(1) = 1
  RDiff(1) = 1.875d0
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

  !Tabs = 168.0d0
  Eps = 1.519d0
  EpsInf = 1.519d0
  DerEps = Zero
  RSolv = 2.07d0
  VMol = 42.71d0
  TCE = 9.82d-2
  !STen = 42.86d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.410d-02
  !NATyp = 1
  ! Atom Kr:
  NTT(1) = 1
  RDiff(1) = 2.07d0
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

  !Tabs = 210.0d0
  Eps = 1.706d0
  EpsInf = 1.706d0
  DerEps = Zero
  RSolv = 2.20d0
  VMol = 50.38d0
  TCE = 9.82d-2
  !STen = 42.86d0
  !DSTen = Zero
  !CMF = Zero
  ! Atomic parameters for dispersion and repulsion
  !Rho = 1.195d-02
  !NATyp = 1
  ! Atom Xe:
  NTT(1) = 1
  RDiff(1) = 2.20d0
  KT(1) = DKList(53)
  RWT(1) = RList(53)

end if

! Use user specified value of the dielectric constant

if (Eps_User /= -One) Eps = Eps_User
if (EpsInf_User /= Zero) EpsInf = EpsInf_User

return

end subroutine DataSol
!====
function RList(IA)
! Assigns Caillet-Claverie's atomic parameters for dispersion and repulsion.

implicit real*8(A-H,O-Z)
parameter(MaxIA=110)
real*8 RW(MaxIA)
save RW
data RW/ &
        !  H       He
        1.20d+00,1.28d+00, &
        ! Li       Be        B        C        N        O        F       Ne
        0.00d+00,0.00d+00,0.00d+00,1.70d+00,1.60d+00,1.50d+00,1.45d+00,1.38d+00, &
        ! Na       Mg       Al       Si        P        S       Cl       Ar
        1.20d+00,0.00d+00,0.00d+00,0.00d+00,1.85d+00,1.80d+00,1.76d+00,1.66d+00, &
        !  K      Ca
        1.46d+00,0.00d+00, &
                 ! Sc       Ti        V       Cr       Mn       Fe       Co       Ni       Cu       Zn
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                          ! Ga       Ge       As       Se       Br       Kr
                          0.00d+00,0.00d+00,0.00d+00,0.00d+00,1.85d+00,1.76d+00, &
        ! Rb       Sr
        0.00d+00,0.00d+00, &
                 !  Y       Zr       Nb       Mo       Tc       Ru       Rh       Pd       Ag       Cd
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                          ! In       Sn       Sb       Te        I       Xe
                          0.00d+00,0.00d+00,0.00d+00,0.00d+00,1.96d+00,1.85d+00, &
        ! Cs       Ba
        0.00d+00,0.00d+00, &
                 ! La       Ce       Pr       Nd       Pm       Sm       Eu
                 ! Gd       Tb       Dy       Ho       Er       Tm       Yb
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                 ! Lu       Hf       Ta        W       Re       Os       Ir       Pt       Au       Hg
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                          ! Tl       Pb       Bi       Po       At       Rn
                          0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
        ! Fr       Ra
        0.00d+00,0.00d+00, &
                 ! Ac       Th       Pa        U       Np       Pu       Am
                 ! Cm       Bk       Cf       Es       Fm       Md
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                 0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,9*0.00d+00/

if ((IA > 0) .and. (IA <= MaxIA)) then
  RList = RW(IA)
else
  RList = 0.0d0
  write(6,'(a)') 'IA out of range in RList.'
  call Abend()
end if

return

end function RList
!====
function DKList(NumAt)
! Assigns Caillet-Claverie's atomic parameters for dispersion and repulsion.

implicit real*8(A-H,O-Z)
parameter(NElem=110)
real*8 K(NElem)
save K
data K/ &
       !  H        He
       1.000d+00,0.594d+00, &
       ! Li        Be         B         C         N        O          F        Ne
       0.000d+00,0.000d+00,0.000d+00,1.000d+00,1.180d+00,1.360d+00,1.500d+00,1.209d+00, &
       ! Na        Mg        Al        Si         P        S         Cl        Ar
       1.400d+00,0.000d+00,0.000d+00,0.000d+00,2.100d+00,2.400d+00,2.100d+00,2.127d+00, &
       !  K        Ca
       2.900d+00,0.000d+00, &
                 ! Sc        Ti         V        Cr        Mn        Fe        Co        Ni        Cu        Zn
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
                           ! Ga        Ge        As        Se        Br        Kr
                           0.000d+00,0.000d+00,0.000d+00,0.000d+00,2.400d+00,2.566d+00, &
       ! Rb        Sr
       0.000d+00,0.000d+00, &
                 !  Y        Zr        Nb        Mo        Tc        Ru        Rh        Pd        Ag        Cd
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
                           ! In        Sn        Sb        Te         I        Xe
                           0.000d+00,0.000d+00,0.000d+00,0.000d+00,3.200d+00,0.000d+00, &
       ! Cs        Ba
       0.000d+00,0.000d+00, &
                 ! La        Ce        Pr        Nd        Pm        Sm        Eu
                 ! Gd        Tb        Dy        Ho        Er        Tm        Yb
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
                 ! Lu        Hf        Ta         W        Re        Os        Ir        Pt        Au        Hg
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
                           ! Tl        Pb        Bi        Po        At        Rn
                           0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
       ! Fr        Ra
       0.000d+00,0.000d+00, &
                 ! Ac        Th        Pa         U        Np        Pu        Am
                 ! Cm        Bk        Cf        Es        Fm        Md
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00, &
                 0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,0.000d+00,9*0.00d+00/

DKList = K(NumAt)

return

end function DKList
