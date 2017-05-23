************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine DataSol(IDSolv)
C
C     Database of optical and physical data for various solvent.
C
      Implicit Real*8 (A-H,O-Z)
#include "rctfld.fh"
      Data One/1.0D0/, Zero/0.0d0/
C
C     Here begins the collection of data for each solvent.
      If(IDSolv.eq. 1) Then
C
C     Water
C
        Eps    = 78.39d0
        EpsInf = 1.776d0
        DerEps = -0.3562d0
        RSolv  = 1.385d0
        VMol   = 18.07d0
        TCE    = 2.57d-4
        STen   = 71.81d0
        DSTen  = 0.650d0
        CMF    = 1.277d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 3.348d-2
        NATyp  = 2
C     Atom O:
        NTT(1)    = 1
        RDiff(1) = 1.5d0
        KT(1)    = DKList(8)
        RWT(1)   = RList(8)
C     Atom H:
        NTT(2)    = 2
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C
      Else If(IDSolv.eq. 2) Then
C
C     Acetonitrile (Warning: RSolv and VMol for this solvent
C                   need to be checked)
C
        Eps    = 36.64d0
        EpsInf = 1.806d0
        DerEps = Zero
        RSolv  = 2.155d0
        VMol   = 53.68d0
        TCE    = 1.192d-3
        STen   = 36.47d0
        DSTen  = 1.373d0
        CMF    = 0.808d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.153d-2
        NATyp  = 3
C     Atom C:
        NTT(1)    = 2
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 3
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom N:
        NTT(3)    = 1
        RDiff(3) = 1.6d0
        KT(3)    = DKList(7)
        RWT(3)   = RList(7)
C
      Else If(IDSolv.eq. 3) Then
C
C     Methanol
C
        Eps    = 32.63d0
        EpsInf = 1.758d0
        DerEps = -0.1984d0
        RSolv  = 1.855d0
        VMol   = 40.7d0
        TCE    = 1.182d-3
        STen   = 22.12d0
        DSTen  = 1.154d0
        CMF    = 1.776d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.495d-2
        NATyp  = 3
C     Atom C:
        NTT(1)    = 1
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 4
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom O:
        NTT(3)    = 1
        RDiff(3) = 1.5d0
        KT(3)    = DKList(8)
        RWT(3)   = RList(8)
C
      Else If(IDSolv.eq. 4) Then
C
C     Ethanol
C
        Eps    = 24.55d0
        EpsInf = 1.847d0
        DerEps = -0.1510d0
        RSolv  = 2.180d0
        VMol   = 58.7d0
        TCE    = 1.103d-3
        STen   = 21.89d0
        DSTen  = 1.146d0
        CMF    = 1.543d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.032d-2
        NATyp  = 3
C     Atom C:
        NTT(1)    = 2
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 6
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom O:
        NTT(3)    = 1
        RDiff(3) = 1.5d0
        KT(3)    = DKList(8)
        RWT(3)   = RList(8)
C
      Else If(IDSolv.eq. 5) Then
C
C     IsoQuinoline
C
        Eps    = 10.43d0
        EpsInf = 1.010d0
        DerEps = Zero
        RSolv  = 3.50d0
        VMol   = 117.27d0
        TCE    = 1.255d-3
        STen   = 26.53d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 5.135d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 9
        RDiff(1) = 1.5d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 7
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom N:
        NTT(3)    = 1
        RDiff(3) = 1.6d0
        KT(3)    = DKList(7)
        RWT(3)   = RList(7)
C
      Else If(IDSolv.eq. 6) Then
C
C     Quinoline
C
        Eps    = 9.03d0
        EpsInf = 1.010d0
        DerEps = Zero
        RSolv  = 3.50d0
        VMol   = 117.27d0
        TCE    = 1.255d-3
        STen   = 26.53d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 5.135d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 9
        RDiff(1) = 1.5d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 7
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom N:
        NTT(3)    = 1
        RDiff(3) = 1.6d0
        KT(3)    = DKList(7)
        RWT(3)   = RList(7)
C
      Else If(IDSolv.eq. 7) Then
C
C     Chloroform
C
        Eps    = 4.90d0
        EpsInf = 2.085d0
        DerEps = Zero
        RSolv  = 2.48d0
        VMol   = 80.7d0
        TCE    = 1.255d-3
        STen   = 26.53d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 7.482d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 1
        RDiff(1) = 2.82d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 1
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom Cl:
        NTT(3)    = 3
        RDiff(3) = 1.79d0
        KT(3)    = DKList(17)
        RWT(3)   = RList(17)
C
      Else If(IDSolv.eq. 8) Then
C
C     EthylEther
C
        Eps    = 4.335d0
        EpsInf = Zero
        DerEps = Zero
        RSolv  = 2.785d0
        VMol   = 103.84d0
        TCE    = 1.617d-3
        STen   = Zero
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 5.799d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 4
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 10
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom O:
        NTT(3)    = 1
        RDiff(3) = 1.50d0
        KT(3)    = DKList(8)
        RWT(3)   = RList(8)
C
      Else If(IDSolv.eq. 9) Then
C
C     MethyleneChloride
C
        Eps    = 8.93d0
        EpsInf = 2.020d0
        DerEps = Zero
        RSolv  = 2.27d0
        VMol   = 64.5d0
        TCE    = 1.367d-3
        STen   = 27.33d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 9.74d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 1
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 2
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom Cl:
        NTT(3)    = 2
        RDiff(3) = 1.79d0
        KT(3)    = DKList(17)
        RWT(3)   = RList(17)
C
      Else If(IDSolv.eq.10) Then
C
C     DiChloroEthane
C
        Eps    = 10.36d0
        EpsInf = 2.085d0
        DerEps = Zero
        RSolv  = 2.505d0
        VMol   = 79.4d0
        TCE    = 1.156d-3
        STen   = 31.54d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 7.576d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 2
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 4
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom Cl:
        NTT(3)    = 2
        RDiff(3) = 1.79d0
        KT(3)    = DKList(17)
        RWT(3)   = RList(17)
C
      Else If(IDSolv.eq.11) Then
C
C     CarbonTetraChloride
C
        Eps    = 2.228d0
        EpsInf = 2.129d0
        DerEps = Zero
        RSolv  = 2.685d0
        VMol   = 96.5d0
        TCE    = 1.270d-3
        STen   = 26.15d0
        DSTen  = 1.436d0
        CMF    = 0.629d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 6.241d-3
        NATyp  = 2
C     Atom C:
        NTT(1)    = 1
        RDiff(1) = 2.82d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom Cl:
        NTT(2)    = 4
        RDiff(2) = 1.79d0
        KT(2)    = DKList(17)
        RWT(2)   = RList(17)
C
      Else If(IDSolv.eq.12) Then
C
C     Benzene
C
        Eps    = 2.247d0
        EpsInf = 2.244d0
        DerEps = Zero
        RSolv  = 2.63d0
        VMol   = 88.91d0
        TCE    = 1.380d-3
        STen   = 28.18d0
        DSTen  = 1.469d0
        CMF    = 0.629d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 6.773d-3
        NATyp  = 2
C     Atom C:
        NTT(1)    = 6
        RDiff(1) = 1.5d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 6
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C
      Else If(IDSolv.eq.13) Then
C
C     Toluene
C
        Eps    = 2.379d0
        EpsInf = 2.232d0
        DerEps = Zero
        RSolv  = 2.82d0
        VMol   = 106.3d0
        TCE    = 1.08d-3
        STen   = 27.92d0
        DSTen  = 1.391d0
        CMF    = 0.679d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 5.665d-3
        NATyp  = 2
C     Atom C:
        NTT(1)    = 7
        RDiff(1) = 1.5d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 8
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C
      Else If(IDSolv.eq.14) Then
C
C     ChloroBenzene
C
        Eps    = 5.621d0
        EpsInf = 2.320d0
        DerEps = Zero
        RSolv  = 2.805d0
        VMol   = 101.79d0
        TCE    = 0.981d-3
        STen   = 32.69d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 5.916d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 6
        RDiff(1) = 1.5d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 5
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom Cl:
        NTT(3)    = 1
        RDiff(3) = 1.79d0
        KT(3)    = DKList(17)
        RWT(3)   = RList(17)
C
      Else If(IDSolv.eq.15) Then
C
C     NitroMethane
C
        Eps    = 38.20d0
        EpsInf = 1.904d0
        DerEps = Zero
        RSolv  = 2.155d0
        VMol   = 53.68d0
        TCE    = 1.192d-3
        STen   = 36.47d0
        DSTen  = 1.373d0
        CMF    = 0.808d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.122d-2
        NATyp  = 4
C     Atom C:
        NTT(1)    = 1
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 3
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom N:
        NTT(3)    = 1
        RDiff(3) = 1.6d0
        KT(3)    = DKList(7)
        RWT(3)   = RList(7)
C     Atom O:
        NTT(4)    = 2
        RDiff(4) = 1.5d0
        KT(4)    = DKList(8)
        RWT(4)   = RList(8)
C
      Else If(IDSolv.eq.16) Then
C
C     Heptane
C
        Eps    = 1.92d0
        EpsInf = 1.918d0
        DerEps = Zero
        RSolv  = 3.125d0
        VMol   = 146.56d0
        TCE    = 1.25d-3
        STen   = 19.80d0
        DSTen  = 1.505d0
        CMF    = 0.687d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 4.109d-3
        NATyp  = 2
C     Atom C:
        NTT(1)    = 7
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 16
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C
      Else If(IDSolv.eq.17) Then
C
C     CycloHexane
C
        Eps    = 2.023d0
        EpsInf = 2.028d0
        DerEps = Zero
        RSolv  = 2.815d0
        VMol   = 108.10d0
        TCE    = 1.20d-3
        STen   = 24.38d0
        DSTen  = 1.467d0
        CMF    = 0.621d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 5.571d-3
        NATyp  = 2
C     Atom C:
        NTT(1)    = 6
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 12
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C
      Else If(IDSolv.eq.18) Then
C
C     Aniline
C
        Eps    = 6.89d0
        EpsInf = 2.506d0
        DerEps = Zero
        RSolv  = 2.80d0
        VMol   = 91.15d0
        TCE    = 0.85d-3
        STen   = 42.79d0
        DSTen  = 0.731d0
        CMF    = 0.972d0
C     Atomic parameters for dispersion and repulsion
        Rho    = 6.607d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 6
        RDiff(1) = 1.5d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 7
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom N:
        NTT(3)    = 1
        RDiff(3) = 1.6d0
        KT(3)    = DKList(7)
        RWT(3)   = RList(7)
C
      Else If(IDSolv.eq.19) Then
C
C     Acetone
C
        Eps    = 20.7d0
        EpsInf = 1.841d0
        DerEps = -9.77d-2
        RSolv  = 2.38d0
        VMol   = 73.52d0
        TCE    = 1.42d-3
        STen   = 22.67d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 8.190d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 3
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 6
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom O:
        NTT(3)    = 1
        RDiff(3) = 1.5d0
        KT(3)    = DKList(8)
        RWT(3)   = RList(8)
C
      Else If(IDSolv.eq.20) Then
C
C     TetraHydroFuran
C
        Eps    = 7.58d0
        EpsInf = 1.971d0
        DerEps = Zero
        RSolv  = 2.56d0
        VMol   = 81.11d0
        TCE    = 1.142d-3
        STen   = 26.40d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 7.425d-3
        NATyp  = 3
C     Atom C:
        NTT(1)    = 4
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 8
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom O:
        NTT(3)    = 1
        RDiff(3) = 1.5d0
        KT(3)    = DKList(8)
        RWT(3)   = RList(8)
C
      Else If(IDSolv.eq.21) Then
C
C     DiMethylSulfoxide
C
        Eps    = 46.7d0
        EpsInf = 2.179d0
        DerEps = -0.1902d0
        RSolv  = 2.455d0
        VMol   = 70.94d0
        TCE    = 9.82d-2
        STen   = 42.86d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 8.49d-3
        NATyp  = 4
C     Atom C:
        NTT(1)    = 2
        RDiff(1) = 1.76d0
        KT(1)    = DKList(6)
        RWT(1)   = RList(6)
C     Atom H:
        NTT(2)    = 6
        RDiff(2) = 1.2d0
        KT(2)    = DKList(1)
        RWT(2)   = RList(1)
C     Atom S:
        NTT(3)    = 1
        RDiff(3) = 1.80d0
        KT(3)    = DKList(16)
        RWT(3)   = RList(16)
C     Atom O:
        NTT(4)    = 1
        RDiff(4) = 1.5d0
        KT(4)    = DKList(8)
        RWT(4)   = RList(8)
C
      Else If(IDSolv.eq.22) Then
C
C     Argon
C
C     Warning: the following data are referred to the absolute
C              temperature of 118 K (the highest available). Change
C              the parameter Tabs to this value in actual calculations.
C              To use other temperatures, one has to provide also the
C              corresponding density (in g cm-3) to compute RHO (i.e.
C              the numeral density) and VMol (i.e. the molar volume)
C              properly. Of course, the proper dielectric constant
C              must be provided too.
C
        Tabs   = 118.0d0
        Eps    = 1.430d0
        EpsInf = 1.430d0
        DerEps = Zero
        RSolv  = 1.875d0
        VMol   = 34.29d0
        TCE    = 9.82d-2
        STen   = 42.86d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.756d-2
        NATyp  = 1
C     Atom Ar:
        NTT(1)    = 1
        RDiff(1) = 1.875d0
        KT(1)    = DKList(18)
        RWT(1)   = RList(18)
C
      Else If(IDSolv.eq.23) Then
C
C     Krypton
C
C     Warning: the following data are referred to the absolute
C              temperature of 168 K (the highest available). Change
C              the parameter Tabs to this value in actual calculations.
C              To use other temperatures, one has to provide also the
C              corresponding density (in g cm-3) to compute RHO (i.e.
C              the numeral density) and VMol (i.e. the molar volume)
C              properly. Of course, the proper dielectric constant
C              must be provided too.
C
        Tabs   = 168.0d0
        Eps    = 1.519d0
        EpsInf = 1.519d0
        DerEps = Zero
        RSolv  = 2.07d0
        VMol   = 42.71d0
        TCE    = 9.82d-2
        STen   = 42.86d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.410d-02
        NATyp  = 1
C     Atom Kr:
        NTT(1)    = 1
        RDiff(1) = 2.07d0
        KT(1)    = DKList(36)
        RWT(1)   = RList(36)
C
      Else If(IDSolv.eq.24) Then
C
C     Xenon
C
C     Warning: the following data are referred to the absolute
C              temperature of 210 K (the highest available). Change
C              the parameter Tabs to this value in actual calculations.
C              To use other temperatures, one has to provide also the
C              corresponding density (in g cm-3) to compute RHO (i.e.
C              the numeral density) and VMol (i.e. the molar volume)
C              properly. Of course, the proper dielectric constant
C              must be provided too.
C
        Tabs   = 210.0d0
        Eps    = 1.706d0
        EpsInf = 1.706d0
        DerEps = Zero
        RSolv  = 2.20d0
        VMol   = 50.38d0
        TCE    = 9.82d-2
        STen   = 42.86d0
        DSTen  = Zero
        CMF    = Zero
C     Atomic parameters for dispersion and repulsion
        Rho    = 1.195d-02
        NATyp  = 1
C     Atom Xe:
        NTT(1)    = 1
        RDiff(1) = 2.20d0
        KT(1)    = DKList(53)
        RWT(1)   = RList(53)
      Endif
C
C---- Use user specified value of the dielectric constant
C
      If (Eps_User.ne.-One) Eps=Eps_User
      If (EpsInf_User.ne.Zero) EpsInf=EpsInf_User
C
      Return
      End
C
      Function RList(IA)
      Implicit Real*8 (A-H,O-Z)
C
C     Assignes Caillet-Claverie's atomic parameters for dispersion and
C     repulsion.
C
      Parameter(MaxIA=110)
      Real*8 RW(MaxIA)
C
      Save RW
      Data RW/
C        H             He
     $ 1.20d+00,      1.28d+00,
C       Li            Be           B               C            N
     $ 0.00d+00,     0.00d+00,    0.00d+00,       1.70d+00,    1.60d+00,
C       O             F            Ne
     $ 1.50d+00,     1.45d+00,    1.38d+00,
C       Na            Mg           Al              Si           P
     $ 1.20d+00,     0.00d+00,    0.00d+00,       0.00d+00,    1.85d+00,
C       S             Cl           Ar
     $ 1.80d+00,     1.76d+00,    1.66d+00,
C       K             Ca
     $ 1.46d+00,     0.00d+00,
C       Sc            Ti           V               Cr           Mn
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Fe            Co           Ni              Cu           Zn
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ga            Ge           As              Se           Br
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    1.85d+00,
C       Kr
     $ 1.76d+00,
C       Rb            Sr
     $ 0.00d+00,     0.00d+00,
C       Y             Zr           Nb              Mo           Tc
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ru            Rh              Pd              Ag           Cd
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       In            Sn           Sb              Te           I
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    1.96d+00,
C       Xe
     $ 1.85d+00,
C       Cs            Ba
     $ 0.00d+00,     0.00d+00,
C       La            Ce           Pr              Nd           Pm
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Sm            Eu           Gd              Tb           Dy
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ho            Er           Tm               Yb           Lu
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Hf
     $ 0.00d+00,
C       Ta            W            Re              Os           Ir
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Pt            Au           Hg              Tl           Pb
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Bi            Po           At              Rn           Fr
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Ra            Ac           Th
     $ 0.00d+00,     0.00d+00,    0.00d+00,
C       Pa            U            Np              Pu           Am
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Cm            Bk           Cf              Es           Fm
     $ 0.00d+00,     0.00d+00,    0.00d+00,       0.00d+00,    0.00d+00,
C       Md
     $ 0.00d+00, 9*0.00d+00 /
C
      If(IA.gt.0.and.IA.le.MaxIA) then
        RList = RW(IA)
      Else
        RList = 0.0D0
        write(6,'(''IA out of range in RList.'')')
        Call Abend()
      EndIf
C
      Return
      End
C
      Function DKList(NumAt)
      Implicit Real*8 (A-H,O-Z)
C
C     Assignes Caillet-Claverie's atomic parameters for dispersion and
C     repulsion.
C
      Parameter(NElem=110)
      Real*8 K(NElem)
      Save K
      Data K/
C        H             He
     $ 1.00d+00,     0.594d+00,
C       Li            Be           B               C            N
     $ 0.000d+00,    0.000d+00,   0.000d+00,     1.000d+00,   1.180d+00,
C       O             F            Ne
     $ 1.360d+00,    1.500d+00,   1.209d+00,
C       Na            Mg           Al              Si           P
     $ 1.400d+00,    0.000d+00,   0.000d+00,     0.000d+00,   2.100d+00,
C       S             Cl           Ar
     $ 2.400d+00,    2.100d+00,   2.127d+00,
C       K             Ca
     $ 2.900d+00,    0.000d+00,
C       Sc            Ti           V               Cr           Mn
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Fe            Co           Ni              Cu           Zn
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Ga            Ge           As              Se           Br
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   2.400d+00,
C       Kr
     $ 2.566d+00,
C       Rb            Sr
     $ 0.000d+00,    0.000d+00,
C       Y             Zr           Nb              Mo           Tc
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Ru            Rh              Pd              Ag           Cd
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       In            Sn           Sb              Te           I
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   3.200d+00,
C       Xe
     $ 0.000d+00,
C       Cs            Ba
     $ 0.000d+00,    0.000d+00,
C       La            Ce           Pr              Nd           Pm
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Sm            Eu           Gd              Tb           Dy
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Ho            Er           Tm               Yb           Lu
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Hf
     $ 0.000d+00,
C       Ta            W            Re              Os           Ir
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Pt            Au           Hg              Tl           Pb
     $ 0.000d+00,    0.000d+00,   0.0005d+00,    0.000d+00,   0.000d+00,
C       Bi            Po           At              Rn           Fr
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Ra            Ac           Th
     $ 0.000d+00,    0.000d+00,   0.000d+00,
C       Pa            U            Np              Pu           Am
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Cm            Bk           Cf              Es           Fm
     $ 0.000d+00,    0.000d+00,   0.000d+00,     0.000d+00,   0.000d+00,
C       Md
     $ 0.000d+00, 9*0.00d+00 /
C
      DKList = K(NumAt)
C
      Return
      End
