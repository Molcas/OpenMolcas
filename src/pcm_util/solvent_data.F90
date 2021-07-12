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

module Solvent_Data
! Database of optical and physical data for various solvents.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: Initialized = .false.
integer(kind=iwp), parameter :: MxT = 4 ! Maximum number of atom types per solvent

type AtSolv
  integer(kind=iwp) :: NTT = 0
  real(kind=wp) :: RDiff = Zero, KT = Zero, RWT = Zero
end type AtSolv

type Solvent
  character(len=32) :: SName
  real(kind=wp) :: Eps, EpsInf, DerEps, RSolv, VMol, TCE !, STen, DSTen, CMF, Rho, Tabs
  type(AtSolv) :: Atoms(MxT)
end type Solvent

!Caillet-Claverie's atomic parameters for dispersion and repulsion.
!    H       He
!   Li       Be        B        C        N       O         F       Ne
!   Na       Mg       Al       Si        P       S        Cl       Ar
!    K       Ca
!            Sc       Ti        V       Cr       Mn       Fe       Co       Ni       Cu       Zn
!                     Ga       Ge       As       Se       Br       Kr
!   Rb       Sr
!             Y       Zr       Nb       Mo       Tc       Ru       Rh       Pd       Ag       Cd
!                     In       Sn       Sb       Te        I       Xe
! All zero (no data) after Xe
real(kind=wp), parameter :: DK(54) = [ &
  1.000_wp,0.594_wp, &
    Zero  ,  Zero  ,  Zero  ,1.000_wp,1.180_wp,1.360_wp,1.500_wp,1.209_wp, &
  1.400_wp,  Zero  ,  Zero  ,  Zero  ,2.100_wp,2.400_wp,2.100_wp,2.127_wp, &
  2.900_wp,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,2.400_wp,2.566_wp, &
    Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,3.200_wp,  Zero &
]
real(kind=wp), parameter :: R(54) = [ &
  1.20_wp ,1.28_wp , &
    Zero  ,  Zero  ,  Zero  ,1.70_wp ,1.60_wp ,1.50_wp ,1.45_wp ,1.38_wp , &
  1.20_wp ,  Zero  ,  Zero  ,  Zero  ,1.85_wp ,1.80_wp ,1.76_wp ,1.66_wp , &
  1.46_wp ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,1.85_wp ,1.76_wp , &
    Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,1.96_wp ,1.85_wp &
]

! Other atomic radii
!    H        He
!   Li        Be         B         C         N         O         F        Ne
!   Na        Mg        Al        Si         P         S        Cl        Ar
!    K        Ca
!             Sc        Ti         V        Cr        Mn        Fe        Co        Ni        Cu        Zn
!                       Ga        Ge        As        Se        Br        Kr
!   Rb        Sr (+2)
!             Y  (+3)   Zr (+4)   Nb (+5)   Mo (+6)   Tc (+5)   Ru (+2)   Rh (+3)   Pd (+2)   Ag (+1)   Cd (+2)
!                       In        Sn        Sb        Te         I        Xe
!   Cs        Ba (+2)
!             La (+3)   Ce (+3)   Pr (+3)   Nd (+3)   Pm (+3)   Sm (+3)   Eu (+3)
!             Gd (+3)   Tb (+3)   Dy (+3)   Ho (+3)   Er (+3)   Tm (+3)   Yb (+3)
!             Lu (+3)   Hf (+4)   Ta (+5)    W(+4,+6) Re(+5,+7) Os (+6)   Ir (+3)   Pt        Au        Hg
!                       Tl        Pb        Bi (+3)   Po (+2)   At        Rn (+4)
!   Fr        Ra (+2)
!             Ac (+3)   Th (+4)   Pa (+4)    U (+4)   Np (+4)   Pu (+4)   Am (+4)
!             Cm (+3)   Bk (+3)   Cf (+3)   Es (+3)   Fm (+3)   Md (+3)   No (+3)
!             Lw (+3)   Ru

! Pauling radii
real(kind=wp), parameter :: PRad(104) = [ &
  1.20_wp ,1.20_wp , &
  1.37_wp ,1.45_wp ,1.45_wp ,1.50_wp ,1.50_wp ,1.40_wp ,1.35_wp ,1.30_wp , &
  1.57_wp ,1.36_wp ,1.24_wp ,1.17_wp ,1.90_wp ,1.85_wp ,1.80_wp ,1.88_wp , &
  2.75_wp ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,1.63_wp ,1.40_wp ,1.39_wp , &
                    1.87_wp ,1.86_wp ,2.00_wp ,2.00_wp ,1.95_wp ,2.02_wp , &
    Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,1.63_wp ,1.72_wp ,1.58_wp , &
                    1.93_wp ,2.17_wp ,2.20_wp ,2.20_wp ,2.15_wp ,2.16_wp , &
    Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,1.72_wp ,1.66_wp ,1.55_wp , &
                    1.96_wp ,1.02_wp ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
    Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,1.86_wp ,  Zero  ,  Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
             Zero  ,  Zero &
]
! Warning: the following are atomic diameters.
real(kind=wp), parameter :: UFFRad(0:104) = [ &
    Zero  , &
  2.886_wp,2.362_wp, &
  2.451_wp,2.745_wp,4.083_wp,3.851_wp,3.660_wp,3.500_wp,3.364_wp,3.243_wp, &
  2.983_wp,3.021_wp,4.499_wp,4.295_wp,4.147_wp,4.035_wp,3.947_wp,3.868_wp, &
  3.812_wp,3.399_wp, &
           3.295_wp,3.175_wp,3.144_wp,3.023_wp,2.961_wp,2.912_wp,2.872_wp,2.834_wp,3.495_wp,2.763_wp, &
                    4.383_wp,4.280_wp,4.230_wp,4.205_wp,4.189_wp,4.141_wp, &
  4.114_wp,3.641_wp, &
           3.345_wp,3.124_wp,3.165_wp,3.052_wp,2.998_wp,2.963_wp,2.929_wp,2.899_wp,3.148_wp,2.848_wp, &
                    4.463_wp,4.392_wp,4.420_wp,4.470_wp,4.50_wp,4.404_wp, &
  4.517_wp,3.703_wp, &
           3.522_wp,3.556_wp,3.606_wp,3.575_wp,3.547_wp,3.520_wp,3.493_wp, &
           3.368_wp,3.451_wp,3.428_wp,3.409_wp,3.391_wp,3.374_wp,3.355_wp, &
           3.640_wp,3.141_wp,3.170_wp,3.069_wp,2.954_wp,3.120_wp,2.840_wp,2.754_wp,3.293_wp,2.705_wp, &
                    4.347_wp,4.297_wp,4.370_wp,4.709_wp,4.750_wp,4.765_wp, &
  4.900_wp,3.677_wp, &
           3.478_wp,3.396_wp,3.424_wp,3.395_wp,3.424_wp,3.424_wp,3.381_wp, &
           3.326_wp,3.339_wp,3.313_wp,3.299_wp,3.286_wp,3.274_wp,3.248_wp, &
           3.236_wp,3.500_wp &
]
! Covalent radii
! Parameters for atoms heavier than At are
! taken from the UFF force field (A.K.Rappe',C.J.Casewit,K.S.Colwell,
! W.A.Goddard III and W.M.Skiff, J.Am.Chem.Soc. 114,10024 (1992)).
real(kind=wp), parameter :: CRad(0:104) = [ &
    Zero  , &
  0.354_wp,0.849_wp, &
  1.336_wp,1.010_wp,0.838_wp,0.757_wp,0.700_wp,0.658_wp,0.668_wp,0.920_wp, &
  1.539_wp,1.421_wp,1.244_wp,1.117_wp,1.101_wp,1.064_wp,1.044_wp,1.032_wp, &
  1.953_wp,1.761_wp, &
           1.513_wp,1.412_wp,1.402_wp,1.345_wp,1.382_wp,1.270_wp,1.241_wp,1.164_wp,1.302_wp,1.193_wp, &
                    1.260_wp,1.197_wp,1.211_wp,1.190_wp,1.192_wp,1.147_wp, &
  2.260_wp,2.052_wp, &
           1.698_wp,1.564_wp,1.473_wp,1.467_wp,1.322_wp,1.478_wp,1.332_wp,1.338_wp,1.386_wp,1.403_wp, &
                    1.459_wp,1.398_wp,1.407_wp,1.386_wp,1.382_wp,1.267_wp, &
  2.570_wp,2.277_wp, &
           1.943_wp,1.841_wp,1.823_wp,1.816_wp,1.801_wp,1.780_wp,1.771_wp, &
           1.735_wp,1.732_wp,1.710_wp,1.696_wp,1.673_wp,1.660_wp,1.637_wp, &
           1.671_wp,1.611_wp,1.511_wp,1.392_wp,1.372_wp,1.372_wp,1.371_wp,1.364_wp,1.262_wp,1.340_wp, &
                    1.518_wp,1.459_wp,1.512_wp,1.500_wp,1.545_wp,1.420_wp, &
  2.880_wp,2.512_wp, &
           1.983_wp,1.721_wp,1.711_wp,1.684_wp,1.666_wp,1.657_wp,1.660_wp, &
           1.801_wp,1.761_wp,1.750_wp,1.724_wp,1.712_wp,1.689_wp,1.679_wp, &
           1.698_wp,1.850_wp &
]

type(Solvent), protected :: SolvData(24)

public :: Init_Solvent_Data, Pauling, RCov97, SolvData

contains

subroutine Init_Solvent_Data()
  ! Skip if already initialized
  if (Initialized) return

  ! For each solvent, the following data is given:
  !     SName
  !     Eps
  !     EpsInf
  !     DerEps
  !     RSolv
  !     VMol
  !     TCE
  !     (commented out: STen, DSTen, CMF, Rho, and Tabs for some)
  ! and then, for each atom type in the solvent:
  !     NTT
  !     RDiff
  !     KT
  !     RWT
  ! (using default initialization for derived types does not play nice with all compilers)

  ! Water
  !                     Eps,     EpsInf,  DerEps,    RSolv,   VMol,    TCE,             STen,    DSTen,   CMF,     Rho
  SolvData(1) = Solvent('WATER', &
                        78.39_wp,1.776_wp,-0.3562_wp,1.385_wp,18.07_wp,2.57e-4_wp, & ! ,71.81_wp,0.650_wp,1.277_wp,3.348e-2_wp)
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(1,1.5_wp,DK(8),R(8)), & ! Atom O
                         AtSolv(2,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(0,Zero,Zero,Zero), &
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! Acetonitrile (Warning: RSolv and VMol for this solvent need to be checked)
  !                     Eps,     EpsInf,  DerEps,RSolv,   VMol,    TCE,              STen,    DSTen,   CMF,     Rho
  SolvData(2) = Solvent('ACETONITRILE', &
                        36.64_wp,1.806_wp, Zero ,2.155_wp,53.68_wp,1.192e-3_wp, & ! ,36.47_wp,1.373_wp,0.808_wp,1.153e-2_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(2,1.76_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(3,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(1,1.6_wp,DK(7),R(7)), & ! Atom N
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! Methanol
  !                     Eps,     EpsInf,  DerEps,    RSolv,   VMol,   TCE,              STen,    DSTen,   CMF,     Rho
  SolvData(3) = Solvent('METHANOL', &
                        32.63_wp,1.758_wp,-0.1984_wp,1.855_wp,40.7_wp,1.182e-3_wp, & ! ,22.12_wp,1.154_wp,1.776_wp,1.495e-2_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(1,1.76_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(4,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(1,1.5_wp,DK(8),R(8)), & ! Atom O
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! Ethanol
  !                     Eps,     EpsInf,  DerEps,    RSolv,   VMol,   TCE,              STen,    DSTen,   CMF,     Rho
  SolvData(4) = Solvent('ETHANOL', &
                        24.55_wp,1.847_wp,-0.1510_wp,2.180_wp,58.7_wp,1.103e-3_wp, & ! ,21.89_wp,1.146_wp,1.543_wp,1.032e-2_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(2,1.76_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(6,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(1,1.5_wp,DK(8),R(8)), & ! Atom O
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! IsoQuinoline
  !                     Eps,     EpsInf,  DerEps,RSolv,  VMol,     TCE,              STen,    DSTen, CMF,   Rho
  SolvData(5) = Solvent('ISOQUINOLINE', &
                        10.43_wp,1.010_wp, Zero ,3.50_wp,117.27_wp,1.255e-3_wp, & ! ,26.53_wp, Zero , Zero ,5.135e-3_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(9,1.5_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(7,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(1,1.6_wp,DK(7),R(7)), & ! Atom N
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! Quinoline
  !                     Eps,    EpsInf,  DerEps,RSolv,  VMol,     TCE,              STen,    DSTen, CMF,   Rho
  SolvData(6) = Solvent('QUINOLINE', &
                        9.03_wp,1.010_wp, Zero ,3.50_wp,117.27_wp,1.255e-3_wp, & ! ,26.53_wp, Zero , Zero ,5.135e-3_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(9,1.5_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(7,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(1,1.6_wp,DK(7),R(7)), & ! Atom N
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! Chloroform
  !                     Eps,    EpsInf,  DerEps,RSolv,  VMol,   TCE,              STen,    DSTen, CMF,   Rho
  SolvData(7) = Solvent('CHLOROFORM', &
                        4.90_wp,2.085_wp, Zero ,2.48_wp,80.7_wp,1.255e-3_wp, & ! ,26.53_wp, Zero , Zero ,7.482e-3_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(1,2.82_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(1,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(3,1.79_wp,DK(17),R(17)), & ! Atom Cl
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! EthylEther
  !                     Eps,     EpsInf,DerEps,RSolv,   VMol,      TCE,             STen,  DSTen, CMF,   Rho
  SolvData(8) = Solvent('ETHYLETHER', &
                        4.335_wp, Zero , Zero ,2.785_wp,103.84_wp,1.617e-3_wp, & ! , Zero , Zero , Zero ,5.799e-3_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(4,1.76_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(10,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(1,1.50_wp,DK(8),R(8)), & ! Atom O
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! MethyleneChloride
  !                     Eps,    EpsInf,  DerEps,RSolv,  VMol,   TCE,              STen,    DSTen, CMF,   Rho
  SolvData(9) = Solvent('METHYLENECHLORIDE', &
                        8.93_wp,2.020_wp, Zero ,2.27_wp,64.5_wp,1.367e-3_wp, & ! ,27.33_wp, Zero , Zero ,9.74e-3_wp
                        ! Atomic parameters for dispersion and repulsion
                        [AtSolv(1,1.76_wp,DK(6),R(6)), & ! Atom C
                         AtSolv(2,1.2_wp,DK(1),R(1)), & ! Atom H
                         AtSolv(2,1.79_wp,DK(17),R(17)), & ! Atom Cl
                         AtSolv(0,Zero,Zero,Zero) &
                        ])

  ! DiChloroEthane
  !                      Eps,     EpsInf,  DerEps,RSolv,   VMol,   TCE,              STen,    DSTen, CMF,   Rho
  SolvData(10) = Solvent('DICHLOROETHANE', &
                         10.36_wp,2.085_wp, Zero ,2.505_wp,79.4_wp,1.156e-3_wp, & ! ,31.54_wp, Zero , Zero ,7.576e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(2,1.76_wp,DK(6),R(6)), &  ! Atom C
                          AtSolv(4,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(2,1.79_wp,DK(17),R(17)), & ! Atom Cl
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! CarbonTetraChloride
  !                      Eps,     EpsInf,  DerEps,RSolv,   VMol,   TCE,              STen,    DSTen,   CMF,     Rho
  SolvData(11) = Solvent('CARBONTETRACHLORIDE', &
                         2.228_wp,2.129_wp, Zero ,2.685_wp,96.5_wp,1.270e-3_wp, & ! ,26.15_wp,1.436_wp,0.629_wp,6.241e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(1,2.82_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(4,1.79_wp,DK(17),R(17)), & ! Atom Cl
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Benzene
  !                      Eps,     EpsInf,  DerEps,RSolv,  VMol,    TCE,              STen,    DSTen,   CMF,     Rho
  SolvData(12) = Solvent('BENZENE', &
                         2.247_wp,2.244_wp, Zero ,2.63_wp,88.91_wp,1.380e-3_wp, & ! ,28.18_wp,1.469_wp,0.629_wp,6.773e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(6,1.5_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(6,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Toluene
  !                      Eps,     EpsInf,  DerEps,RSolv,  VMol,    TCE,             STen,    DSTen,   CMF,     Rho
  SolvData(13) = Solvent('TOLUENE', &
                         2.379_wp,2.232_wp, Zero ,2.82_wp,106.3_wp,1.08e-3_wp, & ! ,27.92_wp,1.391_wp,0.679_wp,5.665e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(7,1.5_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(8,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Chlorobenzene
  !                      Eps,     EpsInf,  DerEps,RSolv,   VMol,     TCE,              STen,    DSTen, CMF,   Rho
  SolvData(14) = Solvent('CHLOROBENZENE', &
                         5.621_wp,2.320_wp, Zero ,2.805_wp,101.79_wp,0.981e-3_wp, & ! ,32.69_wp, Zero , Zero ,5.916e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(6,1.5_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(5,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(1,1.79_wp,DK(17),R(17)), & ! Atom Cl
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! NitroMethane
  !                      Eps,     EpsInf,  DerEps,RSolv,   VMol,    TCE,              STen,    DSTen,   CMF,     Rho
  SolvData(15) = Solvent('NITROMETHANE', &
                         38.20_wp,1.904_wp, Zero ,2.155_wp,53.68_wp,1.192e-3_wp, & ! ,36.47_wp,1.373_wp,0.808_wp,1.122e-2_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(1,1.76_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(3,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(1,1.6_wp,DK(7),R(7)), & ! Atom N
                          AtSolv(2,1.5_wp,DK(8),R(8)) & ! Atom O
                         ])

  ! Heptane
  !                      Eps,    EpsInf,  DerEps,RSolv,   VMol,     TCE,             STen,    DSTen,   CMF,     Rho
  SolvData(16) = Solvent('HEPTANE', &
                         1.92_wp,1.918_wp, Zero ,3.125_wp,146.56_wp,1.25e-3_wp, & ! ,19.80_wp,1.505_wp,0.687_wp,4.109e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(7,1.76_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(16,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Cyclohexane
  !                      Eps,     EpsInf,  DerEps,RSolv,   VMol,     TCE,             STen,    DSTen,   CMF,     Rho
  SolvData(17) = Solvent('CYCLOHEXANE', &
                         2.023_wp,2.028_wp, Zero ,2.815_wp,108.10_wp,1.20e-3_wp, & ! ,24.38_wp,1.467_wp,0.621_wp,5.571e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(6,1.76_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(12,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Aniline
  !                      Eps,    EpsInf,  DerEps,RSolv,  VMol,    TCE,             STen,    DSTen,   CMF,     Rho
  SolvData(18) = Solvent('ANILINE', &
                         6.89_wp,2.506_wp, Zero ,2.80_wp,91.15_wp,0.85e-3_wp, & ! ,42.79_wp,0.731_wp,0.972_wp,6.607e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(6,1.5_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(7,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(1,1.6_wp,DK(7),R(7)), & ! Atom N
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Acetone
  !                      Eps,    EpsInf,  DerEps,     RSolv,  VMol,    TCE,             STen,    DSTen, CMF,   Rho
  SolvData(19) = Solvent('ACETONE', &
                         20.7_wp,1.841_wp,-9.77e-2_wp,2.38_wp,73.52_wp,1.42e-3_wp, & ! ,22.67_wp, Zero , Zero ,8.190e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(3,1.76_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(6,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(1,1.5_wp,DK(8),R(8)), & ! Atom O
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! TetraHydroFuran
  !                      Eps,    EpsInf,  DerEps,RSolv,  VMol,    TCE,              STen,    DSTen, CMF,   Rho
  SolvData(20) = Solvent('TETRAHYDROFURAN', &
                         7.58_wp,1.971_wp, Zero ,2.56_wp,81.11_wp,1.142e-3_wp, & ! ,26.40_wp, Zero , Zero ,7.425e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(4,1.76_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(8,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(1,1.5_wp,DK(8),R(8)), & ! Atom O
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! DiMethylSulfoxide
  !                      Eps,    EpsInf,  DerEps,    RSolv,   VMol,    TCE,             STen,    DSTen, CMF,   Rho
  SolvData(21) = Solvent('DIMETHYLSULFOXIDE', &
                         46.7_wp,2.179_wp,-0.1902_wp,2.455_wp,70.94_wp,9.82e-2_wp, & ! ,42.86_wp, Zero , Zero ,8.49e-3_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(2,1.76_wp,DK(6),R(6)), & ! Atom C
                          AtSolv(6,1.2_wp,DK(1),R(1)), & ! Atom H
                          AtSolv(1,1.80_wp,DK(16),R(16)), & ! Atom S
                          AtSolv(1,1.5_wp,DK(8),R(8)) & ! Atom O
                         ])

  ! Argon
  ! Warning: the following data are referred to the absolute
  !          temperature of 118 K (the highest available). Change
  !          the parameter Tabs to this value in actual calculations.
  !          To use other temperatures, one has to provide also the
  !          corresponding density (in g cm-3) to compute RHO (i.e.
  !          the numeral density) and VMol (i.e. the molar volume)
  !          properly. Of course, the proper dielectric constant
  !          must be provided too.
  !                      Eps,     EpsInf,  DerEps,RSolv,   VMol,    TCE,             STen,    DSTen, CMF,   Rho,        Tabs
  SolvData(22) = Solvent('ARGON', &
                         1.430_wp,1.430_wp, Zero ,1.875_wp,34.29_wp,9.82e-2_wp, & ! ,42.86_wp, Zero , Zero ,1.756e-2_wp,118.0_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(1,1.875_wp,DK(18),R(18)), & ! Atom Ar
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Krypton
  ! Warning: the following data are referred to the absolute
  !          temperature of 168 K (the highest available). Change
  !          the parameter Tabs to this value in actual calculations.
  !          To use other temperatures, one has to provide also the
  !          corresponding density (in g cm-3) to compute RHO (i.e.
  !          the numeral density) and VMol (i.e. the molar volume)
  !          properly. Of course, the proper dielectric constant
  !          must be provided too.
  !                      Eps,     EpsInf,  DerEps,RSolv,  VMol,    TCE,             STen,    DSTen, CMF,   Rho,        Tabs
  SolvData(23) = Solvent('KRYPTON', &
                         1.519_wp,1.519_wp, Zero ,2.07_wp,42.71_wp,9.82e-2_wp, & ! ,42.86_wp, Zero , Zero ,1.410e-2_wp,168.0_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(1,2.07_wp,DK(36),R(36)), & ! Atom Kr
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  ! Xenon
  ! Warning: the following data are referred to the absolute
  !          temperature of 210 K (the highest available). Change
  !          the parameter Tabs to this value in actual calculations.
  !          To use other temperatures, one has to provide also the
  !          corresponding density (in g cm-3) to compute RHO (i.e.
  !          the numeral density) and VMol (i.e. the molar volume)
  !          properly. Of course, the proper dielectric constant
  !          must be provided too.
  !                      Eps,     EpsInf,  DerEps,RSolv,  VMol,    TCE,             STen,    DSTen, CMF,   Rho,        Tabs
  SolvData(24) = Solvent('XENON', &
                         1.706_wp,1.706_wp, Zero ,2.20_wp,50.38_wp,9.82e-2_wp, & ! ,42.86_wp, Zero , Zero ,1.195e-2_wp,210.0_wp
                         ! Atomic parameters for dispersion and repulsion
                         [AtSolv(1,2.20_wp,DK(53),R(53)), & ! Atom Xe
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero), &
                          AtSolv(0,Zero,Zero,Zero) &
                         ])

  Initialized = .true.

end subroutine Init_Solvent_Data

function Pauling(Z)
! Pauling radius for atomic number Z (UFF when not defined)
  real(kind=wp) :: Pauling
  integer(kind=iwp), intent(in) :: Z
  Pauling = PRad(Z)
  if (Pauling == Zero) Pauling = UFFRad(Z)
end function Pauling

function RCov97(Za,Zb)
! This function returns an estimated covalent bond distance (in Ang)
! between atoms of atomic numbers Za and Zb. Setting Zb to 0 returns
! the covalent radius of Za.
  real(kind=wp) :: RCov97
  integer(kind=iwp), intent(in) :: Za, Zb
  integer(kind=iwp) :: i, j
  i = min(max(Za,lbound(CRad,1)),ubound(CRad,1))
  j = min(max(Zb,lbound(CRad,1)),ubound(CRad,1))
  RCov97 = CRad(i)+CRad(j)
end function RCov97

end module Solvent_Data
