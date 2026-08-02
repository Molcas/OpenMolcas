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

module Constants

use Definitions, only: wp

implicit none
private

! Common numbers

real(kind=wp), parameter :: Zero = 0.0_wp, One = 1.0_wp, Two = 2.0_wp, Three  = 3.0_wp, Four = 4.0_wp, Five = 5.0_wp, &
                            Six = 6.0_wp, Seven = 7.0_wp, Eight = 8.0_wp, Nine = 9.0_wp, Ten = 10.0_wp, Eleven = 11.0_wp, &
                            Twelve = 12.0_wp, Half = 0.5_wp, Quart = 0.25_wp, OneHalf = 1.5_wp

complex(kind=wp), parameter :: cZero = (Zero,Zero), cOne = (One,Zero), Onei = (Zero,One)

! Pi and related constants

real(kind=wp), parameter :: Pi = Four*atan2(One,One), TwoP34 = One/(Two*Pi)**(Three/Four), TwoP54 = sqrt(Two)*Pi**(Five/Four)

!***********************************************************************
!                                                                      *
!              Physical constants and conversion factors               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! The constants and conversion factors in this file are taken from     *
! NIST at URL https://physics.nist.gov/cuu/Constants                   *
!                                                                      *
! CODATA 2010: Peter J. Mohr, Barry N. Taylor, David B. Newell,        *
!              CODATA recommended values of the fundamental physical   *
!              constants: 2010                                         *
!              Rev. Mod. Phys. vol. 84(4) 1527-1605 (2012).            *
!              doi:10.1103/RevModPhys.84.1527                          *
!                                                                      *
! CODATA 2014: Peter J. Mohr, David B. Newell, Barry N. Taylor,        *
!              CODATA recommended values of the fundamental physical   *
!              constants: 2014                                         *
!              Rev. Mod. Phys. vol. 88(3) 035009 (2016).               *
!              doi:10.1103/RevModPhys.88.035009                        *
!                                                                      *
! CODATA 2018: Eite Tiesinga, Peter J. Mohr, David B. Newell,          *
!              Barry N. Taylor,                                        *
!              CODATA recommended values of the fundamental physical   *
!              constants: 2018                                         *
!              Rev. Mod. Phys. vol. 93(2) 025010 (2021).               *
!              doi:10.1103/RevModPhys.93.025010                        *
!                                                                      *
! CODATA 2022: Eite Tiesinga, Peter J. Mohr, David B. Newell,          *
!              Barry N. Taylor,                                        *
!              The 2022 CODATA Recommended Values of the Fundamental   *
!              Physical Constants                                      *
!              (Web Version 9.0)                                       *
!                                                                      *
!----------------------------------------------------------------------*

! Constants
! ---------
! CONST_AMU_IN_SI_               -- Atomic mass unit (1/12*m[C12]) in SI units.
! CONST_ATM_IN_SI_               -- 1atm in SI units.
! CONST_AU_TIME_IN_SI_           -- Atomic units time unit in SI units.
! CONST_AU_VELOCITY_IN_SI_       -- Atomic units velocity unit in SI units.
! CONST_AVOGADRO_                -- Avogadro's number in SI units.
! CONST_BOHR_MAGNETON_IN_SI_     -- Bohr magneton in SI units.
! CONST_BOHR_RADIUS_IN_SI_       -- Bohr radius in SI units.
! CONST_BOLTZMANN_               -- Boltzmann's constant in SI units.
! CONST_C_IN_AU_                 -- Speed of light in atomic units.
! CONST_C_IN_SI_                 -- Speed of light in SI units.
! CONST_DIELECTRIC_IN_SI_        -- Electric constant (eps0) in SI units.
! CONST_DIPOLE_IN_SI_            -- Atomic unit of dipole (e*a0) in SI units.
! CONST_ELECTRON_G_FACTOR_       -- Electron g-factor (dimensionless).
! CONST_ELECTRON_MASS_IN_SI_     -- Mass of the electron in SI units.
! CONST_ELEMENTARY_CHARGE_IN_SI_ -- Elementary charge in SI units.
! CONST_MOLAR_GAS_               -- Molar gas constant in SI units.
! CONST_MUON_MASS_IN_SI_         -- Mass of the muon in SI units.
! CONST_PLANCK_                  -- Planck's constant in SI units.
!
! Conversion factors
! ------------------
! CONV_AMU_TO_AU_                -- Convert 1AMU to au.
! CONV_ATM_TO_AU_                -- Convert 1tm to au.
! CONV_AU_TO_CM1_                -- Convert 1Eh to reciprocal cm.
! CONV_AU_TO_DEBYE               -- Convert 1au to debye.
! CONV_AU_TO_EV_                 -- Convert 1Eh to electron-volt.
! CONV_AU_TO_HZ_                 -- Convert 1Eh to hertz.
! CONV_AU_TO_KJ_                 -- Convert 1Eh to kilojoule.
! CONV_AU_TO_KJ_PER_MOLE_        -- Convert 1Eh to kilojoule per mole.
! CONV_AU_TO_T                   -- Convert 1au to tesla.
! CONV_CAL_TO_J_                 -- Convert 1cal to joule.

!***********************************************************************
! Select a set of physical constants, valid values are:
! 2010, 2014, 2018, 2022

#ifndef CODATA_SET
#define CODATA_SET 2022
#endif

! Constants

#if CODATA_SET == 2022
#  define CONST_AMU_IN_SI_                1.66053906892e-27_wp
#  define CONST_AU_TIME_IN_SI_            2.4188843265864e-17_wp
#  define CONST_AU_VELOCITY_IN_SI_        2.18769126216e6_wp
#  define CONST_AVOGADRO_                 6.02214076e23_wp
#  define CONST_BOHR_RADIUS_IN_SI_        5.29177210544e-11_wp
#  define CONST_ELECTRON_MASS_IN_SI_      9.1093837139e-31_wp
#  define CONST_MUON_MASS_IN_SI_          1.883531627e-28_wp
#  define CONST_BOLTZMANN_                1.380649e-23_wp
#  define CONST_ELEMENTARY_CHARGE_IN_SI_  1.602176634e-19_wp
#  define CONST_PLANCK_                   6.62607015e-34_wp
#  define CONST_DIELECTRIC_IN_SI_         8.8541878188e-12_wp
#  define CONST_ELECTRON_G_FACTOR_       -2.00231930436092_wp
#  define CONST_BOHR_MAGNETON_IN_SI_      9.2740100657e-24_wp
#elif CODATA_SET == 2018
! CONST_AVOGADRO_, CONST_BOLTZMANN_, CONST_ELEMENTARY_CHARGE_IN_SI_, CONST_PLANCK_ are now exact
#  define CONST_AMU_IN_SI_                1.66053906660e-27_wp
#  define CONST_AU_TIME_IN_SI_            2.4188843265857e-17_wp
#  define CONST_AU_VELOCITY_IN_SI_        2.18769126364e6_wp
#  define CONST_AVOGADRO_                 6.02214076e23_wp
#  define CONST_BOHR_RADIUS_IN_SI_        5.29177210903e-11_wp
#  define CONST_ELECTRON_MASS_IN_SI_      9.1093837015e-31_wp
#  define CONST_MUON_MASS_IN_SI_          1.883531627e-28_wp
#  define CONST_BOLTZMANN_                1.380649e-23_wp
#  define CONST_ELEMENTARY_CHARGE_IN_SI_  1.602176634e-19_wp
#  define CONST_PLANCK_                   6.62607015e-34_wp
#  define CONST_DIELECTRIC_IN_SI_         8.8541878128e-12_wp
#  define CONST_ELECTRON_G_FACTOR_       -2.00231930436256_wp
#  define CONST_BOHR_MAGNETON_IN_SI_      9.2740100783e-24_wp
#elif CODATA_SET == 2014
#  define CONST_AMU_IN_SI_                1.660539040e-27_wp
#  define CONST_AU_TIME_IN_SI_            2.418884326509e-17_wp
#  define CONST_AU_VELOCITY_IN_SI_        2.18769126277e6_wp
#  define CONST_AVOGADRO_                 6.022140857e23_wp
#  define CONST_BOHR_RADIUS_IN_SI_        5.2917721067e-11_wp
#  define CONST_ELECTRON_MASS_IN_SI_      9.10938356e-31_wp
#  define CONST_MUON_MASS_IN_SI_          1.883531594e-28_wp
#  define CONST_MOLAR_GAS_                8.3144598_wp
#  define CONST_ELEMENTARY_CHARGE_IN_SI_  1.6021766208e-19_wp
#  define CONST_PLANCK_                   6.626070040e-34_wp
#  define CONST_DIELECTRIC_IN_SI_         8.854187817e-12_wp
#  define CONST_ELECTRON_G_FACTOR_       -2.00231930436182_wp
#  define CONST_BOHR_MAGNETON_IN_SI_      9.274009994e-24_wp
#elif CODATA_SET == 2010
#  define CONST_AMU_IN_SI_                1.660538921e-27_wp
#  define CONST_AU_TIME_IN_SI_            2.418884326502e-17_wp
#  define CONST_AU_VELOCITY_IN_SI_        2.18769126379e6_wp
#  define CONST_AVOGADRO_                 6.02214129e23_wp
#  define CONST_BOHR_RADIUS_IN_SI_        5.2917721092e-11_wp
#  define CONST_ELECTRON_MASS_IN_SI_      9.10938291e-31_wp
#  define CONST_MUON_MASS_IN_SI_          1.883531475e-28_wp
#  define CONST_MOLAR_GAS_                8.3144621_wp
#  define CONST_ELEMENTARY_CHARGE_IN_SI_  1.602176565e-19_wp
#  define CONST_PLANCK_                   6.62606957e-34_wp
#  define CONST_DIELECTRIC_IN_SI_         8.85418782e-12_wp
#  define CONST_ELECTRON_G_FACTOR_       -2.00231930436153_wp
#  define CONST_BOHR_MAGNETON_IN_SI_      9.27400968e-24_wp
#endif
#define CONST_C_IN_SI_                    2.99792458e8_wp
#define CONST_ATM_IN_SI_                  101325.0_wp

! Derived constants

#if CODATA_SET >= 2018
#  define CONST_MOLAR_GAS_                (CONST_AVOGADRO_*CONST_BOLTZMANN_)
#else
#  define CONST_BOLTZMANN_                (CONST_MOLAR_GAS_/CONST_AVOGADRO_)
#endif
#define CONST_C_IN_AU_                    (CONST_C_IN_SI_/CONST_AU_VELOCITY_IN_SI_)
#define CONST_DIPOLE_IN_SI_               (CONST_ELEMENTARY_CHARGE_IN_SI_*CONST_BOHR_RADIUS_IN_SI_)

! Conversion factors

#if CODATA_SET == 2022
#  define CONV_AU_TO_EV_                  27.211386245981_wp
#  define CONV_AU_TO_CM1_                 2.1947463136314e5_wp
#  define CONV_AU_TO_KJ_                  4.3597447222060e-21_wp
#  define CONV_AU_TO_HZ_                  6.5796839204999e15_wp
#  define CONV_AU_TO_T_                   2.35051757077e5_wp
#elif CODATA_SET == 2018
#  define CONV_AU_TO_EV_                  27.211386245988_wp
#  define CONV_AU_TO_CM1_                 2.1947463136320e5_wp
#  define CONV_AU_TO_KJ_                  4.3597447222071e-21_wp
#  define CONV_AU_TO_HZ_                  6.579683920502e15_wp
#  define CONV_AU_TO_T_                   2.35051756758e5_wp
#elif CODATA_SET == 2014
#  define CONV_AU_TO_EV_                  27.21138602_wp
#  define CONV_AU_TO_CM1_                 2.194746313702e5_wp
#  define CONV_AU_TO_KJ_                  4.359744650e-21_wp
#  define CONV_AU_TO_HZ_                  6.579683920711e15_wp
#  define CONV_AU_TO_T_                   2.350517550e5_wp
#elif CODATA_SET == 2010
#  define CONV_AU_TO_EV_                  27.21138505_wp
#  define CONV_AU_TO_CM1_                 2.194746313708e5_wp
#  define CONV_AU_TO_KJ_                  4.35974434e-21_wp
#  define CONV_AU_TO_HZ_                  6.579683920729e15_wp
#  define CONV_AU_TO_T_                   2.350517464e5_wp
#endif
#define CONV_CAL_TO_J_                    4.184_wp

! Derived conversion factors

#define CONV_AU_TO_DEBYE_                 (CONST_DIPOLE_IN_SI_*CONST_C_IN_SI_*1.0e21_wp)
#define CONV_AMU_TO_AU_                   (CONST_AMU_IN_SI_/CONST_ELECTRON_MASS_IN_SI_)
#define CONV_AU_TO_KJ_PER_MOLE_           (CONST_AVOGADRO_*CONV_AU_TO_KJ_)
#define CONV_ATM_TO_AU_                   (CONST_ATM_IN_SI_*CONST_BOHR_RADIUS_IN_SI_**3/CONV_AU_TO_KJ_*1.0e-3_wp)

! Fortran constants

real(kind=wp), parameter :: Angstrom = CONST_BOHR_RADIUS_IN_SI_*1.0e10_wp, &
                            atmToau = CONV_ATM_TO_AU_, &
                            atmToPa = CONST_ATM_IN_SI_, &
                            ATokg = CONST_AMU_IN_SI_, &
                            auTocm = CONV_AU_TO_CM1_, &
                            auToeV = CONV_AU_TO_EV_, &
                            auTofs = CONST_AU_TIME_IN_SI_*1.0e15_wp, &
                            auToHz = CONV_AU_TO_HZ_, &
                            auTokcalmol = CONV_AU_TO_KJ_PER_MOLE_/CONV_CAL_TO_J_, &
                            auTokJ = CONV_AU_TO_KJ_, &
                            auTokJmol = CONV_AU_TO_KJ_PER_MOLE_, &
                            auTokJmolnm = 1.0e-9*CONV_AU_TO_KJ_PER_MOLE_/CONST_BOHR_RADIUS_IN_SI_, &
                            auToN = CONV_AU_TO_KJ_/CONST_BOHR_RADIUS_IN_SI_*1.0e3_wp, &
                            auToT = CONV_AU_TO_T_, &
                            c_in_au = CONST_C_IN_AU_, &
                            cal_to_J = CONV_CAL_TO_J_, &
                            cLight = CONST_C_IN_SI_, &
                            cm_s = CONST_C_IN_SI_*1.0e2_wp, &
                            Debye = CONV_AU_TO_DEBYE_, &
                            deg2rad = Pi/180.0_wp, &
                            diel = CONST_DIELECTRIC_IN_SI_, &
                            elcharge = CONST_ELEMENTARY_CHARGE_IN_SI_, &
                            elmass = CONST_ELECTRON_MASS_IN_SI_, &
                            gElectron = CONST_ELECTRON_G_FACTOR_, &
                            hPlanck = CONST_PLANCK_, &
                            kBoltzmann = CONST_BOLTZMANN_, &
                            mBohr = CONST_BOHR_MAGNETON_IN_SI_, &
                            mu2elmass = CONST_MUON_MASS_IN_SI_/CONST_ELECTRON_MASS_IN_SI_, &
                            rBohr = CONST_BOHR_RADIUS_IN_SI_, &
                            Rgas = CONST_MOLAR_GAS_, &
                            rNAVO = CONST_AVOGADRO_, &
                            uToau = CONV_AMU_TO_AU_

public :: Angstrom, atmToau, atmToPa, ATokg, auTocm, auToeV, auTofs, auToHz, auTokcalmol, auTokJ, auTokJmol, auTokJmolnm, auToN, &
          auToT, c_in_au, cal_to_J, cLight, cm_s, cOne, cZero, Debye, deg2rad, diel, Eight, elcharge, Eleven, elmass, Five, Four, &
          gElectron, Half, hPlanck, kBoltzmann, mBohr, mu2elmass, Nine, One, OneHalf, Onei, Pi, Quart, rBohr, Rgas, rNAVO, Seven, &
          Six, Ten, Three, Twelve, Two, TwoP34, TwoP54, uToau, Zero

end module Constants
