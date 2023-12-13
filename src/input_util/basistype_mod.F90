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
! Copyright (C) 2001-2005, Valera Veryazov                             *
!***********************************************************************

module BasisType_Mod

implicit none
private

!
!  'strange' order is due to back compatibility issues
!    UNK: unknown
!
character(len=*), parameter :: BasTypeCon = 'SEG:ANO:RAF:CCC:UNK:UNC:ARC:GEN:SEC:'
! 1: Contraction type
!  - SEG: segmented
!  - SEC: segmented Cartesian
!  - ANO: ANO
!  - RAF: Raffenetti
!  - CCC: CC contraction
!  - UNC: uncontracted
!  - ARC: automatic Raffenetti
!  - GEN: general contraction

character(len=*), parameter :: BasTypeAll = 'AE_:NAE:YES:NO_:UNK:'
! 2: All electron
!  - AE_: all electron
!  - NAE: not all electron
!  - YES: same as AE_
!  - NO_: same as NAE

character(len=*), parameter :: BasTypeRel = 'NRH:RH_:RCP:DKH:UNK:DK2:DK3:DK4:DK5:DK6:DK7:DK8:RYD:X2C:'
! 3: Relativity
!  - RH_: relativistic Hamiltonian
!  - NRH: non relativistic Hamiltonian
!  - RCP: relativistic core potential
!  - DKH: DKH
!  - X2C: exact 2-component
!  - RYD: Rydberg

character(len=*), parameter :: BasTypeNuc = 'PN_:GN_:MGN:UNK:'
! 4: Nucleus
!  - PN_: point nucleus
!  - GN_: Gaussian nucleus
!  - MGN: modified Gaussian

public :: BasTypeCon, BasTypeAll, BasTypeRel, BasTypeNuc

end module BasisType_Mod
