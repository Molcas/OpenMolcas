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

implicit none
private

#include "real.fh"
public :: Zero, One, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten, Eleven, Twelve
public :: Half, Quart, OneHalf, Pi, SqrtP2, TwoP34, TwoP54, One2C2

#include "constants2.fh"
public :: diel, rPi, UTOAU, elmass, ATOKG, elcharge, rNAVO, cLight, auTocm, rPlanck, kBoltzmann, rBohr, cm_s, Debye, Angstrom, RF, &
          auToHz, auTofs, auToN, auToeV, auTokJ, c_in_au

end module Constants
