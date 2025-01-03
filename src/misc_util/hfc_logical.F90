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
Module HFC_logical
use Definitions, only: IWP
Private
!      UHF_HFC: to control the calculation of hyperfine coupling tensor
!      matrix in unrestricted Hartree-Fock scf calculations.
!      It is mainly used in scf program and the related integral_util
!      directory.
!
!      MAG_X2C: true when x2c-transformed hyperfine magnetic integrals
!      are calculated.

Logical(kind=iwp),Public:: UHF_HFC, MAG_X2C
End Module HFC_logical
