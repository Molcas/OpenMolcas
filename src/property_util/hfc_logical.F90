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

module HFC_logical

use Definitions, only: iwp

implicit none
private

logical(kind=iwp) :: MagX2C_Avail, UHF_HFC

! VARIABLE DESCRIPTION
!
!  UHF_HFC     : to control the calculation of hyperfine coupling tensor
!               matrix in unrestricted Hartree-Fock scf calculations.
!               It is mainly used in scf program and the related integral_util
!               directory.
!
! MagX2C_Avail : is controlled by iRdOne.
!             = .TRUE. when MagX2C integral is accessible and vice versa

public :: MagX2C_Avail, UHF_HFC

end module HFC_logical
