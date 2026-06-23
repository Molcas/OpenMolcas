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

integer(kind=iwp) :: MagX2C_Req
logical(kind=iwp) :: MagX2C_Avail, UHF_HFC

! VARIABLE DESCRIPTION
!
!  UHF_HFC     : to control the calculation of hyperfine coupling tensor
!               matrix in unrestricted Hartree-Fock scf calculations.
!               It is mainly used in scf program and the related integral_util
!               directory.
!
! MagX2C_Req   : an integer specifies whether RX2C, MTXC are requested in &SEWARD
!              = -2 (non-relativistic case), both are turned off: no RX2C, no MXTC
!              = -1 RX2C = ON, but users do not specify MXTC
!              =  2 both are turned on
!
! MagX2C_Avail : is controlled by iRdOne.
!             = .TRUE. when MagX2C integral is accessible and vice versa
!
! NOTE: MagX2C_Req is used for handling input instead of calling iRdOne.

public :: MagX2C_Avail, MagX2C_Req, UHF_HFC

end module HFC_logical
