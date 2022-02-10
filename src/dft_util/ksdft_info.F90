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

module KSDFT_Info

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ifav, ifav_n, ifiv, ifiv_n, LuMC, LuMT
real(kind=wp) :: CoefR = One, CoefX = One, FA_time, FI_time, Funcaa = Zero, Funcbb = Zero, Funccc = Zero, PUVX_time, sp_time
logical(kind=iwp) :: do_pdftpot
character(len=80) :: KSDFA

public :: ifav, ifav_n, ifiv, ifiv_n, CoefR, CoefX, FA_time, FI_time, Funcaa, Funcbb, Funccc, KSDFA, LuMC, LuMT, PUVX_time, &
          sp_time, do_pdftpot

end module KSDFT_Info
