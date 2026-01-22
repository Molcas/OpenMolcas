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

module Jobiph_J

use Molcas, only: MxRoot, MxSym
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iadr15_j(15), iroot_j(mxroot), ispin_j, lroots_j, lsym_j, nactel_j, nash_j(mxsym), nbas_j(mxsym), nconf_j, &
                     ndel_j(mxsym), nelec3_j, nfro_j(mxsym), nhole1_j, nish_j(mxsym), nroots_j, nrs1_j(mxsym), nrs2_j(mxsym), &
                     nrs3_j(mxsym), nsym_j
real(kind=wp) :: weight_j(mxroot)
character(len=72) :: title_j(18)

public :: iadr15_j, iroot_j, ispin_j, lroots_j, lsym_j, nactel_j, nash_j, nbas_j, nconf_j, ndel_j, nelec3_j, nfro_j, nhole1_j, &
          nish_j, nroots_j, nrs1_j, nrs2_j, nrs3_j, nsym_j, title_j, weight_j

end module Jobiph_J
