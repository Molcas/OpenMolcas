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

module Spinfo

! stuff from spinfo.fh
!
! MXTYP, MINOP, NTYP, NDTFTP, NCSFTP, NCNFTP (only MCLR's own Lucia copy now; RASSCF's
! tables are in lucia_data)
!
! stuff from ciinfo.fh
!
! ICOMBI (the CI dimensions NDTASM, NCSASM are in CI_Info, CI_St(iState)%nDet,
! %nCSF; NCNASM was not needed)
!
! (the HEXS/DEXS requests are in CI_Info, Space_Inp, the number of kept CSFs
! in CI_St(iState)%nCSF_Kept)

use Molcas, only: MxSym
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MXTYP = 30

integer(kind=iwp) :: MINOP, NCNFTP(MXTYP,mxSym), NCSFTP(MXTYP), NDTFTP(MXTYP), NTYP

public :: MINOP, NCNFTP, NCSFTP, NDTFTP, NTYP

end module Spinfo
