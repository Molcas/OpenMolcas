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
! MXTYP, MINOP, NTYP, NDTFTP, NCSFTP, NCNFTP
!
! stuff from ciinfo.fh
!
! ICOMBI (the CI dimensions NDTASM, NCSASM are in CI_Info, CI_St(iState)%nDet,
! %nCSF; NCNASM was not needed)
!
! stuff from lucia_ini.fh
!
! i = 1, 5: combinations, particle hole(sigma), count_aa, count_ab, a/p_parts
! nSpeed, iSpeed
!
! (the HEXS/DEXS requests are in CI_Info, Space_Inp, the number of kept CSFs
! in CI_St(iState)%nCSF_Kept)

use Molcas, only: MxSym
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MXTYP = 30, nSpeed = 5

integer(kind=iwp) :: iSpeed(nSpeed), MINOP, NCNFTP(MXTYP,mxSym), NCSFTP(MXTYP), NDTFTP(MXTYP), NTYP
logical(kind=iwp) :: DoComb

public :: DoComb, iSpeed, MINOP, NCNFTP, NCSFTP, NDTFTP, NTYP

end module Spinfo
