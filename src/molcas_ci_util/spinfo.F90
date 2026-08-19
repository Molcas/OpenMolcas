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
! MXTYP, MS2, MINOP, NTYP, NDTFTP, NCSFTP, NCNFTP
!
! stuff from ciinfo.fh
!
! ICOMBI, NDET, NDTASM, NCSASM, NCNASM
!
! stuff from lucia_ini.fh
!
! i = 1, 5: combinations, particle hole(sigma), count_aa, count_ab, a/p_parts
! nSpeed, iSpeed
!
! ELIMINATED_IN_GAS, 2ELIMINATED_IN_GAS,
! nactel, ms2, ispin, STSYM, itmax, nroots, ipt2,
! iprci, ngas, INOCALC, ISAVE_EXP, N_ELIMINATED_GAS, N_2ELIMINATED_GAS,
! I_ELIMINATE_GAS, nCSF_HEXS

use Molcas, only: MxGAS, MxSym
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MXTYP = 30, nSpeed = 5

integer(kind=iwp) :: I2ELIMINATED_IN_GAS(MxGAS), I_ELIMINATE_GAS, IELIMINATED_IN_GAS(MxGAS), &
                     INOCALC, iprci, ipt2, ISAVE_EXP, iSpeed(nSpeed), &
                     ispin, itmax, STSYM, MINOP, MS2, N_2ELIMINATED_GAS, &
                     N_ELIMINATED_GAS, nactel, NCNASM(mxSym), NCNFTP(MXTYP,mxSym), NCSASM(mxSym), nCSF_HEXS, &
                     NCSFTP(MXTYP), NDET, NDTASM(mxSym), NDTFTP(MXTYP), ngas, nroots, NTYP
logical(kind=iwp) :: DoComb

public :: DoComb, I2ELIMINATED_IN_GAS, I_ELIMINATE_GAS, IELIMINATED_IN_GAS, &
          INOCALC, iprci, ipt2, ISAVE_EXP, iSpeed, ispin, itmax, STSYM, MINOP, &
          MS2, N_2ELIMINATED_GAS, N_ELIMINATED_GAS, nactel, NCNASM, NCNFTP, NCSASM, nCSF_HEXS, &
          NCSFTP, NDET, NDTASM, NDTFTP, ngas, nroots, NTYP

end module Spinfo
