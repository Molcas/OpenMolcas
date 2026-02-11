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
! Root_Molcas, bas_Molcas, gssh_molcas, igsoccx_molcas, ELIMINATED_IN_GAS_MOLCAS, 2ELIMINATED_IN_GAS_MOLCAS, potnuc_Molcas,
! thre_Molcas, nsym_Molcas, nactel_Molcas, ms2_Molcas, ispin_Molcas, lsym_Molcas, itmax_Molcas, nroots_Molcas, ipt2_Molcas,
! iprci_Molcas, ngas_molcas, INOCALC_MOLCAS, ISAVE_EXP_MOLCAS, IEXPAND_MOLCAS, N_ELIMINATED_GAS_MOLCAS, N_2ELIMINATED_GAS_MOLCAS,
! I_ELIMINATE_GAS_MOLCAS, nCSF_HEXS
!
! stuff from bk_approx.fh
!
! DoBKAP, NGASBK, IOCCPSPC

use Molcas, only: MxGAS, MxSym
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MXTYP = 30, nSpeed = 5

integer(kind=iwp) :: I2ELIMINATED_IN_GAS_MOLCAS(MxGAS), I_ELIMINATE_GAS_MOLCAS, IELIMINATED_IN_GAS_MOLCAS(MxGAS), IEXPAND_MOLCAS, &
                     igsoccx_molcas(MxGAS,2), INOCALC_MOLCAS, IOCCPSPC(20,2), iprci_Molcas, ipt2_Molcas, ISAVE_EXP_MOLCAS, &
                     iSpeed(nSpeed), ispin_Molcas, itmax_Molcas, lsym_Molcas, MINOP, MS2, ms2_Molcas, N_2ELIMINATED_GAS_MOLCAS, &
                     N_ELIMINATED_GAS_MOLCAS, nactel_Molcas, NCNASM(mxSym), NCNFTP(MXTYP,mxSym), NCSASM(mxSym), nCSF_HEXS, &
                     NCSFTP(MXTYP), NDET, NDTASM(mxSym), NDTFTP(MXTYP), ngas_molcas, NGASBK, ngssh_molcas(MxGAS,mxSym), &
                     nroots_Molcas, nsym_Molcas, NTYP
real(kind=wp) :: potnuc_Molcas, thre_Molcas
logical(kind=iwp) :: DoBKAP, DoComb

public :: DoBKAP, DoComb, I2ELIMINATED_IN_GAS_MOLCAS, I_ELIMINATE_GAS_MOLCAS, IELIMINATED_IN_GAS_MOLCAS, IEXPAND_MOLCAS, &
          igsoccx_molcas, INOCALC_MOLCAS, IOCCPSPC, iprci_Molcas, ipt2_Molcas, ISAVE_EXP_MOLCAS, iSpeed, ispin_Molcas, &
          itmax_Molcas, lsym_Molcas, MINOP, MS2, ms2_Molcas, N_2ELIMINATED_GAS_MOLCAS, N_ELIMINATED_GAS_MOLCAS, nactel_Molcas, &
          NCNASM, NCNFTP, NCSASM, nCSF_HEXS, NCSFTP, NDET, NDTASM, NDTFTP, ngas_molcas, NGASBK, ngssh_molcas, nroots_Molcas, &
          nsym_Molcas, NTYP, potnuc_Molcas, thre_Molcas

end module Spinfo
