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
Module Spinfo
Implicit None
!stuff from spinfo.fh
INTEGER, PARAMETER :: MXTYP=30, MXSM=8
INTEGER MULTS,MS2,MINOP,MAXOP,NTYP,NDTFTP(MXTYP),NCSFTP(MXTYP),     &
        NCNFTP(MXTYP,MXSM),NCONF_TOT

!stuff from ciinfo.fh
INTEGER, PARAMETER:: MXCISM=8
INTEGER IORB1F,IORB1L,NEL1MN,NEL1MX,IORB3F,IORB3L,NEL3MN,NEL3MX,    &
        MXSASM,MXVBLK,ICOMBI,NDET,NDTASM(MXCISM),NCSASM(MXCISM),NCNASM(MXCISM)

#include "Molcas.fh"
Private :: MaxBfn,MaxBfn_Aux,MxAO,mxAtom,mxroot,mxNemoAtom,Mxdbsc,lCache,mxact,mxina,mxbas,mxOrb,mxSym,mxGAS, &
           LENIN,LENIN1,LENIN2,LENIN3,LENIN4,LENIN5,LENIN6,LENIN8
#include "lucia_ini.fh"
!stuff from lucia_ini.fh
!!
!! Transfers parameters from inpctl in rasscf to lucia_ini in lucia_util
!!
!Integer RtoI_Molcas
!
!! i = 1, 5: combinations, particle hole(sigma), count_aa, count_ab, a/p_parts
!Integer, Parameter:: nSpeed = 5
!Integer iSpeed(nSpeed)
!
!Integer iRoot_Molcas(mxRoot),norb_Molcas(mxSym)
!Integer nbas_Molcas(mxSym),nish_Molcas(mxSym)
!Integer ngssh_molcas(mxgas,mxsym),igsoccx_molcas(mxgas,2)
!Integer IELIMINATED_IN_GAS_MOLCAS(MXGAS)
!Integer I2ELIMINATED_IN_GAS_MOLCAS(MXGAS)
!
!Real*8   potnuc_Molcas,thre_Molcas
!Integer  nsym_Molcas,nactel_Molcas, ms2_Molcas,                   &
!         ispin_Molcas,lsym_Molcas,nhole1_Molcas,nelec3_Molcas,    &
!         itmax_Molcas,nroots_Molcas,ipt2_Molcas,iprci_Molcas,     &
!         ngas_molcas, INOCALC_MOLCAS,ISAVE_EXP_MOLCAS,            &
!         IEXPAND_MOLCAS,N_ELIMINATED_GAS_MOLCAS,                  &
!         N_2ELIMINATED_GAS_MOLCAS,                                &
!         I_ELIMINATE_GAS_MOLCAS,nCSF_HEXS

End Module Spinfo
