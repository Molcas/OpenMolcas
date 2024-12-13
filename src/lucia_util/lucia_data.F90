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
Module Lucia_Data
Implicit None
#include "mxpdim.fh"
Private :: MXPIRR,MXPOBS,MXPR4T,MXPORB,MXPICI,MXPSTT,MXPCSM,MXPNGAS,MXPNSMST,MXPPTSPC,MXPTSOB

! Stuff fron cicisp.fh
Integer       IDUMMY,NICISP,                                      &
              NELCI(MXPICI),                                      &
              ISMOST(MXPCSM,MXPCSM),MXSB,MXSOOB,                  &
              NBLKIC(MXPCSM,MXPICI),LCOLIC(MXPCSM,MXPICI),        &
              MXNTTS,MXSOOB_AS
REAL*8        XISPSM(MXPCSM,MXPICI)

! Stuff fron cgas.fh
INTEGER     IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),                    &
            NGSOB(MXPOBS,MXPNGAS),                                &
            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,      &
            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,                     &
            NCMBSPC, LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),      &
            NMXOCCLS,IPHGAS(MXPNGAS),                             &
            IPHGAS1(MXPNGAS)

! Stuff fron cecore.fh
REAL*8        ECORE,ECORE_ORIG,ECORE_HEX

! Stuff from spinfo_lucia.fh
!./SPINFO/ new of Nov. 2001
! NCONF_TOT added april 29
Integer MINOP,MAXOP,NCONF_PER_SYM(MXPCSM),                   &
        NCONF_PER_OPEN(MXPORB+1,MXPCSM),                     &
        NPDTCNF(MXPORB+1),NPCSCNF(MXPORB+1),                 &
        NPCMCNF(MXPORB+1),                                   &
        IB_CONF_REO(MXPORB+1),IB_CONF_OCC(MXPORB+1),         &
        IB_SD_FOR_OPEN(MXPORB+1),                            &
        NCSF_PER_SYM(MXPCSM), NSD_PER_SYM(MXPCSM),           &
        NCONF_ALL_SYM, NCONF_ALL_SYM_FOR_OCCLS(MXPCSM),      &
        IBCONF_ALL_SYM_FOR_OCCLS(MXPCSM),nconf_tot,nCSF_HEXS

End Module Lucia_data
