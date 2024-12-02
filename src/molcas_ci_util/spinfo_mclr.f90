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
Module spinfo_mclr
#include "detdim.fh"
Private MXPIRR,MXPOBS,MXPR4T,MXINKA,MXPORB,MXPXOT,MXPXST,MXPSHL, &
        MXPL,MXPXT,MXPICI,MXPSTT,MXPCSM,MXPCTP,MXCNSM,MXPWRD, &
        MXNMS,MTYP,MXPNGAS,MXPNSMST,MXPPTSPC
#include "spinfo_mclr.fh"
save
!             MULTSP                        : Spin multiplicity
!        MS2P                        : 2*MS
!             MINOP                        : Minum open orbitals
!        MAXOP                        : Maximum open orbitals
!        NTYP                        : Maxop-MinOp+1
!        NDPCNT(MXPCTP)                : Number of det. / conf type
!        NCPCNT(MXPCTP)                : Number of csf / conf type
!             NCNATS(MXPCTP,MXPCSM)        :
!        NDTASM(MXPCSM)                :       Combinations
!        NCSASM(MXPCSM)                :        CSF
!             NCNASM(MXPCSM)                :i Configurations
!
!     Integer       MULTSP,MS2P,                                        &
!    &              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),     &
!    &              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),&
!    &              NCNASM(MXPCSM)
End Module spinfo_mclr
