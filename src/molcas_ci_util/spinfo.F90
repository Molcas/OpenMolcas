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

End Module Spinfo
