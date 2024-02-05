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
Module GLBBAS
Private
#include "mxpdim.fh"
      Integer ::                                    KLSM2,KRHO1,        &
     &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR,             &
     &              KSBIA,KSBIB,KPNIJ,KIJKK,    KINH1,                  &
     &              KMOAOIN,KMOAOUT,NOCCLS_G,                           &
     &              KMOMO,KSRHO1,                                       &
     &              KINT1_SIMTRH,KINT2_SIMTRH,                          &
     &              KPINT1_SIMTRH,KPINT2_SIMTRH,KINH1_NOCCSYM
      Real*8, Allocatable:: INT1(:), INT1O(:)
      Integer, Allocatable:: PINT1(:), PINT2(:)
      Real*8, Allocatable:: VEC3(:)

! DETERMINE BASE ADRESSES
!             DFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
!             CFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
!             DTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
!             CONF_OCC(I) : SPACE FOR STORING  NCNSM
!                            CONFIGURATION EXPANSIONS
      Integer, Allocatable:: DFTP(:)
      Integer, Allocatable:: CFTP(:)
      Real*8,  Allocatable:: DTOC(:)
      Type iArray
         Integer, Allocatable:: I(:)
      End Type  iArray
      Type (iArray) ::  CONF_OCC(8), CONF_REO(8), SDREO_I(8),           &
     &                  PGINT1A(MXPOBS), PGINT1(MXPOBS)
      Type (iArray), Allocatable:: Z_PTDT(:), REO_PTDT(:)
      Integer, Allocatable:: LSM1(:)

      Integer, Allocatable:: KLOCCLS(:)

      Public         INT1,       PINT1, PINT2, LSM1,KLSM2,KRHO1,        &
     &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR,             &
     &              KSBIA,KSBIB, VEC3,KPNIJ,KIJKK,    KINH1,            &
     &              KMOAOIN,KMOAOUT,KLOCCLS,NOCCLS_G, PGINT1,           &
     &               INT1O,      PGINT1A,                               &
     &              KMOMO,KSRHO1,                                       &
     &              KINT1_SIMTRH,KINT2_SIMTRH,                          &
     &              KPINT1_SIMTRH,KPINT2_SIMTRH,KINH1_NOCCSYM,          &
     &              CONF_OCC,CONF_REO,                                  &
     &              DFTP,CFTP,DTOC,                                     &
     &              SDREO_I,Z_PTDT,REO_PTDT
End Module GLBBAS
