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
Module STRBAS
      Private
#include "mxpdim.fh"
      Integer       KOCSTR(MXPSTT),KNSTSO(MXPSTT),KISTSO(MXPSTT),       &
     &              KSTSTM(MXPSTT,2),KZ(MXPSTT),                        &
     &              KSTREO(MXPSTT),                                     &
     &              KCOBSM,KNIFSJ,KIFSJ,KIFSJO,KSTSTX,                  &
     &              KNSTSGP(MXPNGAS),KISTSGP(MXPNGAS) ,                 &
     &              KIOCLS,KSPGPAN,KSPGPCR

      Public        KOCSTR,KNSTSO,KISTSO,                               &
     &              KSTSTM,KZ,                                          &
     &              KSTREO,                                             &
     &              KCOBSM,KNIFSJ,KIFSJ,KIFSJO,KSTSTX,                  &
     &              KNSTSGP,KISTSGP,                                    &
     &              KIOCLS,KSPGPAN,KSPGPCR
End Module STRBAS
