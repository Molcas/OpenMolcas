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
      Type iArray
         Integer, Allocatable:: I(:)
      End Type iArray
      Type (iArray):: OCSTR(MXPSTT), STREO(MXPSTT), STSTM(MXPSTT,2),    &
     &                NSTSGP(MXPNGAS), ISTSGP(MXPNGAS), NSTSO(MXPSTT),  &
     &                ISTSO(MXPSTT), Zmat(MXPSTT)
      Integer, Allocatable:: IOCLS(:), SPGPAN(:), SPGPCR(:)

      Public ::     OCSTR,NSTSO,ISTSO,STSTM,Zmat,STREO,                 &
     &              NSTSGP,ISTSGP,IOCLS,SPGPAN,SPGPCR
End Module STRBAS
