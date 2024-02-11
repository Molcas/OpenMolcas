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
Module GUGX
Private
INTEGER, Public:: NLEV,IA0,IB0,IC0,                               &
                  NVERT0,NDRT0,LDRT0,NDOWN0,LDOWN0,                 &
                  IFCAS,LV1RAS,LM1RAS,LV3RAS,LM3RAS,                &
                  NVERT,NDRT,LDRT,NDOWN,LDOWN,                      &
                  LUP,NUP,LRAW,NRAW,LDAW,NDAW,                      &
                  MIDLEV,NMIDV,MIDV1,MIDV2,NUW,NLW,MXUP,MXDWN,      &
                  NWALK,NNOW,NIOW,                        &
                  NIPWLK,NICASE,LICASE,NCSF(8),                     &
                  NNOCSF,LNOCSF,NIOCSF,LIOCSF,                      &
                  LLSGN,LUSGN
Integer, Allocatable, Public:: NOW1(:), IOW1(:)
End Module GUGX
