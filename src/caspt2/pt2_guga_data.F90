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
!     This module is a temporary tweak while the whole include file
!     pt2_guga.fh is converted to a proper module which is address
!     using the "only" construct in the "use" statement of the module.
      module pt2_guga_data
      use gugx, only: IA0, IB0, IC0, NLEV, NVERT0, NVERT,               &
     &                LV1RAS,LM1RAS,LV3RAS,LM3RAS,                      &
     &                NOW1, NNOW, IOW1, NIOW,                           &
     &                MIDLEV,NMIDV,MIDV1,MIDV2,                         &
     &                NOCSF, NNOCSF, IOCSF, NIOCSF
      implicit none
#include "pt2_guga.fh"
      Integer, Allocatable:: MVL(:), MVR(:)
      Integer, Allocatable:: NOCP(:), IOCP(:)
      Integer, Allocatable:: ICASE(:), ICOUP(:)
      Real*8,  Allocatable:: VTAB(:)
      save
      end module pt2_guga_data
