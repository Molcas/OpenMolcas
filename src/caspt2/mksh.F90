!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      Subroutine MKSH()
      use caspt2_module, only: NSYM, NINDEP
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use constants, only: One
      use definitions, only: iwp, wp
      Implicit none
      real(kind=wp) Dum(1)
      integer(kind=iwp) ISYM,ICASE,IDISK,NIN
! For completeness, even case H has formally S and B
! matrices. This costs nothing, and saves conditional
! looping, etc in the rest  of the routines.
      DUM(1)=One
      DO ISYM=1,NSYM
        DO ICASE=12,13
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.GT.0) THEN
            IDISK=IDSMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,1,DUM,1,IDISK)
          END IF
        END DO
      END DO
      End Subroutine MKSH
