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
      SUBROUTINE MKSE(DREF,NDREF)
      use definitions, only: iwp, wp
      use constants, only: Two
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NASH,NAES
      IMPLICIT NONE

      INTEGER(kind=iwp), intent(in)::  NDREF
      REAL(kind=wp), intent(in)::  DREF(NDREF)

      REAL(kind=wp), ALLOCATABLE:: SE(:)
      INTEGER(kind=iwp) ISYM,NINP,NINM,NAS,NSE,IT,ITABS,IX,IXABS,ISE,   &
     &                  ID,IDISK
! Set up the matrix SE(t,x)
! Formula used:
!    SE(t,x)=2*dtx - Dtx


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,6)
        IF(NINP.EQ.0) CYCLE
        NINM=NINDEP(ISYM,7)
        NAS=NASH(ISYM)
        NSE=(NAS*(NAS+1))/2
        IF(NSE.GT.0) CALL mma_allocate(SE,NSE,Label='SE')
        DO IT=1,NAS
          ITABS=IT+NAES(ISYM)
          DO IX=1,IT
            IXABS=IX+NAES(ISYM)
            ISE=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
            IF(ITABS.EQ.IXABS) THEN
              SE(ISE)=Two-DREF(ID)
            ELSE
              SE(ISE)=-DREF(ID)
            END IF
          END DO
        END DO

! Write to disk
        IF(NSE.GT.0.and.NINDEP(ISYM,6).GT.0) THEN
          IDISK=IDSMAT(ISYM,6)
          CALL DDAFILE(LUSBT,1,SE,NSE,IDISK)
          IF(NINM.GT.0.and.NINDEP(ISYM,7).GT.0) THEN
            IDISK=IDSMAT(ISYM,7)
            CALL DDAFILE(LUSBT,1,SE,NSE,IDISK)
          END IF
          CALL mma_deallocate(SE)
        END IF
      END DO

      END SUBROUTINE MKSE
