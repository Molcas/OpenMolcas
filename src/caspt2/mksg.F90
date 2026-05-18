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
      SUBROUTINE MKSG(DREF,NDREF)
      use definitions, only: iwp, wp
      use caspt2_global, only: LUSBT
      use EQSOLV, only: IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NASH,NAES
      IMPLICIT None

      INTEGER(kind=iwp), intent(in)::  NDREF
      REAL(kind=wp), intent(in)::  DREF(NDREF)

      REAL(kind=wp), ALLOCATABLE:: SG(:)
      INTEGER(kind=iwp) ISYM,NINP,NINM,NAS,NSG,IT,ITABS,IX,IXABS,ISG,   &
     &                  ID,IDISK
! Set up the matrix SG(t,x)
! Formula used:
!    SG(t,x)= Dtx


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,10)
        IF(NINP.EQ.0) CYCLE
        NINM=NINDEP(ISYM,11)
        NAS=NASH(ISYM)
        NSG=(NAS*(NAS+1))/2
        IF(NSG.GT.0) CALL mma_allocate(SG,NSG,Label='SG')
        DO IT=1,NAS
          ITABS=IT+NAES(ISYM)
          DO IX=1,IT
            IXABS=IX+NAES(ISYM)
            ISG=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
            SG(ISG)= DREF(ID)
          END DO
        END DO

! Write to disk
        IF(NSG.GT.0.and.NINDEP(ISYM,10).GT.0) THEN
          IDISK=IDSMAT(ISYM,10)
          CALL DDAFILE(LUSBT,1,SG,NSG,IDISK)
          IF(NINM.GT.0.and.NINDEP(ISYM,11).GT.0) THEN
            IDISK=IDSMAT(ISYM,11)
            CALL DDAFILE(LUSBT,1,SG,NSG,IDISK)
          END IF
          CALL mma_deallocate(SG)
        END IF
      END DO

      END SUBROUTINE MKSG
