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
      SUBROUTINE MKBE(DREF,NDREF,FD)
      use definitions, only: iwp, wp
      use constants, only: Half, Two
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV, only: IDBMAT, IDSMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NASH,NAES,EPSA,EASUM
      IMPLICIT NONE

      INTEGER(KIND=IWP), INTENT(IN):: NDREF
      REAL(KIND=WP), INTENT(IN):: DREF(NDREF),FD(NDREF)

      Real(KIND=WP), ALLOCATABLE:: BE(:), S(:), SD(:)
      INTEGER(KIND=IWP) ISYM,NINP,NINM,NAS,NBE,NS,IDS,IDIAG,I,IT,ITABS, &
     &                  IX,IXABS,IBE,ID,IDISK,IDT
      Real(KIND=WP) VALUE,ET,EX

! Set up the matrix BE(t,x)
! Formula used:
!    BE(t,x)=-Ftx + (EASUM-Ex-Et)*Dtx
!            + 2dtx Ex


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,6)
        IF(NINP.EQ.0) CYCLE
        NINM=NINDEP(ISYM,7)
        NAS=NASH(ISYM)
        NBE=(NAS*(NAS+1))/2
        IF(NBE.GT.0) THEN
          CALL mma_allocate(BE,NBE,LABEL='BE')
!GG.Nov03  Load in SD the diagonal elements of SE matrix:
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,Label='S')
          CALL mma_allocate(SD,NAS,Label='SD')
          IDS=IDSMAT(ISYM,6)
          CALL DDAFILE(LUSBT,2,S,NS,IDS)
          IDIAG=0
          DO I=1,NAS
            IDIAG=IDIAG+I
            SD(I)=S(IDIAG)
          END DO
          CALL mma_deallocate(S)
        ENDIF
!GG End
        DO IT=1,NAS
          ITABS=IT+NAES(ISYM)
          ET=EPSA(ITABS)
          DO IX=1,IT
            IXABS=IX+NAES(ISYM)
            EX=EPSA(IXABS)
            IBE=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
            VALUE=-FD(ID)+(EASUM-EX-ET)*DREF(ID)
            IF(ITABS.EQ.IXABS) THEN
              VALUE=VALUE+Two*EX
            END IF
!GG.Nov03
            IF (IT.eq.IX) THEN
              IDT=(ITABS*(ITABS+1))/2
              VALUE=VALUE+ipea_shift*Half*DREF(IDT)*SD(IT)
            ENDIF
!GG End
            BE(IBE)=VALUE
          END DO
        END DO

! Write to disk
        IF(NBE.GT.0.and.NINDEP(ISYM,6).GT.0) THEN
          IDISK=IDBMAT(ISYM,6)
          CALL DDAFILE(LUSBT,1,BE,NBE,IDISK)
          IF(NINM.GT.0.and.NINDEP(ISYM,7).GT.0) THEN
            IDISK=IDBMAT(ISYM,7)
            CALL DDAFILE(LUSBT,1,BE,NBE,IDISK)
          END IF
          CALL mma_deallocate(BE)
!GG.Nov03 DisAlloc SD
          CALL mma_deallocate(SD)
!GG End
        END IF
      END DO

      END SUBROUTINE MKBE
