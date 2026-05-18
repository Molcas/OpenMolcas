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
      SUBROUTINE MKBG(DREF,NDREF,FD)
      use definitions, only: iwp, wp
      use constants, only: half, two
      use caspt2_global, only:ipea_shift
      use caspt2_global, only:LUSBT
      use EQSOLV, only: IDSMAT,IDBMAT
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NINDEP,NASH,NAES,EASUM
      IMPLICIT NONE

      INTEGER(KIND=IWP), INTENT(IN):: NDREF
      REAL(KIND=WP), INTENT(IN):: DREF(NDREF),FD(NDREF)

      REAL(KIND=WP), ALLOCATABLE:: BG(:), S(:), SD(:)
      INTEGER(KIND=IWP) ISYM,NINP,NINM,NAS,NBG,NS,IDS,IDIAG,I,IT,ITABS, &
     &                  IX,IXABS,IBG,ID,IDISK,IDT
      REAL(KIND=WP) VALUE

!     Set up the matrix BG(t,x)
!     Formula used:
!     BG(t,x)= Ftx -EASUM*Dtx


      DO ISYM=1,NSYM
        NINP=NINDEP(ISYM,10)
        IF(NINP.EQ.0) CYCLE
        NINM=NINDEP(ISYM,11)
        NAS=NASH(ISYM)
        NBG=(NAS*(NAS+1))/2
        IF(NBG.GT.0) THEN
          CALL mma_Allocate(BG,NBG,LABEL='BG')
!GG.Nov03  Load in SD the diagonal elements of SG matrix:
          NS=(NAS*(NAS+1))/2
          CALL mma_allocate(S,NS,Label='S')
          CALL mma_allocate(SD,NAS,Label='SD')
          IDS=IDSMAT(ISYM,10)
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
          DO IX=1,IT
            IXABS=IX+NAES(ISYM)
            IBG=(IT*(IT-1))/2+IX
            ID=(ITABS*(ITABS-1))/2+IXABS
!GG.Nov03
            VALUE=FD(ID)-EASUM*DREF(ID)
            IF (IT.eq.IX) THEN
              IDT=(ITABS*(ITABS+1))/2
              VALUE = VALUE +                                           &
     &                ipea_shift*half*(two-DREF(IDT))*SD(IT)
            ENDIF
            BG(IBG)=VALUE
!           BG(BG)=FD(ID)-EASUM*DREF(ID)
!GG End
          END DO
        END DO

! Write to disk, and save size and address.
        IF(NBG.GT.0) THEN
         IF(NINDEP(ISYM,10).GT.0) THEN
          IDISK=IDBMAT(ISYM,10)
          CALL DDAFILE(LUSBT,1,BG,NBG,IDISK)
         END IF
         IF(NINM.GT.0.and.NINDEP(ISYM,11).GT.0) THEN
           IDISK=IDBMAT(ISYM,11)
           CALL DDAFILE(LUSBT,1,BG,NBG,IDISK)
         END IF
!GG.Nov03 DisAlloc SD
         CALL mma_deallocate(SD)
!GG End
         CALL mma_deallocate(BG)
        END IF
      END DO

      END SUBROUTINE MKBG
