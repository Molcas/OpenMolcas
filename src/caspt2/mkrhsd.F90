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
!--------------------------------------------
! 1998  PER-AAKE MALMQUIST
! DEPARTMENT OF THEORETICAL CHEMISTRY
! UNIVERSITY OF LUND
! SWEDEN
!--------------------------------------------
      SUBROUTINE MKRHSD(IVEC,FIMO,NFIMO,ERI1,nERI1,ERI2,nERI2,SCR,nSCR)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Zero
      USE SUPERINDEX, only: KTU
      use fake_GA, only: GA_Arrays, Allocate_GA_Array,                  &
     &                            Deallocate_GA_Array
      use caspt2_module, only: NSYM,NINDEP,NSSH,NISH,NTU,NISUP,         &
     &                         NACTEL,NORB,NAES,NASH,NTUES
      IMPLICIT None
      integer(kind=iwp), intent(in):: IVEC, NFIMO,nERI1,nERI2,nSCR
      real(kind=wp), intent(inout):: FIMO(NFIMO)
      real(kind=wp), intent(inout):: ERI1(nERI1),ERI2(nERI2), SCR(nSCR)

      integer(kind=iwp) IOFF(8)
      integer(kind=iwp) ISYM,IO,ISYMI,NAS1,NAS,NIS,NV,LW,NFSUM,         &
     &                  NFIMOES,ISYMA,ISYMU,ISYMT,II,IU,IUABS,          &
     &                  IUTOT,IA,IATOT,IT,ITABS,ITTOT,IWA,IWI,IW1,IW2,  &
     &                  IBUF1,IBUF2,ICASE
      real(kind=wp) ONEADD,FAI,WAITU
! Set up RHS vector of PT2 Linear Equation System, in vector
! number IVEC of LUSOLV, for case 5, AIVX.

      DO ISYM=1,NSYM
        IF(NINDEP(ISYM,5).EQ.0) CYCLE
! Set up offset table:
          IO=0
          DO ISYMA=1,NSYM
            IOFF(ISYMA)=IO
            ISYMI=Mul(ISYMA,ISYM)
            IO=IO+NSSH(ISYMA)*NISH(ISYMI)
          END DO
!   Allocate W; W subdivided into W1,W2.
          NAS1=NTU(ISYM)
          NAS=2*NAS1
          NIS=NISUP(ISYM,5)
          NV=NAS*NIS
          IF(NV.EQ.0) CYCLE
! Compute W1(tu,ai)=(ai,tu) + FIMO(a,i)*delta(t,u)/NACTEL
! Compute W2(tu,ai)=(ti,au)
          LW=Allocate_GA_Array(NV,'WD')
          NFSUM=0
          DO ISYMI=1,NSYM
            NFIMOES=NFSUM
            NFSUM=NFSUM+(NORB(ISYMI)*(NORB(ISYMI)+1))/2
            ISYMA=Mul(ISYMI,ISYM)
            DO ISYMU=1,NSYM
              ISYMT=Mul(ISYMU,ISYM)
              DO II=1,NISH(ISYMI)
                DO IU=1,NASH(ISYMU)
                  IUABS=IU+NAES(ISYMU)
                  IUTOT=IU+NISH(ISYMU)
                  CALL EXCH(ISYMA,ISYMI,ISYMT,ISYMU,                    &
     &                      II,IUTOT,ERI1,SCR)
                  CALL EXCH(ISYMT,ISYMI,ISYMA,ISYMU,                    &
     &                      II,IUTOT,ERI2,SCR)
                  DO IA=1,NSSH(ISYMA)
                    IATOT=IA+NISH(ISYMA)+NASH(ISYMA)
                    ONEADD=Zero
                    IF(ISYM.EQ.1) THEN
                      FAI=FIMO(NFIMOES+(IATOT*(IATOT-1))/2+II)
                      ONEADD=FAI/DBLE(MAX(1,NACTEL))
                    END IF
                    DO IT=1,NASH(ISYMT)
                      ITABS=IT+NAES(ISYMT)
                      ITTOT=IT+NISH(ISYMT)
                      IWA=KTU(ITABS,IUABS)-NTUES(ISYM)
                      IWI=II+NISH(ISYMI)*(IA-1)+IOFF(ISYMA)
                      IW1=IWA+NAS*(IWI-1)
                      IW2=IW1+NAS1
                      IBUF1=IATOT+NORB(ISYMA)*(ITTOT-1)
                      IBUF2=ITTOT+NORB(ISYMT)*(IATOT-1)
                      WAITU=ERI1(IBUF1)
                      IF(ITABS.EQ.IUABS) WAITU=WAITU+ONEADD
                      GA_Arrays(LW)%A(IW1)=WAITU
                      GA_Arrays(LW)%A(IW2)=ERI2(IBUF2)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
!   Put W on disk.
          ICASE=5
          CALL MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
          Call Deallocate_GA_Array(LW)
      END DO

      END SUBROUTINE MKRHSD
