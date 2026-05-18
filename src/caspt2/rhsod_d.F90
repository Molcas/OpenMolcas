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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************

      SUBROUTINE RHSOD_D(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: One
      USE SUPERINDEX, only: MIA, MAREL, MIREL, MTU, MTREL, KTU
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use caspt2_global, only: FIMO
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NACTEL, NASHT, NSYM, NORB, NASUP, NISUP, &
     &                         NIAES, NTUES, NSSH, NASH, NISH

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOBRA1(8,8), IOKET1(8,8), IOBRA2(8,8),          &
     &                  IOKET2(8,8)
      real(kind=wp), ALLOCATABLE:: BRABUF1(:), KETBUF1(:),              &
     &                      BRABUF2(:), KETBUF2(:)
      real(kind=wp) ACTINV, AJTV, AVTJ, FAJ, ONEADD
      integer(kind=iwp) NAS, NIS, lg_W, IASTA, IAEND, IISTA, IIEND, MW, &
     &                  IA, IAABS, IAEND1, IAEND2, IAJ, IAJTOT, IASTA1, &
     &                  IASTA2, IATOT, iCASE, IDX, IFIMOES,  IJ, IJABS, &
     &                  IOAJ, IOAV, IOFFAJ, IOFFAV, IOFFTJ, IOFFTV,     &
     &                  IOTJ, IOTV, ISYA, ISYJ, ISYM, ISYT, ISYV, IT,   &
     &                  ITABS, ITV, IUABS, IUU, IV, IVABS, NAS1,        &
     &                  NBRABUF1, NBRABUF2, NKETBUF1, NKETBUF2, NV, NW
      real(kind=wp), External :: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp) NFIMOES(8)

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case D'
      END IF

!***********************************************************************
! Case D (5,6):
! D1(tv,aj)=(aj,tv) + FIMO(a,j)*Kron(t,v)/NACTEL
! D2(tv,aj)=(tj,av)
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
      CALL CHOVEC_SIZE(4,NBRABUF1,IOBRA1)
      CALL CHOVEC_SIZE(2,NKETBUF1,IOKET1)

      CALL mma_allocate(BRABUF1,NBRABUF1,LABEL='BRABUF1')
      CALL mma_allocate(KETBUF1,NKETBUF1,LABEL='KETBUF1')

      CALL CHOVEC_READ(4,BRABUF1,NBRABUF1)
      CALL CHOVEC_READ(2,KETBUF1,NKETBUF1)

      CALL CHOVEC_SIZE(3,NBRABUF2,IOBRA2)
      CALL CHOVEC_SIZE(1,NKETBUF2,IOKET2)

      CALL mma_allocate(BRABUF2,NBRABUF2,LABEL='BRABUF2')
      CALL mma_allocate(KETBUF2,NKETBUF2,LABEL='KETBUF2')

      CALL CHOVEC_READ(3,BRABUF2,NBRABUF2)
      CALL CHOVEC_READ(1,KETBUF2,NKETBUF2)

      iCASE=5
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
      ! set up FIMO access
      ACTINV=One/DBLE(MAX(1,NACTEL))
      IFIMOES=0
      DO ISYM=1,NSYM
        NFIMOES(ISYM)=IFIMOES
        IFIMOES=IFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2
      END DO

      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

! cases D1, D2 share the RHS along the tu superindex
        NAS1=NAS/2
        IASTA1=IASTA
        IAEND1=IAEND/2
        IASTA2=IAEND1+1
        IAEND2=IAEND

!***********************************************************************
! inner loop over RHS elements in symmetry ISYM
!***********************************************************************
        DO IAJ=IISTA,IIEND
          IAJTOT=IAJ+NIAES(ISYM)
          IJABS=MIA(1,IAJTOT)
          IAABS=MIA(2,IAJTOT)
          IA  =MAREL(1,IAABS)
          ISYA=MAREL(2,IAABS)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          DO ITV=IASTA1,IAEND1 ! these are always all elements
            ITABS=MTU(1,ITV+NTUES(ISYM))
            IVABS=MTU(2,ITV+NTUES(ISYM))
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integral (aj,tv)
            NV=NVTOT_CHOSYM(Mul(ISYA,ISYJ))
            IOAJ=IA-1+NSSH(ISYA)*(IJ-1)
            IOTV=IT-1+NASH(ISYT)*(IV-1)
            IOFFAJ=1+IOBRA1(ISYA,ISYJ)+NV*IOAJ
            IOFFTV=1+IOKET1(ISYT,ISYV)+NV*IOTV
            AJTV=DDOT_(NV,BRABUF1(IOFFAJ),1,KETBUF1(IOFFTV),1)

! D1(tv,aj)=(aj,tv) + FIMO(a,j)*Kron(t,v)/NACTEL
! integrals only
            IDX=ITV+NAS*(IAJ-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=AJTV
#else
            GA_Arrays(lg_w)%A(IDX)=AJTV
#endif
          END DO
! now, dress with FIMO(a,j), only if T==V, so ISYT==ISYV, so if ISYM==1
          IF (ISYM==1) THEN
            IATOT=IA+NISH(ISYA)+NASH(ISYA)
            FAJ=FIMO(NFIMOES(ISYA)+(IATOT*(IATOT-1))/2+IJ)
            ONEADD=FAJ*ACTINV
            DO IUABS=1,NASHT
              IUU=KTU(IUABS,IUABS)
              IDX=IUU+NAS*(IAJ-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=DBL_MB(MW+IDX-1)+ONEADD
#else
              GA_Arrays(lg_w)%A(IDX)=GA_Arrays(lg_w)%A(IDX)             &
     &                                  +ONEADD
#endif
            END DO
          END IF
          DO ITV=IASTA2,IAEND2 ! these are always all elements
            ITABS=MTU(1,ITV-NAS1+NTUES(ISYM))
            IVABS=MTU(2,ITV-NAS1+NTUES(ISYM))
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integral (av,tj)
            NV=NVTOT_CHOSYM(Mul(ISYA,ISYV))
            IOAV=IA-1+NSSH(ISYA)*(IV-1)
            IOTJ=IT-1+NASH(ISYT)*(IJ-1)
            IOFFAV=1+IOBRA2(ISYA,ISYV)+NV*IOAV
            IOFFTJ=1+IOKET2(ISYT,ISYJ)+NV*IOTJ
            AVTJ=DDOT_(NV,BRABUF2(IOFFAV),1,KETBUF2(IOFFTJ),1)

! D2(tv,aj)=(av,tj) + FIMO(a,j)*Kron(t,v)/NACTEL
            IDX=ITV+NAS*(IAJ-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=AVTJ
#else
            GA_Arrays(lg_w)%A(IDX)=AVTJ
#endif
          END DO
        END DO
!***********************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)

      END DO
!***********************************************************************

      CALL mma_deallocate(BRABUF1)
      CALL mma_deallocate(KETBUF1)

      CALL mma_deallocate(BRABUF2)
      CALL mma_deallocate(KETBUF2)

!***********************************************************************

      END SUBROUTINE RHSOD_D
