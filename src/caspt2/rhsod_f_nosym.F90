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

      SUBROUTINE RHSOD_F_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Half
      USE SUPERINDEX, only: MAGEB, MAREL, MTGEU, MTREL, MAGTB, MTGTU
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NSYM, NASUP, NISUP, NAGEBES, NTGEUES,    &
     &                         NSSH, NAGTBES, NTGTUES

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOSYM(8,8)
      real(kind=wp), ALLOCATABLE:: CHOBUF(:)
      real(kind=wp), Parameter:: SQRTH=SQRT(Half)
      real(kind=wp) ATCV, AVCT, FMTVAC, FPTVAC, SCL
      integer(kind=iwp) IA, IAABS, IAEND, IASTA, IIEND, IISTA, MW,      &
     &                  IAGEB, IAGEBTOT, IAGTB, IAGTBTOT, IAT, IAV, IC, &
     &                  ICABS, iCASE, ICT, ICV, IDX, IOFFAT, IOFFAV,    &
     &                  IOFFCT, IOFFCV, ISYA, ISYC, ISYM, ISYT, ISYV,   &
     &                  IT, ITABS, ITGEU, ITGEUTOT, ITGTU, ITGTUTOT,    &
     &                  IV, IVABS, lg_W, NAS, NCHOBUF, NIS, NV, NW
      real(kind=wp), external:: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case F'
      END IF

!***********************************************************************
! Case F (8,9):
! FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
! FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
      CALL CHOVEC_SIZE(3,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,Label='CHOBUF')

      CALL CHOVEC_READ(3,CHOBUF,NCHOBUF)

      iCASE=8
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

!***********************************************************************
! inner loop over RHS elements in symmetry ISYM
!***********************************************************************
        DO IAGEB=IISTA,IIEND
          IAGEBTOT=IAGEB+NAGEBES(ISYM)
          IAABS=MAGEB(1,IAGEBTOT)
          ICABS=MAGEB(2,IAGEBTOT)
          IA  =MAREL(1,IAABS)
          ISYA=MAREL(2,IAABS)
          IC  =MAREL(1,ICABS)
          ISYC=MAREL(2,ICABS)
          DO ITGEU=IASTA,IAEND ! these are always all elements
            ITGEUTOT=ITGEU+NTGEUES(ISYM)
            ITABS=MTGEU(1,ITGEUTOT)
            IVABS=MTGEU(2,ITGEUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (ta,vc) and (tc,va)
            NV=NVTOT_CHOSYM(Mul(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            ICV=IC-1+NSSH(ISYC)*(IV-1)
            IOFFAT=1+IOSYM(ISYA,ISYT)+NV*IAT
            IOFFCV=1+IOSYM(ISYC,ISYV)+NV*ICV
            ATCV=DDOT_(NV,CHOBUF(IOFFAT),1,CHOBUF(IOFFCV),1)

            NV=NVTOT_CHOSYM(Mul(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAV=IA-1+NSSH(ISYA)*(IV-1)
            ICT=IC-1+NSSH(ISYC)*(IT-1)
            IOFFAV=1+IOSYM(ISYA,ISYV)+NV*IAV
            IOFFCT=1+IOSYM(ISYC,ISYT)+NV*ICT
            AVCT=DDOT_(NV,CHOBUF(IOFFAV),1,CHOBUF(IOFFCT),1)

! FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
            SCL=Half
            IF (ITABS==IVABS) SCL=SCL*Half
            IF (IAABS==ICABS) SCL=SCL*SQRTH
            FPTVAC=SCL*(ATCV+AVCT)
! write element FP(tv,ac)
            IDX=ITGEU+NAS*(IAGEB-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=FPTVAC
#else
            GA_Arrays(lg_W)%A(IDX)=FPTVAC
#endif
          END DO
        END DO
!***********************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
!***********************************************************************



      iCASE=9
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

!***********************************************************************
! inner loop over RHS elements in symmetry ISYM
!***********************************************************************
        DO IAGTB=IISTA,IIEND
          IAGTBTOT=IAGTB+NAGTBES(ISYM)
          IAABS=MAGTB(1,IAGTBTOT)
          ICABS=MAGTB(2,IAGTBTOT)
          IA  =MAREL(1,IAABS)
          ISYA=MAREL(2,IAABS)
          IC  =MAREL(1,ICABS)
          ISYC=MAREL(2,ICABS)
          DO ITGTU=IASTA,IAEND ! these are always all elements
            ITGTUTOT=ITGTU+NTGTUES(ISYM)
            ITABS=MTGTU(1,ITGTUTOT)
            IVABS=MTGTU(2,ITGTUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (at,cv) and (av,ct)
            NV=NVTOT_CHOSYM(Mul(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            ICV=IC-1+NSSH(ISYC)*(IV-1)
            IOFFAT=1+IOSYM(ISYA,ISYT)+NV*IAT
            IOFFCV=1+IOSYM(ISYC,ISYV)+NV*ICV
            ATCV=DDOT_(NV,CHOBUF(IOFFAT),1,CHOBUF(IOFFCV),1)

            NV=NVTOT_CHOSYM(Mul(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAV=IA-1+NSSH(ISYA)*(IV-1)
            ICT=IC-1+NSSH(ISYC)*(IT-1)
            IOFFAV=1+IOSYM(ISYA,ISYV)+NV*IAV
            IOFFCT=1+IOSYM(ISYC,ISYT)+NV*ICT
            AVCT=DDOT_(NV,CHOBUF(IOFFAV),1,CHOBUF(IOFFCT),1)

! FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
            SCL=Half
            !IF (ITABS==IVABS) SCL=SCL*0.5D0
            !IF (IAABS==ICABS) SCL=SCL*SQRTH
            FMTVAC=SCL*(AVCT-ATCV)
! write element FM(tv,ac)
            IDX=ITGTU+NAS*(IAGTB-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=FMTVAC
#else
            GA_Arrays(lg_w)%A(IDX)=FMTVAC
#endif
          END DO
        END DO
!***********************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
!***********************************************************************

      CALL mma_deallocate(CHOBUF)

      END SUBROUTINE RHSOD_F_NOSYM
