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

      SUBROUTINE RHSOD_A_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      USE SUPERINDEX, only: MTUV, MTREL
      USE CHOVEC_IO, only: NVTOT_ChoSym, ChoVec_Size, ChoVec_read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use caspt2_global, only: FIMO
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: nSym, NTUV, nIsh, nTUVES, nAsh,          &
     &                         nOrb, nActEl

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOBRA(8,8), IOKET(8,8)
      real(kind=wp), ALLOCATABLE:: BRA(:), KET(:)
      real(kind=wp) ATVXJ, FTJ, TJVX
      real(kind=wp), External:: DDot_
      integer(kind=iwp) iAEnd, iASta, iCase, IDX, IIEnd, IJ, IOFFTJ,    &
     &                  IOFFVX, ISYJ, iSym, iSYT, iSYV, iSYX, IT, ITABS,&
     &                  ITJ, ITTOT, ITVX, ITVXTOT, IV, IVABS,           &
     &                  iVX, iX, iXABS, lg_W, mW, nAS, nBra, nFIMOES,   &
     &                  nIS, nKet, nV, nW, IISTA
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case A'
      END IF

!***********************************************************************
! Case A:
!   RHS(tvx,j)=(tj,vx)+FIMO(t,j)*kron(v,x)/NACTEL
!***********************************************************************

!***********************************************************************
!SVC: read in all the cholesky vectors (need all symmetries)
!***********************************************************************
      CALL CHOVEC_SIZE(1,NBRA,IOBRA)
      CALL CHOVEC_SIZE(2,NKET,IOKET)

      CALL mma_allocate(BRA,NBRA,LABEL='BRA')
      CALL mma_allocate(KET,NKET,LABEL='KET')

      CALL CHOVEC_READ(1,BRA,NBRA)
      CALL CHOVEC_READ(2,KET,NKET)

      ICASE=1
!***********************************************************************
! outer loop over symmetry blocks in the RHS
!***********************************************************************
      NFIMOES=0
      DO ISYM=1,NSYM

        NAS=NTUV(ISYM) !NASUP(ISYM,ICASE)
        NIS=NISH(ISYM) !NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Then
          NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2
          Cycle
        End If

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)

!***********************************************************************
! inner loop over RHS elements in symmetry ISYM
!***********************************************************************
        DO IJ=IISTA,IIEND
          ISYJ=ISYM
          DO ITVX=IASTA,IAEND ! these are always all elements
            ITVXTOT=ITVX+NTUVES(ISYM)
            ITABS=MTUV(1,ITVXTOT)
            IVABS=MTUV(2,ITVXTOT)
            IXABS=MTUV(3,ITVXTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
            IX  =MTREL(1,IXABS)
            ISYX=MTREL(2,IXABS)
! compute integrals (tiuv)
            NV=NVTOT_CHOSYM(Mul(ISYT,ISYJ)) ! JSYM=ISYT*ISYI=ISYU*ISYV
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVX=IV-1+NASH(ISYV)*(IX-1)
            IOFFTJ=1+IOBRA(ISYT,ISYJ)+NV*ITJ
            IOFFVX=1+IOKET(ISYV,ISYX)+NV*IVX
            TJVX=DDOT_(NV,BRA(IOFFTJ),1,KET(IOFFVX),1)
! A(tvx,j) = (tjvx) + FIMO(t,j)*delta(v,x)/NACTEL
            IF (ISYT==ISYJ.AND.IVABS==IXABS) THEN
              ITTOT=IT+NISH(ISYT)
              FTJ=FIMO(NFIMOES+(ITTOT*(ITTOT-1))/2+IJ)
              ATVXJ=TJVX+FTJ/DBLE(MAX(1,NACTEL))
            ELSE
              ATVXJ=TJVX
            END IF
! write element A(tvx,j)
            IDX=ITVX+NAS*(IJ-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=ATVXJ
#else
            GA_Arrays(lg_W)%A(IDX)=ATVXJ
#endif
          END DO
        END DO
!***********************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
!***********************************************************************

      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)

      END SUBROUTINE RHSOD_A_NOSYM
