************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Steven Vancoillie                                      *
************************************************************************
*SVC: compute RHS elements "on demand". If we have access to all the
* Cholesky vectors, we can just instruct a process to compute it's own
* block of RHS elements, computing the integrals directly. This is much
* more computationally intensive, but should scale much better since we
* go from a badly scaling scatter algorithm to no communication at all.
* This also eliminates the need for the GA library in creating the RHS.

* This is a special optimized version for non-symmetric molecules, as
* this allows for convenient sub-blocking of cholesky vectors.

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_NOSYM(IVEC)
      use definitions, only: iwp
#ifdef _DEBUGPRINT_
      use definitions, only: wp
      use caspt2_module, only: nASup, nISup, nSym
#endif
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: VERBOSE
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), Intent(In):: IVEC
#ifdef _DEBUGPRINT_
      integer(kind=iwp) iCase, iSym, NAS, NIS, lg_W
      real(kind=wp) DNRM2
      real(kind=wp), external :: RHS_DDot
#endif


      IF (IPRGLB>=VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using special RHS on-demand algorithm,'
        WRITE(6,'(1X,A)') ' optimized for non-symmetric molecules'
      END IF

#ifdef _MOLCAS_MPP_
      IF (.NOT.Is_Real_Par()) THEN
        WRITE(6,'(1X,A)') 'RHSOD_NOSYM: error: '//
     &                    'fake parallel not supported'
        CALL AbEnd()
      END IF
#endif

      CALL RHSOD_A_NOSYM(IVEC)
      CALL RHSOD_B_NOSYM(IVEC)
      CALL RHSOD_C_NOSYM(IVEC)
      CALL RHSOD_D_NOSYM(IVEC)
      CALL RHSOD_E_NOSYM(IVEC)
      CALL RHSOD_F_NOSYM(IVEC)
      CALL RHSOD_G_NOSYM(IVEC)
      CALL RHSOD_H_NOSYM(IVEC)

#ifdef _DEBUGPRINT_
* compute and print RHS fingerprints
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') 'Case','Sym','Fingerprint'
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') '====','===','==========='
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NAS*NIS/=0) THEN
            CALL RHS_ALLO (NAS,NIS,lg_W)
            CALL RHS_READ (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
            WRITE(6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
          END IF
        END DO
      END DO
#endif

      END SUBROUTINE RHSOD_NOSYM


************************************************************************
* SUBROUTINES FOR THE SEPARATE CASES
************************************************************************

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
      use caspt2_module, only: nSym, NTUV, nIsh, nTUVES, nAsh,
     &                         nOrb, nActEl

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOBRA(8,8), IOKET(8,8)
      real(kind=wp), ALLOCATABLE:: BRA(:), KET(:)
      real(kind=wp) ATVXJ, FTJ, TJVX
      real(kind=wp), External:: DDot_
      integer(kind=iwp) iAEnd, iASta, iCase, IDX, IIEnd, IJ, IOFFTJ,
     &                  IOFFVX, ISYJ, iSym, iSYT, iSYV, iSYX, IT, ITABS,
     &                  ITJ, ITTOT, ITVX, ITVXTOT, IV, IVABS,
     &                  iVX, iX, iXABS, lg_W, mW, nAS, nBra, nFIMOES,
     &                  nIS, nKet, nV, nW, IISTA
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case A'
      END IF

************************************************************************
* Case A:
C   RHS(tvx,j)=(tj,vx)+FIMO(t,j)*kron(v,x)/NACTEL
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(1,NBRA,IOBRA)
      CALL CHOVEC_SIZE(2,NKET,IOKET)

      CALL mma_allocate(BRA,NBRA,LABEL='BRA')
      CALL mma_allocate(KET,NKET,LABEL='KET')

      CALL CHOVEC_READ(1,BRA,NBRA)
      CALL CHOVEC_READ(2,KET,NKET)

      ICASE=1
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
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

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
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
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
************************************************************************

      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)

      END SUBROUTINE RHSOD_A_NOSYM

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_C_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use Constants, only: Zero
      USE SUPERINDEX, only: MTUV, MTREL, KTUV
      USE CHOVEC_IO, only: NVTOT_CHOSYM, CHOVEC_SIZE, CHOVEC_READ
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use caspt2_global, only: FIMO
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: nActEl, nAshT, nSym, nASup, nISup,
     &                         NTUVES, nSsh, nAsh, nIsh, NAES, nOrb
      IMPLICIT None
      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOBRA(8,8), IOKET(8,8)
      real(kind=wp), ALLOCATABLE:: BRA(:), KET(:)
      real(kind=wp) AddOne, ATVX, FAT, SUMU
      real(kind=wp), External:: DDot_
      integer(kind=iwp) IA, IAEND, IASTA, IAT, IATOT, iCASE, IDX,
     &                  IIEND, IISTA, IOFFAT, IOFFVX, ISYA, iSym,
     &                  ISYT, ISYV, ISYX, IT, ITABS, ITTOT, ITVV,
     &                  ITVX, ITVXTOT, IUABS, IUUT, IV, IVABS, IVX,
     &                  IX, IXABS, lg_W, mW, NAS, NBRA, NFIMOES, NIS,
     &                  NKET, NV, NW
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case C'
      END IF

************************************************************************
* Case C:
C   RHS(tvx,a)=(at,vx)+(FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NBRA,IOBRA)
      CALL CHOVEC_SIZE(2,NKET,IOKET)

      CALL mma_allocate(BRA,NBRA,LABEL='BRA')
      CALL mma_allocate(KET,NKET,LABEL='KET')

      CALL CHOVEC_READ(3,BRA,NBRA)
      CALL CHOVEC_READ(2,KET,NKET)

      ICASE=4
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      NFIMOES=0
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE) !NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE) !NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Then
           NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2
           Cycle
        End If

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IA=IISTA,IIEND
          ISYA=ISYM
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
! compute integrals (at,vx)
            NV=NVTOT_CHOSYM(Mul(ISYA,ISYT)) ! JSYM=ISYT*ISYI=ISYU*ISYV
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            IVX=IV-1+NASH(ISYV)*(IX-1)
            IOFFAT=1+IOBRA(ISYA,ISYT)+NV*IAT
            IOFFVX=1+IOKET(ISYV,ISYX)+NV*IVX
            ATVX=DDOT_(NV,BRA(IOFFAT),1,KET(IOFFVX),1)

! W(tvx,a) = (at,vx) + (FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
! write element W(tvx,j), only the (at,vx) part
            IDX=ITVX+NAS*(IA-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=ATVX
#else
            GA_Arrays(lg_W)%A(IDX)=ATVX
#endif
          END DO
! now, add in the part with corrections to the integrals
          IATOT=IA+NISH(ISYM)+NASH(ISYM)
          DO IT=1,NASH(ISYM)
            ITTOT=IT+NISH(ISYM)
            FAT=FIMO(NFIMOES+(IATOT*(IATOT-1))/2+ITTOT)
            SUMU=Zero
            ITABS=NAES(ISYM)+IT
            DO IUABS=1,NASHT
              IUUT=KTUV(IUABS,IUABS,ITABS)-NTUVES(ISYM)
              IDX=IUUT+NAS*(IA-IISTA)
#ifdef _MOLCAS_MPP_
              SUMU=SUMU+DBL_MB(MW+IDX-1)
#else
              SUMU=SUMU+GA_Arrays(lg_w)%A(IDX)
#endif
            END DO
            ADDONE=(FAT-SUMU)/DBLE(MAX(1,NACTEL))
            DO IVABS=1,NASHT
              ITVV=KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
              IDX=ITVV+NAS*(IA-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=DBL_MB(MW+IDX-1)+ADDONE
#else
              GA_Arrays(lg_w)%A(IDX)=GA_Arrays(lg_w)%A(IDX)
     &                                  +ADDONE
#endif
            END DO
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
************************************************************************

      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)

      END SUBROUTINE RHSOD_C_NOSYM

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_B_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Half
      USE SUPERINDEX, only: MIGEJ, MIREL, MTGEU, MTREL, MIGTJ, MTGTU
      USE CHOVEC_IO, only: NVTOT_ChoSym, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: nSym, nASup, nISup, NIGEJES, NTGEUES,
     &                         nAsh, NIGTJES, NTGTUES

      IMPLICIT None

      integer(kind=iwp), intent(in) :: iVec

      integer(kind=iwp) IOSYM(8,8)
      real(kind=wp), ALLOCATABLE:: CHOBUF(:)
      real(kind=wp), parameter:: SQRTH=SQRT(Half)
      real(kind=wp) BMTVJL, BPTVJL, SCL, TJVL, TLVJ
      integer(kind=iwp) IAEND, IASTA, iCASE, IDX, IIEND,
     &              IISTA, IJ, IJABS, IJGEL, IJGELTOT, IJGTL, IJGTLTOT,
     &              IL, ILABS, IOFFTJ, IOFFTL, IOFFVJ, IOFFVL, ISYJ,
     &              ISYL, ISYM, ISYT, ISYV, IT, ITABS, ITGEU,
     &              ITGEUTOT, ITGTU, ITGTUTOT, ITJ, ITL, IV, IVABS,
     &              IVJ, IVL, lg_W, MW, NAS, NCHOBUF, NIS, NV, NW
      real(kind=wp), external:: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case B'
      END IF

************************************************************************
* Case B (2,3):
C   Let  W(tv,j,l)=(jt,lv):
C   BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
C   BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(1,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

      CALL CHOVEC_READ(1,CHOBUF,NCHOBUF)

      iCASE=2
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJGEL=IISTA,IIEND
          IJGELTOT=IJGEL+NIGEJES(ISYM)
          IJABS=MIGEJ(1,IJGELTOT)
          ILABS=MIGEJ(2,IJGELTOT)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          IL  =MIREL(1,ILABS)
          ISYL=MIREL(2,ILABS)
          DO ITGEU=IASTA,IAEND ! these are always all elements
            ITGEUTOT=ITGEU+NTGEUES(ISYM)
            ITABS=MTGEU(1,ITGEUTOT)
            IVABS=MTGEU(2,ITGEUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (ajcl) and (alcj)
            NV=NVTOT_CHOSYM(Mul(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVL=IV-1+NASH(ISYV)*(IL-1)
            IOFFTJ=1+IOSYM(ISYT,ISYJ)+NV*ITJ
            IOFFVL=1+IOSYM(ISYV,ISYL)+NV*IVL
            TJVL=DDOT_(NV,CHOBUF(IOFFTJ),1,CHOBUF(IOFFVL),1)

            NV=NVTOT_CHOSYM(Mul(ISYT,ISYL))
            ITL=IT-1+NASH(ISYT)*(IL-1)
            IVJ=IV-1+NASH(ISYV)*(IJ-1)
            IOFFTL=1+IOSYM(ISYT,ISYL)+NV*ITL
            IOFFVJ=1+IOSYM(ISYV,ISYJ)+NV*IVJ
            TLVJ=DDOT_(NV,CHOBUF(IOFFTL),1,CHOBUF(IOFFVJ),1)

! BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
            SCL=Half
            IF (ITABS==IVABS) SCL=SCL*Half
            IF (ILABS==IJABS) SCL=SCL*SQRTH
            BPTVJL=SCL*(TJVL+TLVJ)
! write element HP(ac,jl)
            IDX=ITGEU+NAS*(IJGEL-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=BPTVJL
#else
            GA_Arrays(lg_w)%A(IDX)=BPTVJL
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************



      iCASE=3
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJGTL=IISTA,IIEND
          IJGTLTOT=IJGTL+NIGTJES(ISYM)
          IJABS=MIGTJ(1,IJGTLTOT)
          ILABS=MIGTJ(2,IJGTLTOT)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          IL  =MIREL(1,ILABS)
          ISYL=MIREL(2,ILABS)
          DO ITGTU=IASTA,IAEND ! these are always all elements
            ITGTUTOT=ITGTU+NTGTUES(ISYM)
            ITABS=MTGTU(1,ITGTUTOT)
            IVABS=MTGTU(2,ITGTUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (tj,vl) and (tlvj)
            NV=NVTOT_CHOSYM(Mul(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVL=IV-1+NASH(ISYV)*(IL-1)
            IOFFTJ=1+IOSYM(ISYT,ISYJ)+NV*ITJ
            IOFFVL=1+IOSYM(ISYV,ISYL)+NV*IVL
            TJVL=DDOT_(NV,CHOBUF(IOFFTJ),1,CHOBUF(IOFFVL),1)

            NV=NVTOT_CHOSYM(Mul(ISYT,ISYL))
            ITL=IT-1+NASH(ISYT)*(IL-1)
            IVJ=IV-1+NASH(ISYV)*(IJ-1)
            IOFFTL=1+IOSYM(ISYT,ISYL)+NV*ITL
            IOFFVJ=1+IOSYM(ISYV,ISYJ)+NV*IVJ
            TLVJ=DDOT_(NV,CHOBUF(IOFFTL),1,CHOBUF(IOFFVJ),1)

! BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
            SCL=Half
            !IF (ITABS==IVABS) SCL=SCL*0.5D0
            !IF (ILABS==IJABS) SCL=SCL*SQRTH
            BMTVJL=SCL*(TJVL-TLVJ)
! write element BM(tv,jl)
            IDX=ITGTU+NAS*(IJGTL-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=BMTVJL
#else
            GA_Arrays(lg_w)%A(IDX)=BMTVJL
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************

      CALL mma_deallocate(CHOBUF)

      END SUBROUTINE RHSOD_B_NOSYM

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
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
      use caspt2_module, only: NSYM, NASUP, NISUP, NAGEBES, NTGEUES,
     &                         NSSH, NAGTBES, NTGTUES

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOSYM(8,8)
      real(kind=wp), ALLOCATABLE:: CHOBUF(:)
      real(kind=wp), Parameter:: SQRTH=SQRT(Half)
      real(kind=wp) ATCV, AVCT, FMTVAC, FPTVAC, SCL
      integer(kind=iwp) IA, IAABS, IAEND, IASTA, IIEND, IISTA, MW,
     &                  IAGEB, IAGEBTOT, IAGTB, IAGTBTOT, IAT, IAV, IC,
     &                  ICABS, iCASE, ICT, ICV, IDX, IOFFAT, IOFFAV,
     &                  IOFFCT, IOFFCV, ISYA, ISYC, ISYM, ISYT, ISYV,
     &                  IT, ITABS, ITGEU, ITGEUTOT, ITGTU, ITGTUTOT,
     &                  IV, IVABS, lg_W, NAS, NCHOBUF, NIS, NV, NW
      real(kind=wp), external:: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case F'
      END IF

************************************************************************
* Case F (8,9):
C FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
C FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,Label='CHOBUF')

      CALL CHOVEC_READ(3,CHOBUF,NCHOBUF)

      iCASE=8
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
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
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************



      iCASE=9
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
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
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************

      CALL mma_deallocate(CHOBUF)

      END SUBROUTINE RHSOD_F_NOSYM

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_H_NOSYM(IVEC)
      use definitions, only: iwp, wp
      use constants, only: Zero, Half, One, Three
      USE SUPERINDEX, only: MIGEJ, MAGEB, MIGTJ, MAGTB
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NSSHT, NAGEB, NIGEJ, NAGTB, NIGTJ

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOSYM(8,8)
      real(kind=wp), ALLOCATABLE:: CHOBUF(:)
      real(kind=wp) AJCL, ALCJ, HMACJL, HPACJL, SCL
      integer(kind=iwp) IA, IASTA, IAEND, IISTA, IIEND, MW, IAGEB,
     &                  IAGTB, IC, iCASE, IDX, IJ, IJGEL, IJGTL, IL,
     &                  lg_W, LIJOFF, LILOFF, NAS, NBLOCK, NCHOBUF,
     &                  NIS, NV, NW
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), PARAMETER :: NOSYM = 1
      real(kind=wp), ALLOCATABLE :: AIBJ(:,:)
      real(kind=wp), parameter:: SQRT3=SQRT(Three), SQRTH=SQRT(Half)

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case H'
      END IF

************************************************************************
* Case H:
C   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
C   WM(jl,ac)=((ajcl)-(alcj))*SQRT(3.0D0)
************************************************************************

      NV=NVTOT_CHOSYM(NOSYM)
      Call mma_ALLOCATE(AIBJ,NSSHT,NSSHT,Label='AIBJ')
      NBLOCK=NV*NSSHT

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,Label='CHOBUF')

      CALL CHOVEC_READ(4,CHOBUF,NCHOBUF)

      iCASE=12

      NAS=NAGEB(NOSYM)
      NIS=NIGEJ(NOSYM)
      NW=NAS*NIS

      IF(NW/=0) Then

      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
      NW=NAS*(IIEND-IISTA+1)

      DO IJGEL=IISTA,IIEND
        IJ=MIGEJ(1,IJGEL)
        IL=MIGEJ(2,IJGEL)
        LIJOFF=1+NBLOCK*(IJ-1)
        LILOFF=1+NBLOCK*(IL-1)
        ! precompute integral blocks
        CALL DGEMM_('T','N',NSSHT,NSSHT,NV,
     &              One,CHOBUF(LIJOFF),NV,CHOBUF(LILOFF),NV,
     &              Zero,AIBJ,NSSHT)
        DO IAGEB=IASTA,IAEND ! these are always all elements
          IA=MAGEB(1,IAGEB)
          IC=MAGEB(2,IAGEB)
          AJCL=AIBJ(IA,IC)
          ALCJ=AIBJ(IC,IA)
! HP(ac,jl)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
          SCL=One
          IF (IA==IC) SCL=SCL*SQRTH
          IF (IL==IJ) SCL=SCL*SQRTH
          HPACJL=SCL*(AJCL+ALCJ)
! write element HP(ac,jl)
          IDX=IAGEB+NAS*(IJGEL-IISTA)
#ifdef _MOLCAS_MPP_
          DBL_MB(MW+IDX-1)=HPACJL
#else
          GA_Arrays(lg_w)%A(IDX)=HPACJL
#endif
        END DO
      END DO
************************************************************************
      CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
      CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,NOSYM,iVEC)
      CALL RHS_FREE (lg_W)
************************************************************************
      End If

      iCASE=13

      NAS=NAGTB(NOSYM)
      NIS=NIGTJ(NOSYM)
      NW=NAS*NIS

      IF(NW/=0) Then

      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
      NW=NAS*(IIEND-IISTA+1)

      DO IJGTL=IISTA,IIEND
        IJ=MIGTJ(1,IJGTL)
        IL=MIGTJ(2,IJGTL)
        LIJOFF=1+NBLOCK*(IJ-1)
        LILOFF=1+NBLOCK*(IL-1)
        ! precompute integral blocks
        CALL DGEMM_('T','N',NSSHT,NSSHT,NV,
     &              One,CHOBUF(LIJOFF),NV,CHOBUF(LILOFF),NV,
     &              Zero,AIBJ,NSSHT)
        DO IAGTB=IASTA,IAEND ! these are always all elements
          IA=MAGTB(1,IAGTB)
          IC=MAGTB(2,IAGTB)
          AJCL=AIBJ(IA,IC)
          ALCJ=AIBJ(IC,IA)
! HM(ac,jl)=((ajcl)-(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
          SCL=SQRT3
          HMACJL=SCL*(AJCL-ALCJ)
! write element HM(ac,jl)
          IDX=IAGTB+NAS*(IJGTL-IISTA)
#ifdef _MOLCAS_MPP_
          DBL_MB(MW+IDX-1)=HMACJL
#else
          GA_Arrays(lg_W)%A(IDX)=HMACJL
#endif
        END DO
      END DO
************************************************************************

      CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
      CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,NOSYM,iVEC)
      CALL RHS_FREE (lg_W)
      End If
************************************************************************

      CALL mma_deallocate(CHOBUF)

      call mma_DEALLOCATE(AIBJ)

      END SUBROUTINE RHSOD_H_NOSYM


*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_D_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: One
      USE SUPERINDEX, only: MIA, MAREL, MIREL, MTU, MTREL, KTU
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use caspt2_global, only: FIMO
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NACTEL, NASHT, NSYM, NORB, NASUP, NISUP,
     &                         NIAES, NTUES, NSSH, NASH, NISH, NISH

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOBRA1(8,8), IOKET1(8,8), IOBRA2(8,8),
     &                  IOKET2(8,8)
      real(kind=wp), ALLOCATABLE:: BRABUF1(:), KETBUF1(:),
     &                      BRABUF2(:), KETBUF2(:)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp) NFIMOES(8)
      real(kind=wp) ACTINV, AJTV, AVTJ, FAJ, ONEADD
      integer(kind=iwp) IA, IAABS, NAS, NIS, lg_W, IASTA, IAEND, IISTA,
     &                  IIEND, MW, IAEND1, IAEND2, IAJ, IAJTOT, IASTA1,
     &                  IASTA2, IATOT, iCASE, IDX, IFIMOES, IJ, IJABS,
     &                  IOAJ, IOAV, IOFFAJ, IOFFAV, IOFFTJ, IOFFTV,
     &                  IOTJ, IOTV, ISYA, ISYJ, ISYM, ISYT, ISYV, IT,
     &                  ITABS, ITV, IUABS, IUU, IV, IVABS, NAS1,
     &                  NBRABUF1, NBRABUF2, NKETBUF1, NKETBUF2, NV, NW
      real(kind=wp), external:: DDot_

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case D'
      END IF

************************************************************************
* Case D (5,6):
C D1(tv,aj)=(aj,tv) + FIMO(a,j)*Kron(t,v)/NACTEL
C D2(tv,aj)=(tj,av)
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
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
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
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

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
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
              GA_Arrays(lg_w)%A(IDX)=GA_Arrays(lg_w)%A(IDX)
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
            GA_Arrays(lg_W)%A(IDX)=AVTJ
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)

      END DO
************************************************************************

      CALL mma_deallocate(BRABUF1)
      CALL mma_deallocate(KETBUF1)

      CALL mma_deallocate(BRABUF2)
      CALL mma_deallocate(KETBUF2)

************************************************************************

      END SUBROUTINE RHSOD_D_NOSYM

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_E_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use Constants, only: Half, OneHalf
      USE SUPERINDEX, only: MIGEJ, MIREL, MIGTJ, MIREL
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NSYM, NASUP, NISUP, NSSH, NIGEJ,
     &                         NASH, NIGEJES, NIGTJ, NIGTJES

      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC

      integer(kind=iwp) IOBRA(8,8), IOKET(8,8)
      real(kind=wp), ALLOCATABLE:: BRABUF(:), KETBUF(:)
      real(kind=wp), parameter:: SQRTH=SQRT(Half), SQRTA=SQRT(OneHalf)
      real(kind=wp) AJVL, ALVJ, EM, EP, SCL
      integer(kind=iwp) IA, NAS, NIS, lg_W, IASTA, IAEND, IISTA, IIEND,
     &                  MW, IAJ, IAJGEL, IAJGELEND, IAJGELSTA, IAJGTL,
     &                  IAJGTLEND, IAJGTLSTA, IAL, iCASE, IDX, IJ,
     &                  IJABS, IJGEL, IJGELTOT, IJGTL, IJGTLTOT, IL,
     &                  ILABS, IOFF, IOFFAJ, IOFFAL, IOFFVJ, IOFFVL,
     &                  ISYA, ISYJ, ISYJL, ISYL, ISYM, ISYV, IV, IVJ,
     &                  IVL, NA, NBRABUF, NJL, NKETBUF, NV, NW
      real(kind=wp), External:: DDot_
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case E'
      END IF

************************************************************************
* Case E (6,7):
C EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
C EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
************************************************************************

* -SVC- Case E is slightly special, in that the inactive superindices are
* so large, that it is suboptimal to have a direct translation table for
* them. Instead, the code loops over symmetry blocks of A-JL and figures
* out if the indices on the processor fall within a block or not. Within
* a A-JL symmetry block, NA(ISYA) and NIGEJ(ISYJL) are known, so they can
* be determined by integer division. This could be optimized by combining
* it with loop peeling (on the todo list?).


************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NBRABUF,IOBRA)
      CALL CHOVEC_SIZE(1,NKETBUF,IOKET)

      CALL mma_allocate(BRABUF,NBRABUF,LABEL='BRABUF')
      CALL mma_allocate(KETBUF,NKETBUF,LABEL='KETBUF')

      CALL CHOVEC_READ(4,BRABUF,NBRABUF)
      CALL CHOVEC_READ(1,KETBUF,NKETBUF)

      iCASE=6
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYA=1,NSYM
          ISYJL=Mul(ISYA,ISYM)
          ISYV=ISYM

          NA=NSSH(ISYA)
          NJL=NIGEJ(ISYJL)
          ! what is start/end in this block?
          IAJGELSTA=MAX(IISTA-IOFF,1)
          IAJGELEND=MIN(IIEND-IOFF,NA*NJL)

          DO IAJGEL=IAJGELSTA,IAJGELEND
            IJGEL=(IAJGEL-1)/NA+1
            IA=IAJGEL-NA*(IJGEL-1)
            IJGELTOT=IJGEL+NIGEJES(ISYJL)
            IJABS=MIGEJ(1,IJGELTOT)
            ILABS=MIGEJ(2,IJGELTOT)
            IJ  =MIREL(1,IJABS)
            ISYJ=MIREL(2,IJABS)
            IL  =MIREL(1,ILABS)
            ISYL=MIREL(2,ILABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(Mul(ISYA,ISYJ))
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IVL=IV-1+NASH(ISYV)*(IL-1)
              IOFFAJ=1+IOBRA(ISYA,ISYJ)+NV*IAJ
              IOFFVL=1+IOKET(ISYV,ISYL)+NV*IVL
              AJVL=DDOT_(NV,BRABUF(IOFFAJ),1,KETBUF(IOFFVL),1)

              NV=NVTOT_CHOSYM(Mul(ISYA,ISYL))
              IAL=IA-1+NSSH(ISYA)*(IL-1)
              IVJ=IV-1+NASH(ISYV)*(IJ-1)
              IOFFAL=1+IOBRA(ISYA,ISYL)+NV*IAL
              IOFFVJ=1+IOKET(ISYV,ISYJ)+NV*IVJ
              ALVJ=DDOT_(NV,BRABUF(IOFFAL),1,KETBUF(IOFFVJ),1)

! EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
              IF (ILABS==IJABS) THEN
                SCL=Half
              ELSE
                SCL=SQRTH
              END IF
              EP=SCL*(AJVL+ALVJ)
! write element EP
              IDX=IV+NAS*(IAJGEL+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=EP
#else
              GA_Arrays(lg_w)%A(IDX)=EP
#endif
            END DO
          END DO

          IOFF=IOFF+NA*NJL
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************



      iCASE=7
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYA=1,NSYM
          ISYJL=Mul(ISYA,ISYM)
          ISYV=ISYM

          NA=NSSH(ISYA)
          NJL=NIGTJ(ISYJL)
          ! what is start/end in this block?
          IAJGTLSTA=MAX(IISTA-IOFF,1)
          IAJGTLEND=MIN(IIEND-IOFF,NA*NJL)

          DO IAJGTL=IAJGTLSTA,IAJGTLEND
            IJGTL=(IAJGTL-1)/NA+1
            IA=IAJGTL-NA*(IJGTL-1)
            IJGTLTOT=IJGTL+NIGTJES(ISYJL)
            IJABS=MIGTJ(1,IJGTLTOT)
            ILABS=MIGTJ(2,IJGTLTOT)
            IJ  =MIREL(1,IJABS)
            ISYJ=MIREL(2,IJABS)
            IL  =MIREL(1,ILABS)
            ISYL=MIREL(2,ILABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(Mul(ISYA,ISYJ))
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IVL=IV-1+NASH(ISYV)*(IL-1)
              IOFFAJ=1+IOBRA(ISYA,ISYJ)+NV*IAJ
              IOFFVL=1+IOKET(ISYV,ISYL)+NV*IVL
              AJVL=DDOT_(NV,BRABUF(IOFFAJ),1,KETBUF(IOFFVL),1)

              NV=NVTOT_CHOSYM(Mul(ISYA,ISYL))
              IAL=IA-1+NSSH(ISYA)*(IL-1)
              IVJ=IV-1+NASH(ISYV)*(IJ-1)
              IOFFAL=1+IOBRA(ISYA,ISYL)+NV*IAL
              IOFFVJ=1+IOKET(ISYV,ISYJ)+NV*IVJ
              ALVJ=DDOT_(NV,BRABUF(IOFFAL),1,KETBUF(IOFFVJ),1)

! EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
              EM=SQRTA*(AJVL-ALVJ)
! write element EM
              IDX=IV+NAS*(IAJGTL+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=EM
#else
              GA_Arrays(lg_w)%A(IDX)=EM
#endif
            END DO
          END DO

          IOFF=IOFF+NA*NJL
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************

      CALL mma_deallocate(BRABUF)
      CALL mma_deallocate(KETBUF)

      END SUBROUTINE RHSOD_E_NOSYM

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_G_NOSYM(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Half, OneHalf
      USE SUPERINDEX, only: MAGEB, MAREL, MAGTB
      USE CHOVEC_IO, only: NVTOT_CHOSYM, ChoVec_Size, ChoVec_Read
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      use caspt2_module, only: NSYM, NASUP, NISUP, NISH, NAGEB,
     &                         NAGEBES, NSSH, NAGTB, NAGTBES

      IMPLICIT None

      integer(iwp), intent(in):: IVEC

      integer(iwp) IOBRA(8,8), IOKET(8,8)
      real(kind=wp), ALLOCATABLE:: BRABUF(:), KETBUF(:)
      real(kind=wp), parameter:: SQRTH=SQRT(Half), SQRTA=SQRT(OneHalf)
      real(kind=wp) AVCJ, CVAJ, GM, GP, SCL
      integer(iwp) IA, IAABS, NAS, NIS, lg_W, IASTA, IAEND, IISTA,
     &             IIEND, MW, IAGEC, IAGECTOT, IAGTC, IAGTCTOT, IAJ,
     &             IAV, IC, ICABS, iCASE, ICJ, ICV, IDX, IJ, IJAGEC,
     &             IJAGECEND, IJAGECSTA, IJAGTC, IJAGTCEND, IJAGTCSTA,
     &             IOFF, IOFFAJ, IOFFAV, IOFFCJ, IOFFCV, ISYA, ISYAC,
     &             ISYC, ISYJ, ISYM, ISYV, IV, NAC, NBRABUF, NJ,
     &             NKETBUF, NV, NW
      real(kind=wp), external:: DDOT_
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb>=DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case G'
      END IF

************************************************************************
* Case G (10,11):
C GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
C GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
************************************************************************

* -SVC- Case G is slightly special, in that the inactive superindices are
* so large, that it is suboptimal to have a direct translation table for
* them. Instead, the code loops over symmetry blocks of J-AC and figures
* out if the indices on the processor fall within a block or not. Within
* a J-AC symmetry block, NJ(ISYJ) and NAGEB(ISYAC) are known, so they can
* be determined by integer division. This could be optimized by combining
* it with loop peeling (on the todo list?).

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NBRABUF,IOBRA)
      CALL CHOVEC_SIZE(4,NKETBUF,IOKET)

      CALL mma_allocate(BRABUF,NBRABUF,LABEL='BRABUF')
      CALL mma_allocate(KETBUF,NKETBUF,LABEL='KETBUF')

      CALL CHOVEC_READ(3,BRABUF,NBRABUF)
      CALL CHOVEC_READ(4,KETBUF,NKETBUF)

      iCASE=10
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYJ=1,NSYM
          ISYAC=Mul(ISYJ,ISYM)
          ISYV=ISYM

          NJ=NISH(ISYJ)
          NAC=NAGEB(ISYAC)
          ! what is start/end in this block?
          IJAGECSTA=MAX(IISTA-IOFF,1)
          IJAGECEND=MIN(IIEND-IOFF,NJ*NAC)

          DO IJAGEC=IJAGECSTA,IJAGECEND
            IAGEC=(IJAGEC-1)/NJ+1
            IJ=IJAGEC-NJ*(IAGEC-1)
            IAGECTOT=IAGEC+NAGEBES(ISYAC)
            IAABS=MAGEB(1,IAGECTOT)
            ICABS=MAGEB(2,IAGECTOT)
            IA  =MAREL(1,IAABS)
            ISYA=MAREL(2,IAABS)
            IC  =MAREL(1,ICABS)
            ISYC=MAREL(2,ICABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(Mul(ISYA,ISYV))
              IAV=IA-1+NSSH(ISYA)*(IV-1)
              ICJ=IC-1+NSSH(ISYC)*(IJ-1)
              IOFFAV=1+IOBRA(ISYA,ISYV)+NV*IAV
              IOFFCJ=1+IOKET(ISYC,ISYJ)+NV*ICJ
              AVCJ=DDOT_(NV,BRABUF(IOFFAV),1,KETBUF(IOFFCJ),1)

              NV=NVTOT_CHOSYM(Mul(ISYC,ISYV))
              ICV=IC-1+NSSH(ISYC)*(IV-1)
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IOFFCV=1+IOBRA(ISYC,ISYV)+NV*ICV
              IOFFAJ=1+IOKET(ISYA,ISYJ)+NV*IAJ
              CVAJ=DDOT_(NV,BRABUF(IOFFCV),1,KETBUF(IOFFAJ),1)

C GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
              IF (IAABS==ICABS) THEN
                SCL=Half
              ELSE
                SCL=SQRTH
              END IF
              GP=SCL*(AVCJ+CVAJ)
! write element EP
              IDX=IV+NAS*(IJAGEC+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=GP
#else
              GA_Arrays(lg_w)%A(IDX)=GP
#endif
            END DO
          END DO

          IOFF=IOFF+NJ*NAC
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************



      iCASE=11
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW==0) Cycle

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYJ=1,NSYM
          ISYAC=Mul(ISYJ,ISYM)
          ISYV=ISYM

          NJ=NISH(ISYJ)
          NAC=NAGTB(ISYAC)
          ! what is start/end in this block?
          IJAGTCSTA=MAX(IISTA-IOFF,1)
          IJAGTCEND=MIN(IIEND-IOFF,NJ*NAC)

          DO IJAGTC=IJAGTCSTA,IJAGTCEND
            IAGTC=(IJAGTC-1)/NJ+1
            IJ=IJAGTC-NJ*(IAGTC-1)
            IAGTCTOT=IAGTC+NAGTBES(ISYAC)
            IAABS=MAGTB(1,IAGTCTOT)
            ICABS=MAGTB(2,IAGTCTOT)
            IA  =MAREL(1,IAABS)
            ISYA=MAREL(2,IAABS)
            IC  =MAREL(1,ICABS)
            ISYC=MAREL(2,ICABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(Mul(ISYA,ISYV))
              IAV=IA-1+NSSH(ISYA)*(IV-1)
              ICJ=IC-1+NSSH(ISYC)*(IJ-1)
              IOFFAV=1+IOBRA(ISYA,ISYV)+NV*IAV
              IOFFCJ=1+IOKET(ISYC,ISYJ)+NV*ICJ
              AVCJ=DDOT_(NV,BRABUF(IOFFAV),1,KETBUF(IOFFCJ),1)

              NV=NVTOT_CHOSYM(Mul(ISYC,ISYV))
              ICV=IC-1+NSSH(ISYC)*(IV-1)
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IOFFCV=1+IOBRA(ISYC,ISYV)+NV*ICV
              IOFFAJ=1+IOKET(ISYA,ISYJ)+NV*IAJ
              CVAJ=DDOT_(NV,BRABUF(IOFFCV),1,KETBUF(IOFFAJ),1)

C GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
              GM=SQRTA*(AVCJ-CVAJ)
! write element GM
              IDX=IV+NAS*(IJAGTC+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=GM
#else
              GA_Arrays(lg_w)%A(IDX)=GM
#endif
            END DO
          END DO

          IOFF=IOFF+NJ*NAC
        END DO
************************************************************************

        CALL RHS_Release_Update (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
      END DO
************************************************************************

      CALL mma_deallocate(BRABUF)
      CALL mma_deallocate(KETBUF)

      END SUBROUTINE RHSOD_G_NOSYM
