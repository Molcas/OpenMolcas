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
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "para_info.fh"


      IF (IPRGLB.GE.VERBOSE) THEN
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
          IF (NAS*NIS.NE.0) THEN
            CALL RHS_ALLO (NAS,NIS,lg_W)
            CALL RHS_READ (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
            WRITE(6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
          END IF
        END DO
      END DO
#endif


      END


************************************************************************
* SUBROUTINES FOR THE SEPARATE CASES
************************************************************************

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      SUBROUTINE RHSOD_A_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOBRA(8,8), IOKET(8,8)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif

      IF (iPrGlb.GE.DEBUG) THEN
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

      CALL GETMEM('BRABUF','ALLO','REAL',LBRA,NBRA)
      CALL GETMEM('KETBUF','ALLO','REAL',LKET,NKET)

      CALL CHOVEC_READ(1,LBRA)
      CALL CHOVEC_READ(2,LKET)

      ICASE=1
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      NFIMOES=0
      DO ISYM=1,NSYM

        NAS=NTUV(ISYM) !NASUP(ISYM,ICASE)
        NIS=NISH(ISYM) !NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 1

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
            NV=NVTOT_CHOSYM(MUL(ISYT,ISYJ)) ! JSYM=ISYT*ISYI=ISYU*ISYV
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVX=IV-1+NASH(ISYV)*(IX-1)
            IOFFTJ=LBRA+IOBRA(ISYT,ISYJ)+NV*ITJ
            IOFFVX=LKET+IOKET(ISYV,ISYX)+NV*IVX
            TJVX=DDOT_(NV,WORK(IOFFTJ),1,WORK(IOFFVX),1)
! A(tvx,j) = (tjvx) + FIMO(t,j)*delta(v,x)/NACTEL
            IF (ISYT.EQ.ISYJ.AND.IVABS.EQ.IXABS) THEN
              ITTOT=IT+NISH(ISYT)
              FTJ=WORK(LFIMO+NFIMOES+(ITTOT*(ITTOT-1))/2+IJ-1)
              ATVXJ=TJVX+FTJ/DBLE(MAX(1,NACTEL))
            ELSE
              ATVXJ=TJVX
            END IF
! write element A(tvx,j)
            IDX=ITVX+NAS*(IJ-IISTA)
            DBL_MB(MW+IDX-1)=ATVXJ
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 1      CONTINUE

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
************************************************************************

      CALL GETMEM('BRABUF','FREE','REAL',LBRA,NBRA)
      CALL GETMEM('KETBUF','FREE','REAL',LKET,NKET)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_C_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOBRA(8,8), IOKET(8,8)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif

      IF (iPrGlb.GE.DEBUG) THEN
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

      CALL GETMEM('BRABUF','ALLO','REAL',LBRA,NBRA)
      CALL GETMEM('KETBUF','ALLO','REAL',LKET,NKET)

      CALL CHOVEC_READ(3,LBRA)
      CALL CHOVEC_READ(2,LKET)

      ICASE=4
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      NFIMOES=0
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE) !NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE) !NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 4

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
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYT)) ! JSYM=ISYT*ISYI=ISYU*ISYV
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            IVX=IV-1+NASH(ISYV)*(IX-1)
            IOFFAT=LBRA+IOBRA(ISYA,ISYT)+NV*IAT
            IOFFVX=LKET+IOKET(ISYV,ISYX)+NV*IVX
            ATVX=DDOT_(NV,WORK(IOFFAT),1,WORK(IOFFVX),1)

! W(tvx,a) = (at,vx) + (FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
! write element W(tvx,j), only the (at,vx) part
            IDX=ITVX+NAS*(IA-IISTA)
            DBL_MB(MW+IDX-1)=ATVX
          END DO
! now, add in the part with corrections to the integrals
          IATOT=IA+NISH(ISYM)+NASH(ISYM)
          DO IT=1,NASH(ISYM)
            ITTOT=IT+NISH(ISYM)
            FAT=WORK(LFIMO+NFIMOES+(IATOT*(IATOT-1))/2+ITTOT-1)
            SUMU=0.0D0
            ITABS=NAES(ISYM)+IT
            DO IUABS=1,NASHT
              IUUT=KTUV(IUABS,IUABS,ITABS)-NTUVES(ISYM)
              IDX=IUUT+NAS*(IA-IISTA)
              SUMU=SUMU+DBL_MB(MW+IDX-1)
            END DO
            ADDONE=(FAT-SUMU)/DBLE(MAX(1,NACTEL))
            DO IVABS=1,NASHT
              ITVV=KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
              IDX=ITVV+NAS*(IA-IISTA)
              DBL_MB(MW+IDX-1)=DBL_MB(MW+IDX-1)+ADDONE
            END DO
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 4      CONTINUE

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
************************************************************************

      CALL GETMEM('BRABUF','FREE','REAL',LBRA,NBRA)
      CALL GETMEM('KETBUF','FREE','REAL',LKET,NKET)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_B_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOSYM(8,8)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case B'
      END IF

************************************************************************
* Case B (2,3):
C   Let  W(tv,j,l)=(jt,lv):
C   BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
C   BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
************************************************************************

      SQRTH=SQRT(0.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(1,NCHOBUF,IOSYM)

      CALL GETMEM('CHOBUF','ALLO','REAL',LCHOBUF,NCHOBUF)

      CALL CHOVEC_READ(1,LCHOBUF)

      iCASE=2
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 2

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
            NV=NVTOT_CHOSYM(MUL(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVL=IV-1+NASH(ISYV)*(IL-1)
            IOFFTJ=LCHOBUF+IOSYM(ISYT,ISYJ)+NV*ITJ
            IOFFVL=LCHOBUF+IOSYM(ISYV,ISYL)+NV*IVL
            TJVL=DDOT_(NV,WORK(IOFFTJ),1,WORK(IOFFVL),1)

            NV=NVTOT_CHOSYM(MUL(ISYT,ISYL))
            ITL=IT-1+NASH(ISYT)*(IL-1)
            IVJ=IV-1+NASH(ISYV)*(IJ-1)
            IOFFTL=LCHOBUF+IOSYM(ISYT,ISYL)+NV*ITL
            IOFFVJ=LCHOBUF+IOSYM(ISYV,ISYJ)+NV*IVJ
            TLVJ=DDOT_(NV,WORK(IOFFTL),1,WORK(IOFFVJ),1)

! BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
            SCL=0.5D0
            IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            IF (ILABS.EQ.IJABS) SCL=SCL*SQRTH
            BPTVJL=SCL*(TJVL+TLVJ)
! write element HP(ac,jl)
            IDX=ITGEU+NAS*(IJGEL-IISTA)
            DBL_MB(MW+IDX-1)=BPTVJL
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 2      CONTINUE
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

        IF(NW.EQ.0) GOTO 3

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
            NV=NVTOT_CHOSYM(MUL(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVL=IV-1+NASH(ISYV)*(IL-1)
            IOFFTJ=LCHOBUF+IOSYM(ISYT,ISYJ)+NV*ITJ
            IOFFVL=LCHOBUF+IOSYM(ISYV,ISYL)+NV*IVL
            TJVL=DDOT_(NV,WORK(IOFFTJ),1,WORK(IOFFVL),1)

            NV=NVTOT_CHOSYM(MUL(ISYT,ISYL))
            ITL=IT-1+NASH(ISYT)*(IL-1)
            IVJ=IV-1+NASH(ISYV)*(IJ-1)
            IOFFTL=LCHOBUF+IOSYM(ISYT,ISYL)+NV*ITL
            IOFFVJ=LCHOBUF+IOSYM(ISYV,ISYJ)+NV*IVJ
            TLVJ=DDOT_(NV,WORK(IOFFTL),1,WORK(IOFFVJ),1)

! BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
            SCL=0.5D0
            !IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            !IF (ILABS.EQ.IJABS) SCL=SCL*SQRTH
            BMTVJL=SCL*(TJVL-TLVJ)
! write element BM(tv,jl)
            IDX=ITGTU+NAS*(IJGTL-IISTA)
            DBL_MB(MW+IDX-1)=BMTVJL
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 3      CONTINUE
      END DO
************************************************************************

      CALL GETMEM('CHOBUF','FREE','REAL',LCHOBUF,NCHOBUF)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_F_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOSYM(8,8)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case F'
      END IF

************************************************************************
* Case F (8,9):
C FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
C FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
************************************************************************

      SQRTH=SQRT(0.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NCHOBUF,IOSYM)

      CALL GETMEM('CHOBUF','ALLO','REAL',LCHOBUF,NCHOBUF)

      CALL CHOVEC_READ(3,LCHOBUF)

      iCASE=8
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 8

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
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            ICV=IC-1+NSSH(ISYC)*(IV-1)
            IOFFAT=LCHOBUF+IOSYM(ISYA,ISYT)+NV*IAT
            IOFFCV=LCHOBUF+IOSYM(ISYC,ISYV)+NV*ICV
            ATCV=DDOT_(NV,WORK(IOFFAT),1,WORK(IOFFCV),1)

            NV=NVTOT_CHOSYM(MUL(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAV=IA-1+NSSH(ISYA)*(IV-1)
            ICT=IC-1+NSSH(ISYC)*(IT-1)
            IOFFAV=LCHOBUF+IOSYM(ISYA,ISYV)+NV*IAV
            IOFFCT=LCHOBUF+IOSYM(ISYC,ISYT)+NV*ICT
            AVCT=DDOT_(NV,WORK(IOFFAV),1,WORK(IOFFCT),1)

! FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
            SCL=0.5D0
            IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            IF (IAABS.EQ.ICABS) SCL=SCL*SQRTH
            FPTVAC=SCL*(ATCV+AVCT)
! write element FP(tv,ac)
            IDX=ITGEU+NAS*(IAGEB-IISTA)
            DBL_MB(MW+IDX-1)=FPTVAC
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 8      CONTINUE
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

        IF(NW.EQ.0) GOTO 9

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
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            ICV=IC-1+NSSH(ISYC)*(IV-1)
            IOFFAT=LCHOBUF+IOSYM(ISYA,ISYT)+NV*IAT
            IOFFCV=LCHOBUF+IOSYM(ISYC,ISYV)+NV*ICV
            ATCV=DDOT_(NV,WORK(IOFFAT),1,WORK(IOFFCV),1)

            NV=NVTOT_CHOSYM(MUL(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAV=IA-1+NSSH(ISYA)*(IV-1)
            ICT=IC-1+NSSH(ISYC)*(IT-1)
            IOFFAV=LCHOBUF+IOSYM(ISYA,ISYV)+NV*IAV
            IOFFCT=LCHOBUF+IOSYM(ISYC,ISYT)+NV*ICT
            AVCT=DDOT_(NV,WORK(IOFFAV),1,WORK(IOFFCT),1)

! FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
            SCL=0.5D0
            !IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            !IF (IAABS.EQ.ICABS) SCL=SCL*SQRTH
            FMTVAC=SCL*(AVCT-ATCV)
! write element FM(tv,ac)
            IDX=ITGTU+NAS*(IAGTB-IISTA)
            DBL_MB(MW+IDX-1)=FMTVAC
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 9      CONTINUE
      END DO
************************************************************************

      CALL GETMEM('CHOBUF','FREE','REAL',LCHOBUF,NCHOBUF)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_H_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOSYM(8,8)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif
      INTEGER, PARAMETER :: NOSYM = 1
      REAL*8, ALLOCATABLE :: AIBJ(:,:)

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case H'
      END IF

************************************************************************
* Case H:
C   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
C   WM(jl,ac)=((ajcl)-(alcj))*SQRT(3.0D0)
************************************************************************

      SQRT3=SQRT(3.0D0)
      SQRTH=SQRT(0.5D0)

      NV=NVTOT_CHOSYM(NOSYM)
      ALLOCATE(AIBJ(NSSHT,NSSHT))
      NBLOCK=NV*NSSHT

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NCHOBUF,IOSYM)

      CALL GETMEM('CHOBUF','ALLO','REAL',LCHOBUF,NCHOBUF)

      CALL CHOVEC_READ(4,LCHOBUF)

      iCASE=12

      NAS=NAGEB(NOSYM)
      NIS=NIGEJ(NOSYM)
      NW=NAS*NIS

      IF(NW.EQ.0) GOTO 12

      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
      NW=NAS*(IIEND-IISTA+1)

      DO IJGEL=IISTA,IIEND
        IJ=MIGEJ(1,IJGEL)
        IL=MIGEJ(2,IJGEL)
        LIJOFF=LCHOBUF+NBLOCK*(IJ-1)
        LILOFF=LCHOBUF+NBLOCK*(IL-1)
        ! precompute integral blocks
        CALL DGEMM_('T','N',NSSHT,NSSHT,NV,
     &              1.0D0,WORK(LIJOFF),NV,WORK(LILOFF),NV,
     &              0.0D0,AIBJ,NSSHT)
        DO IAGEB=IASTA,IAEND ! these are always all elements
          IA=MAGEB(1,IAGEB)
          IC=MAGEB(2,IAGEB)
          AJCL=AIBJ(IA,IC)
          ALCJ=AIBJ(IC,IA)
! HP(ac,jl)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
          SCL=1.0D0
          IF (IA.EQ.IC) SCL=SCL*SQRTH
          IF (IL.EQ.IJ) SCL=SCL*SQRTH
          HPACJL=SCL*(AJCL+ALCJ)
! write element HP(ac,jl)
          IDX=IAGEB+NAS*(IJGEL-IISTA)
          DBL_MB(MW+IDX-1)=HPACJL
        END DO
      END DO
************************************************************************
      CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
      CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,NOSYM,iVEC)
      CALL RHS_FREE (NAS,NIS,lg_W)
************************************************************************
 12   CONTINUE

      iCASE=13

      NAS=NAGTB(NOSYM)
      NIS=NIGTJ(NOSYM)
      NW=NAS*NIS

      IF(NW.EQ.0) GOTO 13

      CALL RHS_ALLO (NAS,NIS,lg_W)
      CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
      NW=NAS*(IIEND-IISTA+1)

      DO IJGTL=IISTA,IIEND
        IJ=MIGTJ(1,IJGTL)
        IL=MIGTJ(2,IJGTL)
        LIJOFF=LCHOBUF+NBLOCK*(IJ-1)
        LILOFF=LCHOBUF+NBLOCK*(IL-1)
        ! precompute integral blocks
        CALL DGEMM_('T','N',NSSHT,NSSHT,NV,
     &              1.0D0,WORK(LIJOFF),NV,WORK(LILOFF),NV,
     &              0.0D0,AIBJ,NSSHT)
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
          DBL_MB(MW+IDX-1)=HMACJL
        END DO
      END DO
************************************************************************

      CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
      CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,NOSYM,iVEC)
      CALL RHS_FREE (NAS,NIS,lg_W)
 13   CONTINUE
************************************************************************

      CALL GETMEM('CHOBUF','FREE','REAL',LCHOBUF,NCHOBUF)

      DEALLOCATE(AIBJ)

      END


*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_D_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOBRA1(8,8), IOKET1(8,8), IOBRA2(8,8), IOKET2(8,8)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif
      DIMENSION NFIMOES(8)

      IF (iPrGlb.GE.DEBUG) THEN
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

      CALL GETMEM('BRABUF1','ALLO','REAL',LBRABUF1,NBRABUF1)
      CALL GETMEM('KETBUF1','ALLO','REAL',LKETBUF1,NKETBUF1)

      CALL CHOVEC_READ(4,LBRABUF1)
      CALL CHOVEC_READ(2,LKETBUF1)

      CALL CHOVEC_SIZE(3,NBRABUF2,IOBRA2)
      CALL CHOVEC_SIZE(1,NKETBUF2,IOKET2)

      CALL GETMEM('BRABUF2','ALLO','REAL',LBRABUF2,NBRABUF2)
      CALL GETMEM('KETBUF2','ALLO','REAL',LKETBUF2,NKETBUF2)

      CALL CHOVEC_READ(3,LBRABUF2)
      CALL CHOVEC_READ(1,LKETBUF2)

      iCASE=5
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      ! set up FIMO access
      ACTINV=1.0D0/DBLE(MAX(1,NACTEL))
      IFIMOES=0
      DO ISYM=1,NSYM
        NFIMOES(ISYM)=IFIMOES
        IFIMOES=IFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2
      END DO

      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 8

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
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ))
            IOAJ=IA-1+NSSH(ISYA)*(IJ-1)
            IOTV=IT-1+NASH(ISYT)*(IV-1)
            IOFFAJ=LBRABUF1+IOBRA1(ISYA,ISYJ)+NV*IOAJ
            IOFFTV=LKETBUF1+IOKET1(ISYT,ISYV)+NV*IOTV
            AJTV=DDOT_(NV,WORK(IOFFAJ),1,WORK(IOFFTV),1)

! D1(tv,aj)=(aj,tv) + FIMO(a,j)*Kron(t,v)/NACTEL
! integrals only
            IDX=ITV+NAS*(IAJ-IISTA)
            DBL_MB(MW+IDX-1)=AJTV
          END DO
! now, dress with FIMO(a,j), only if T==V, so ISYT==ISYV, so if ISYM==1
          IF (ISYM.EQ.1) THEN
            IATOT=IA+NISH(ISYA)+NASH(ISYA)
            FAJ=WORK(LFIMO+NFIMOES(ISYA)+(IATOT*(IATOT-1))/2+IJ-1)
            ONEADD=FAJ*ACTINV
            DO IUABS=1,NASHT
              IUU=KTU(IUABS,IUABS)
              IDX=IUU+NAS*(IAJ-IISTA)
              DBL_MB(MW+IDX-1)=DBL_MB(MW+IDX-1)+ONEADD
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
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYV))
            IOAV=IA-1+NSSH(ISYA)*(IV-1)
            IOTJ=IT-1+NASH(ISYT)*(IJ-1)
            IOFFAV=LBRABUF2+IOBRA2(ISYA,ISYV)+NV*IOAV
            IOFFTJ=LKETBUF2+IOKET2(ISYT,ISYJ)+NV*IOTJ
            AVTJ=DDOT_(NV,WORK(IOFFAV),1,WORK(IOFFTJ),1)

! D2(tv,aj)=(av,tj) + FIMO(a,j)*Kron(t,v)/NACTEL
            IDX=ITV+NAS*(IAJ-IISTA)
            DBL_MB(MW+IDX-1)=AVTJ
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 8      CONTINUE

      END DO
************************************************************************

      CALL GETMEM('BRABUF1','FREE','REAL',LBRABUF1,NBRABUF1)
      CALL GETMEM('KETBUF1','FREE','REAL',LKETBUF1,NKETBUF1)

      CALL GETMEM('BRABUF2','FREE','REAL',LBRABUF2,NBRABUF2)
      CALL GETMEM('KETBUF2','FREE','REAL',LKETBUF2,NKETBUF2)

************************************************************************

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_E_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOBRA(8,8), IOKET(8,8)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif

      IF (iPrGlb.GE.DEBUG) THEN
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

      SQRTH=SQRT(0.5D0)
      SQRTA=SQRT(1.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NBRABUF,IOBRA)
      CALL CHOVEC_SIZE(1,NKETBUF,IOKET)

      CALL GETMEM('BRABUF','ALLO','REAL',LBRABUF,NBRABUF)
      CALL GETMEM('KETBUF','ALLO','REAL',LKETBUF,NKETBUF)

      CALL CHOVEC_READ(4,LBRABUF)
      CALL CHOVEC_READ(1,LKETBUF)

      iCASE=6
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 6

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYA=1,NSYM
          ISYJL=MUL(ISYA,ISYM)
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
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ))
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IVL=IV-1+NASH(ISYV)*(IL-1)
              IOFFAJ=LBRABUF+IOBRA(ISYA,ISYJ)+NV*IAJ
              IOFFVL=LKETBUF+IOKET(ISYV,ISYL)+NV*IVL
              AJVL=DDOT_(NV,WORK(IOFFAJ),1,WORK(IOFFVL),1)

              NV=NVTOT_CHOSYM(MUL(ISYA,ISYL))
              IAL=IA-1+NSSH(ISYA)*(IL-1)
              IVJ=IV-1+NASH(ISYV)*(IJ-1)
              IOFFAL=LBRABUF+IOBRA(ISYA,ISYL)+NV*IAL
              IOFFVJ=LKETBUF+IOKET(ISYV,ISYJ)+NV*IVJ
              ALVJ=DDOT_(NV,WORK(IOFFAL),1,WORK(IOFFVJ),1)

! EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
              IF (ILABS.EQ.IJABS) THEN
                SCL=0.5D0
              ELSE
                SCL=SQRTH
              END IF
              EP=SCL*(AJVL+ALVJ)
! write element EP
              IDX=IV+NAS*(IAJGEL+IOFF-IISTA)
              DBL_MB(MW+IDX-1)=EP
            END DO
          END DO

          IOFF=IOFF+NA*NJL
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 6      CONTINUE
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

        IF(NW.EQ.0) GOTO 7

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYA=1,NSYM
          ISYJL=MUL(ISYA,ISYM)
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
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ))
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IVL=IV-1+NASH(ISYV)*(IL-1)
              IOFFAJ=LBRABUF+IOBRA(ISYA,ISYJ)+NV*IAJ
              IOFFVL=LKETBUF+IOKET(ISYV,ISYL)+NV*IVL
              AJVL=DDOT_(NV,WORK(IOFFAJ),1,WORK(IOFFVL),1)

              NV=NVTOT_CHOSYM(MUL(ISYA,ISYL))
              IAL=IA-1+NSSH(ISYA)*(IL-1)
              IVJ=IV-1+NASH(ISYV)*(IJ-1)
              IOFFAL=LBRABUF+IOBRA(ISYA,ISYL)+NV*IAL
              IOFFVJ=LKETBUF+IOKET(ISYV,ISYJ)+NV*IVJ
              ALVJ=DDOT_(NV,WORK(IOFFAL),1,WORK(IOFFVJ),1)

! EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
              EM=SQRTA*(AJVL-ALVJ)
! write element EM
              IDX=IV+NAS*(IAJGTL+IOFF-IISTA)
              DBL_MB(MW+IDX-1)=EM
            END DO
          END DO

          IOFF=IOFF+NA*NJL
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 7      CONTINUE
      END DO
************************************************************************

      CALL GETMEM('BRABUF','FREE','REAL',LBRABUF,NBRABUF)
      CALL GETMEM('KETBUF','FREE','REAL',LKETBUF,NKETBUF)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_G_NOSYM(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION IOBRA(8,8), IOKET(8,8)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#else
      REAL*8 DBL_MB(1:IWORKLEN)
      EQUIVALENCE (DBL_MB,WORK)
#endif

      IF (iPrGlb.GE.DEBUG) THEN
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

      SQRTH=SQRT(0.5D0)
      SQRTA=SQRT(1.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NBRABUF,IOBRA)
      CALL CHOVEC_SIZE(4,NKETBUF,IOKET)

      CALL GETMEM('BRABUF','ALLO','REAL',LBRABUF,NBRABUF)
      CALL GETMEM('KETBUF','ALLO','REAL',LKETBUF,NKETBUF)

      CALL CHOVEC_READ(3,LBRABUF)
      CALL CHOVEC_READ(4,LKETBUF)

      iCASE=10
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 10

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYJ=1,NSYM
          ISYAC=MUL(ISYJ,ISYM)
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
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYV))
              IAV=IA-1+NSSH(ISYA)*(IV-1)
              ICJ=IC-1+NSSH(ISYC)*(IJ-1)
              IOFFAV=LBRABUF+IOBRA(ISYA,ISYV)+NV*IAV
              IOFFCJ=LKETBUF+IOKET(ISYC,ISYJ)+NV*ICJ
              AVCJ=DDOT_(NV,WORK(IOFFAV),1,WORK(IOFFCJ),1)

              NV=NVTOT_CHOSYM(MUL(ISYC,ISYV))
              ICV=IC-1+NSSH(ISYC)*(IV-1)
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IOFFCV=LBRABUF+IOBRA(ISYC,ISYV)+NV*ICV
              IOFFAJ=LKETBUF+IOKET(ISYA,ISYJ)+NV*IAJ
              CVAJ=DDOT_(NV,WORK(IOFFCV),1,WORK(IOFFAJ),1)

C GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
              IF (IAABS.EQ.ICABS) THEN
                SCL=0.5D0
              ELSE
                SCL=SQRTH
              END IF
              GP=SCL*(AVCJ+CVAJ)
! write element EP
              IDX=IV+NAS*(IJAGEC+IOFF-IISTA)
              DBL_MB(MW+IDX-1)=GP
            END DO
          END DO

          IOFF=IOFF+NJ*NAC
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 10     CONTINUE
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

        IF(NW.EQ.0) GOTO 11

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYJ=1,NSYM
          ISYAC=MUL(ISYJ,ISYM)
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
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYV))
              IAV=IA-1+NSSH(ISYA)*(IV-1)
              ICJ=IC-1+NSSH(ISYC)*(IJ-1)
              IOFFAV=LBRABUF+IOBRA(ISYA,ISYV)+NV*IAV
              IOFFCJ=LKETBUF+IOKET(ISYC,ISYJ)+NV*ICJ
              AVCJ=DDOT_(NV,WORK(IOFFAV),1,WORK(IOFFCJ),1)

              NV=NVTOT_CHOSYM(MUL(ISYC,ISYV))
              ICV=IC-1+NSSH(ISYC)*(IV-1)
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IOFFCV=LBRABUF+IOBRA(ISYC,ISYV)+NV*ICV
              IOFFAJ=LKETBUF+IOKET(ISYA,ISYJ)+NV*IAJ
              CVAJ=DDOT_(NV,WORK(IOFFCV),1,WORK(IOFFAJ),1)

C GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
              GM=SQRTA*(AVCJ-CVAJ)
! write element GM
              IDX=IV+NAS*(IJAGTC+IOFF-IISTA)
              DBL_MB(MW+IDX-1)=GM
            END DO
          END DO

          IOFF=IOFF+NJ*NAC
        END DO
************************************************************************

        CALL RHS_Release_Update (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (NAS,NIS,lg_W)
 11     CONTINUE
      END DO
************************************************************************

      CALL GETMEM('BRABUF','FREE','REAL',LBRABUF,NBRABUF)
      CALL GETMEM('KETBUF','FREE','REAL',LKETBUF,NKETBUF)

      RETURN
      END
