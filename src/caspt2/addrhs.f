************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE ADDRHSA(IVEC,JSYM,ISYJ,ISYX,NT,NJ,NV,NX,TJVX,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
      DIMENSION TJVX(NT,NJ,NV,NX)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
*      Logical Incore
*                                                                      *
************************************************************************
*                                                                      *
*
      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYJ
      IF(NINDEP(ISYM,1).EQ.0) GOTO 900
      NAS=NTUV(ISYM)
      NIS=NISH(ISYM)
      NWA=NAS*NIS
      IF(NWA.EQ.0) GOTO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option.
*
*     Incore=nBuff.ge.NWA+NT*NJ*NV*NX
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSA'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
      CALL DGEMM_('N','T',NT*NJ,NV*NX,NCHO,
     &            1.0D0,Cho_Bra,NT*NJ,
     &                  Cho_Ket,NV*NX,
     &            0.0D0,TJVX,NT*NJ)
*
C Compute W(tvx,j)=(tj,vx) + FIMO(t,j)*delta(v,x)/NACTEL
      ICASE=1
*     LWA=1+NT*NJ*NV*NX
      LDA=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_A)
      CALL RHS_READ (NAS,NIS,lg_A,iCASE,iSYM,iVEC)

      IBUF=0
      DO IT=1,NT
       ITABS=IT+NAES(ISYT)
       DO IJ=1,NJ
        DO IV=1,NV
         IVABS=IV+NAES(ISYV)
         DO IX=1,NX
          IXABS=IX+NAES(ISYX)
          IW1=KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
          IW2=IJ
          IW=IW1+NAS*(IW2-1)
*SVC      Buff(LWA-1+IW)=Buff(LWA-1+IW)+TJVX(IT,IJ,IV,IX)
          IBUF=IBUF+1
          idxBuf(IBUF)=IW
          Buff(IBUF)=TJVX(IT,IJ,IV,IX)
          IF (IBUF.EQ.NBUFF) THEN
            CALL RHS_SCATTER(LDA,lg_A,Buff,idxBuf,IBUF)
            IBUF=0
          END IF
         END DO
        END DO
       END DO
      END DO
      IF (IBUF.NE.0) THEN
        CALL RHS_SCATTER(LDA,lg_A,Buff,idxBuf,IBUF)
      END IF
C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_A,iCASE,iSYM,iVEC)
      CALL RHS_FREE(NAS,NIS,lg_A)
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSB(IVEC,JSYM,ISYJ,ISYL,NT,NJ,NV,NL,TJVL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
Case B:
* t>v j>l WP(tv,jl)=add ((tj,vl))*(1/2)
* t>v j=l WP(tv,jl)=add ((tj,vl))*(1/2)*(SQRT(2))
* t>v j<l WP(tv,lj)=add ((tj,vl))*(1/2)

* t=v j>l WP(tv,jl)=add ((tj,vl))*(1/4)
* t=v j=l WP(tv,jl)=add ((tj,vl))*(1/4)*(SQRT(2))
* t=v j<l WP(tv,lj)=add ((tj,vl))*(1/4)

      DIMENSION TJVL(NT,NJ,NV,NL)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
*                                                                      *
************************************************************************
*                                                                      *
*
      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      IF(ISYT.LT.ISYV) GOTO 900
      SQ2=SQRT(2.0D0)
      ISYM=MUL(ISYJ,ISYL)
*
      IF(NINDEP(ISYM,2).GT.0) THEN
* The plus combination:
       ICASE=2
       NASP=NTGEU(ISYM)
       NISP=NIGEJ(ISYM)
       NWBP=NASP*NISP
      ELSE
       NWBP=0
      ENDIF
      IF(NINDEP(ISYM,3).GT.0) THEN
* The minus combination:
       ICASE=3
       NASM=NTGTU(ISYM)
       NISM=NIGTJ(ISYM)
       NWBM=NASM*NISM
      ELSE
       NWBM=0
      ENDIF
      If (Max(NWBP,NWBM).le.0) GO TO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option.
*
*     Incore=nBuff.ge.Max(NWBP,NWBM)+NT*NJ*NV*NL
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSB'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
      CALL DGEMM_('N','T',NT*NJ,NV*NL,NCHO,
     &            1.0D0,Cho_Bra,NT*NJ,
     &                  Cho_Ket,NV*NL,
     &            0.0D0,TJVL,NT*NJ)
      If (NWBP.le.0) GO TO 800
*
      IF(NINDEP(ISYM,2).GT.0) THEN
* The plus combination:
       ICASE=2
       NASP=NTGEU(ISYM)
       NISP=NIGEJ(ISYM)
       NWBP=NASP*NISP
*      LWBP=1+NT*NJ*NV*NL
       LDBP=NASP
C Read WP:
       CALL RHS_ALLO (NASP,NISP,lg_BP)
       CALL RHS_READ (NASP,NISP,lg_BP,iCASE,iSYM,iVEC)

       IBUF=0
        DO IT=1,NT
        ITABS=IT+NAES(ISYT)
        IVMAX=NV
        IF(ISYV.EQ.ISYT) IVMAX=IT
        DO IV=1,IVMAX
         IVABS=IV+NAES(ISYV)
         SCL1=0.5D0
         IW1=KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
         IF(ITABS.EQ.IVABS) SCL1=0.25D0
         DO IJ=1,NJ
          IJABS=IJ+NIES(ISYJ)
          DO IL=1,NL
           ILABS=IL+NIES(ISYL)
           SCL=SCL1
           IF(IJABS.GE.ILABS) THEN
            IW2=KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
            IF(IJABS.EQ.ILABS) SCL=SQ2*SCL1
           ELSE
            IW2=KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
           END IF
           IW=IW1+NASP*(IW2-1)
*          Buff(LWBP-1+IW)=Buff(LWBP-1+IW)+SCL*TJVL(IT,IJ,IV,IL)
           IBUF=IBUF+1
           idxBuf(IBUF)=IW
           Buff(IBUF)=SCL*TJVL(IT,IJ,IV,IL)
           IF (IBUF.EQ.NBUFF) THEN
             CALL RHS_SCATTER(LDBP,lg_BP,Buff,idxBuf,IBUF)
             IBUF=0
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDBP,lg_BP,Buff,idxBuf,IBUF)
       END IF
C Put WBP on disk:
       CALL RHS_SAVE (NASP,NISP,lg_BP,iCASE,iSYM,iVEC)
       CALL RHS_FREE(NASP,NISP,lg_BP)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
 800  CONTINUE
      IF(NINDEP(ISYM,3).GT.0) THEN
* The minus combination:
       ICASE=3
       NASM=NTGTU(ISYM)
       NISM=NIGTJ(ISYM)
       NWBM=NASM*NISM
*      LWBM=1+NT*NJ*NV*NL
       LDBM=NASM
C Read WM:
       CALL RHS_ALLO (NASM,NISM,lg_BM)
       CALL RHS_READ (NASM,NISM,lg_BM,iCASE,iSYM,iVEC)

       IBUF=0
       DO IT=1,NT
        ITABS=IT+NAES(ISYT)
        IVMAX=NV
        IF(ISYV.EQ.ISYT) IVMAX=IT-1
        DO IV=1,IVMAX
         IVABS=IV+NAES(ISYV)
         IW1=KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
         DO IJ=1,NJ
          IJABS=IJ+NIES(ISYJ)
          DO IL=1,NL
           ILABS=IL+NIES(ISYL)
           IF(IJABS.GT.ILABS) THEN
            IW2=KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
            IW=IW1+NASM*(IW2-1)
*           Buff(LWBM-1+IW)=Buff(LWBM-1+IW)+0.5D0*TJVL(IT,IJ,IV,IL)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=0.5D0*TJVL(IT,IJ,IV,IL)
           ELSE IF (IJABS.LT.ILABS) THEN
            IW2=KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
            IW=IW1+NASM*(IW2-1)
*           Buff(LWBM-1+IW)=Buff(LWBM-1+IW)-0.5D0*TJVL(IT,IJ,IV,IL)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=-0.5D0*TJVL(IT,IJ,IV,IL)
           END IF
           IF (IBUF.EQ.NBUFF) THEN
             CALL RHS_SCATTER(LDBM,lg_BM,Buff,idxBuf,IBUF)
             IBUF=0
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDBM,lg_BM,Buff,idxBuf,IBUF)
       END IF
C Put WBM on disk:
       CALL RHS_SAVE (NASM,NISM,lg_BM,iCASE,iSYM,iVEC)
       CALL RHS_FREE (NASM,NISM,lg_BM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSC(IVEC,JSYM,ISYU,ISYX,NA,NU,NV,NX,AUVX,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION AUVX(NA,NU,NV,NX)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NX,NCHO)
*      Logical Incore
Case C:
C   Allocate W. Put in W(uvx,a)=(au,vx) +
C             (FIMO(a,t)-sum(y)(ay,yt))*delta(u,v)/NACTEL.
*                                                                      *
************************************************************************
*                                                                      *
*
      ISYA=MUL(JSYM,ISYU)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYA
      IF(NINDEP(ISYM,4).EQ.0) GOTO 900
      NAS=NTUV(ISYM)
      NIS=NSSH(ISYM)
      NWC=NAS*NIS
      IF(NWC.EQ.0) GOTO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option.
*
*     Incore=nBuff.ge.NWC+NA*NU*NV*NX
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSC'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
      CALL DGEMM_('N','T',NA*NU,NV*NX,NCHO,
     &            1.0D0,Cho_Bra,NA*NU,
     &                  Cho_Ket,NV*NX,
     &            0.0D0,AUVX,NA*NU)
*
      ICASE=4
*     LWC=1+NA*NU*NV*NX
      LDC=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_C)
      CALL RHS_READ (NAS,NIS,lg_C,iCASE,iSYM,iVEC)

      IBUF=0
      DO IA=1,NA
       DO IU=1,NU
        IUABS=IU+NAES(ISYU)
        DO IV=1,NV
         IVABS=IV+NAES(ISYV)
         DO IX=1,NX
          IXABS=IX+NAES(ISYX)
          IW1=KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
          IW2=IA
          IW=IW1+NAS*(IW2-1)
*         Buff(LWC-1+IW)=Buff(LWC-1+IW)+AUVX(IA,IU,IV,IX)
          IBUF=IBUF+1
          idxBuf(IBUF)=IW
          Buff(IBUF)=AUVX(IA,IU,IV,IX)
          IF (IBUF.EQ.NBUFF) THEN
            CALL RHS_SCATTER(LDC,lg_C,Buff,idxBuf,IBUF)
            IBUF=0
          END IF
         END DO
        END DO
       END DO
      END DO
      IF (IBUF.NE.0) THEN
        CALL RHS_SCATTER(LDC,lg_C,Buff,idxBuf,IBUF)
      END IF
C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_C,iCASE,iSYM,iVEC)
      CALL RHS_FREE (NAS,NIS,lg_C)

*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSD1(IVEC,JSYM,ISYJ,ISYX,NA,NJ,NV,NX,AJVX,
     &                    nBuff,Buff,idxBuf,
     &                    Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION AJVX(NV,NX,*)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV,NX,NCHO)
*      Logical Incore
      DIMENSION IOFFD(8,8)
Case D:
C Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
C Compute W2(vu,al)=(au,vl)
*                                                                      *
************************************************************************
*                                                                      *
*
      DO ISW=1,NSYM
       IO=0
       DO ISA=1,NSYM
        IOFFD(ISA,ISW)=IO
        ISI=MUL(ISA,ISW)
        IO=IO+NSSH(ISA)*NISH(ISI)
       END DO
      END DO

      ISYA=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=JSYM
      IF(NINDEP(ISYM,5).EQ.0) GOTO 900
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) GOTO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option?
*
*     Incore=nBuff.ge.NWD+NA*NJ*NV*NX
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSD1'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
C     CALL DGEMM_('N','T',NA*NJ,NV*NX,NCHO,
C    &            1.0D0,Cho_Bra,NA*NJ,
C    &                  Cho_Ket,NV*NX,
C    &            0.0D0,AJVX,NA*NJ)
*

C Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
      ICASE=5
*     LWD=1+NA*NJ*NV*NX
      LDD=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_D)
      CALL RHS_READ (NAS,NIS,lg_D,iCASE,iSYM,iVEC)

      NBXSZA=NSECBX
      NBXSZJ=NINABX

      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        NASZ=IAEND-IASTA+1
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          NJSZ=IJEND-IJSTA+1

          IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
          CALL DGEMM_('N','T',NV*NX,NASZ*NJSZ,NCHO,
     &         1.0D0,Cho_Ket,NV*NX,
     &         Cho_Bra(IAJSTA,1),NA*NJ,
     &         0.0D0,AJVX,NV*NX)

      IAJ=0
      IBUF=0
      DO IJ=IJSTA,IJEND
       DO IA=IASTA,IAEND
       IAJ=IAJ+1

        DO IX=1,NX
         IXABS=IX+NAES(ISYX)
          DO IV=1,NV
          IVABS=IV+NAES(ISYV)

          IW1=KTU(IVABS,IXABS)-NTUES(ISYM)
          IW2=IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
          IW=IW1+NAS*(IW2-1)
          IBUF=IBUF+1
*         WD=Buff(LWD-1+IW)+AJVX(IV,IX,IAJ)
*         Buff(LWD-1+IW)=WD
          idxBuf(IBUF)=IW
          Buff(IBUF)=AJVX(IV,IX,IAJ)
          IF (IBUF.EQ.NBUFF) THEN
            CALL RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
            IBUF=0
          END IF
         END DO
        END DO
       END DO
      END DO
      IF (IBUF.NE.0) THEN
        CALL RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
      END IF

        ENDDO
      ENDDO

C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_D,iCASE,iSYM,iVEC)
      CALL RHS_FREE (NAS,NIS,lg_D)
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSD2(IVEC,JSYM,ISYU,ISYL,NA,NU,NV,NL,AUVL,
     &                    nBuff,Buff,idxBuf,
     &                    Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION AUVL(NA,NU,NV,NL)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
      DIMENSION IOFFD(8,8)
Case D:
C Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
C Compute W2(vu,al)=(au,vl)
*                                                                      *
************************************************************************
*                                                                      *
*
      DO ISYW=1,NSYM
       IO=0
       DO ISYA=1,NSYM
        IOFFD(ISYA,ISYW)=IO
        ISYI=MUL(ISYA,ISYW)
        IO=IO+NSSH(ISYA)*NISH(ISYI)
       END DO
      END DO

      ISYA=MUL(JSYM,ISYU)
      ISYV=MUL(JSYM,ISYL)
      ISYM=MUL(ISYU,ISYV)
      IF(NINDEP(ISYM,5).EQ.0) GOTO 900
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) GOTO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option?
*
*     Incore=nBuff.ge.NWD+NA*NU*NV*NL
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSD2'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
*
      CALL DGEMM_('N','T',NA*NU,NV*NL,NCHO,
     &            1.0D0,Cho_Bra,NA*NU,
     &                  Cho_Ket,NV*NL,
     &            0.0D0,AUVL,NA*NU)
*
C Compute W2(vu,al)=(au,vl)
      ICASE=5
*     LWD=1+NA*NU*NV*NL
      LDD=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_D)
      CALL RHS_READ (NAS,NIS,lg_D,iCASE,iSYM,iVEC)

      IBUF=0
      DO IA=1,NA
       DO IU=1,NU
        IUABS=IU+NAES(ISYU)
        DO IV=1,NV
         IVABS=IV+NAES(ISYV)
         DO IL=1,NL
          IW1=NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
          IW2=IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
          IW=IW1+NAS*(IW2-1)
*         Buff(LWD-1+IW)=Buff(LWD-1+IW)+AUVL(IA,IU,IV,IL)
          IBUF=IBUF+1
          idxBuf(IBUF)=IW
          Buff(IBUF)=AUVL(IA,IU,IV,IL)
          IF (IBUF.EQ.NBUFF) THEN
            CALL RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
            IBUF=0
          END IF
         END DO
        END DO
       END DO
      END DO
      IF (IBUF.NE.0) THEN
        CALL RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
      END IF
*
C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_D,iCASE,iSYM,iVEC)
      CALL RHS_FREE (NAS,NIS,lg_D)
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSE(IVEC,JSYM,ISYJ,ISYL,NA,NJ,NV,NL,AJVL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

      DIMENSION AJVL(NV,NL,*)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
      DIMENSION IOFF1(8),IOFF2(8)
Case E:
*                                                                      *
************************************************************************
*                                                                      *
*
      SQ32=SQRT(1.5D0)
      ISYA=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      ISYM=ISYV
      ISYJL=MUL(ISYJ,ISYL)

C Set up offset table:
      IO1=0
      IO2=0
      DO ISA=1,NSYM
        IOFF1(ISA)=IO1
        IOFF2(ISA)=IO2
        ISIJ=MUL(ISA,ISYM)
        IO1=IO1+NSSH(ISA)*NIGEJ(ISIJ)
        IO2=IO2+NSSH(ISA)*NIGTJ(ISIJ)
      END DO

      NAS=NASH(ISYM)
      NISP=NISUP(ISYM,6)
      NISM=NISUP(ISYM,7)
      NWP=NAS*NISP
      NWM=NAS*NISM
      NW=NWP+NWM
      If (NW.eq.0) GO TO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option?
*
*     Incore=nBuff.ge.Max(NWP,NWM)+NA*NJ*NV*NL
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSE'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
C     CALL DGEMM_('N','T',NA*NJ,NV*NL,NCHO,
C    &            1.0D0,Cho_Bra,NA*NJ,
C    &                  Cho_Ket,NV*NL,
C    &            0.0D0,AJVL,NA*NJ)
*
*     LWE=1+NA*NJ*NV*NL
*     LWEP=LWE
*     LWEM=LWE
      LDEP=NAS
      LDEM=NAS
* The plus combination:
      IF (NWP.GT.0) THEN
       ICASE=6
C Read WP:
      CALL RHS_ALLO (NAS,NISP,lg_EP)
      CALL RHS_READ (NAS,NISP,lg_EP,iCASE,iSYM,iVEC)

       NBXSZA=NSECBX
       NBXSZJ=NINABX

       DO IASTA=1,NA,NBXSZA
         IAEND=MIN(IASTA-1+NBXSZA,NA)
         NASZ=IAEND-IASTA+1
         DO IJSTA=1,NJ,NBXSZJ
           IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
           NJSZ=IJEND-IJSTA+1

           IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
           CALL DGEMM_('N','T',NV*NL,NASZ*NJSZ,NCHO,
     &          1.0D0,Cho_Ket,NV*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          0.0D0,AJVL,NV*NL)

       IAJ=0
       IBUF=0
       DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        DO IA=IASTA,IAEND
        IAJ=IAJ+1

         DO IV=1,NV
          DO IL=1,NL
           ILABS=IL+NIES(ISYL)
           SCL=SQRT(0.5D0)
           IF(IJABS.GE.ILABS) THEN
            JGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
            IF(IJABS.EQ.ILABS) SCL=1.0D0
           ELSE
            JGEL=KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
           END IF
           IW1=IV
           IW2=IA+NA*(JGEL-1)+IOFF1(ISYA)
           IW=IW1+NAS*(IW2-1)
C   WP(v,a,jl)=  ((ajvl)+(alvj))/SQRT(2+2*Kron(jl))
*          Buff(LWEP-1+IWEP)=Buff(LWEP-1+IWEP)+SCL*AJVL(IV,IL,IAJ)
           IBUF=IBUF+1
           idxBuf(IBUF)=IW
           Buff(IBUF)=SCL*AJVL(IV,IL,IAJ)
           IF (IBUF.EQ.NBUFF) THEN
            CALL RHS_SCATTER(LDEP,lg_EP,Buff,idxBuf,IBUF)
            IBUF=0
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDEP,lg_EP,Buff,idxBuf,IBUF)
       END IF

         ENDDO
       ENDDO

      CALL RHS_SAVE (NAS,NISP,lg_EP,iCASE,iSYM,iVEC)
      CALL RHS_FREE (NAS,NISP,lg_EP)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
* The minus combination:
      IF (NWM.GT.0) THEN
       ICASE=7
C Read WM:
      CALL RHS_ALLO (NAS,NISM,lg_EM)
      CALL RHS_READ (NAS,NISM,lg_EM,iCASE,iSYM,iVEC)

       NBXSZA=NSECBX
       NBXSZJ=NINABX

       DO IASTA=1,NA,NBXSZA
         IAEND=MIN(IASTA-1+NBXSZA,NA)
         NASZ=IAEND-IASTA+1
         DO IJSTA=1,NJ,NBXSZJ
           IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
           NJSZ=IJEND-IJSTA+1

           IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
           CALL DGEMM_('N','T',NV*NL,NASZ*NJSZ,NCHO,
     &          1.0D0,Cho_Ket,NV*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          0.0D0,AJVL,NV*NL)

       IAJ=0
       IBUF=0
       DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        DO IA=IASTA,IAEND
        IAJ=IAJ+1

         DO IV=1,NV
          DO IL=1,NL
           ILABS=IL+NIES(ISYL)
           IF(IJABS.NE.ILABS) THEN
            IF(IJABS.GT.ILABS) THEN
             SCL=SQ32
             JGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
            ELSE
             SCL=-SQ32
             JGTL=KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
            END IF
            IW1=IV
            IW2=IA+NA*(JGTL-1)+IOFF2(ISYA)
            IW=IW1+NAS*(IW2-1)
*           Buff(LWEM-1+IWEM)=Buff(LWEM-1+IWEM)+SCL*AJVL(IV,IL,IAJ)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=SCL*AJVL(IV,IL,IAJ)
            IF (IBUF.EQ.NBUFF) THEN
             CALL RHS_SCATTER(LDEM,lg_EM,Buff,idxBuf,IBUF)
             IBUF=0
            END IF
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDEM,lg_EM,Buff,idxBuf,IBUF)
       END IF

         ENDDO
       ENDDO

      CALL RHS_SAVE (NAS,NISM,lg_EM,iCASE,iSYM,iVEC)
      CALL RHS_FREE (NAS,NISM,lg_EM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSF(IVEC,JSYM,ISYU,ISYX,NA,NU,NC,NX,AUCX,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION AUCX(NA,NU,NC,NX)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO)
*      Logical Incore
Case F:
C   WP(ux,ac)=((aucx)+(axcu))*(1-Kron(x,u)/2) /2
C With new normalisation, replace /2 with /(2*SQRT(1+Kron(ac))
C   WM(ux,ac)= -((aucx)-(axcu))/2
*                                                                      *
************************************************************************
*                                                                      *
*
      IF(ISYU.LT.ISYX) GOTO 900

      ISYA=MUL(JSYM,ISYU)
      ISYC=MUL(JSYM,ISYX)
      ISYM=MUL(ISYU,ISYX)
*
      IF(NINDEP(ISYM,8).GT.0) THEN
* The plus combination:
       NASP=NTGEU(ISYM)
       NISP=NAGEB(ISYM)
       NWFP=NASP*NISP
      ELSE
       NWFP=0
      ENDIF
      IF(NINDEP(ISYM,9).GT.0) THEN
       ICASE=9
* The minus combination:
       NASM=NTGTU(ISYM)
       NISM=NAGTB(ISYM)
       NWFM=NASM*NISM
      ELSE
       NWFM=0
      ENDIF
      If (NWFP+NWFM.le.0) GO TO 900
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option?
*
*     Incore=nBuff.ge.Max(NWFP,NWFM)+NA*NU*NC*NX
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSF'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
      CALL DGEMM_('N','T',NA*NU,NC*NX,NCHO,
     &            1.0D0,Cho_Bra,NA*NU,
     &                  Cho_Ket,NC*NX,
     &            0.0D0,AUCX,NA*NU)
*
      IF (NWFP.le.0) GO TO 800
      IF(NINDEP(ISYM,8).GT.0) THEN
* The plus combination:
       ICASE=8
       NASP=NTGEU(ISYM)
       NISP=NAGEB(ISYM)
       NWFP=NASP*NISP
*      LWFP=1+NA*NU*NC*NX
       LDFP=NASP
C Read WP:
       CALL RHS_ALLO (NASP,NISP,lg_FP)
       CALL RHS_READ (NASP,NISP,lg_FP,iCASE,iSYM,iVEC)

       IBUF=0
       DO IU=1,NU
        IUABS=IU+NAES(ISYU)
        IXMAX=NX
        IF(ISYU.EQ.ISYX) IXMAX=IU
        DO IX=1,IXMAX
         IXABS=IX+NAES(ISYX)
         SCL1=0.5D0
         IF(IUABS.EQ.IXABS) SCL1=0.25D0
         IW1=KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
         DO IA=1,NA
          IAABS=IA+NSES(ISYA)
          DO IC=1,NC
           ICABS=IC+NSES(ISYC)
           SCL=SCL1
           IF(IAABS.GE.ICABS) THEN
            IW2=KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
            IF(IAABS.EQ.ICABS) SCL=SQRT(2.0D0)*SCL1
           ELSE
            IW2=KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
           END IF
           IW=IW1+NASP*(IW2-1)
*          WFP=Buff(LWFP-1+IW)+SCL*AUCX(IA,IU,IC,IX)
*          Buff(LWFP-1+IW)=WFP
           IBUF=IBUF+1
           idxBuf(IBUF)=IW
           Buff(IBUF)=SCL*AUCX(IA,IU,IC,IX)
           IF (IBUF.EQ.NBUFF) THEN
             CALL RHS_SCATTER(LDFP,lg_FP,Buff,idxBuf,IBUF)
             IBUF=0
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDFP,lg_FP,Buff,idxBuf,IBUF)
       END IF
       CALL RHS_SAVE (NASP,NISP,lg_FP,iCASE,iSYM,iVEC)
       CALL RHS_FREE (NASP,NISP,lg_FP)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
 800  CONTINUE
      IF (NWFM.le.0) GO TO 900
      IF(NINDEP(ISYM,9).GT.0) THEN
       ICASE=9
* The minus combination:
       NASM=NTGTU(ISYM)
       NISM=NAGTB(ISYM)
       NWFM=NASM*NISM
*      LWFM=1+NA*NU*NC*NX
       LDFM=NASM
C Read WM:
       CALL RHS_ALLO (NASM,NISM,lg_FM)
       CALL RHS_READ (NASM,NISM,lg_FM,iCASE,iSYM,iVEC)

       IBUF=0
       DO IU=1,NU
        IUABS=IU+NAES(ISYU)
        IXMAX=NX
        IF(ISYU.EQ.ISYX) IXMAX=IU-1
        DO IX=1,IXMAX
         IXABS=IX+NAES(ISYX)
         IW1=KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
         DO IA=1,NA
          IAABS=IA+NSES(ISYA)
          DO IC=1,NC
           ICABS=IC+NSES(ISYC)
           IF(IAABS.GT.ICABS) THEN
            IW2=KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
            IW=IW1+NASM*(IW2-1)
*           Buff(LWFM-1+IW)=Buff(LWFM-1+IW)-0.5D0*AUCX(IA,IU,IC,IX)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=-0.5D0*AUCX(IA,IU,IC,IX)
           ELSE IF(IAABS.LT.ICABS) THEN
            IW2=KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
            IW=IW1+NASM*(IW2-1)
*           Buff(LWFM-1+IW)=Buff(LWFM-1+IW)+0.5D0*AUCX(IA,IU,IC,IX)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=0.5D0*AUCX(IA,IU,IC,IX)
           END IF
           IF (IBUF.EQ.NBUFF) THEN
             CALL RHS_SCATTER(LDFM,lg_FM,Buff,idxBuf,IBUF)
             IBUF=0
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDFM,lg_FM,Buff,idxBuf,IBUF)
       END IF
C Put WFM on disk:
       CALL RHS_SAVE (NASM,NISM,lg_FM,iCASE,iSYM,iVEC)
       CALL RHS_FREE (NASM,NISM,lg_FM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSG(IVEC,JSYM,ISYU,ISYL,NA,NU,NC,NL,AUCL,NAUCL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

      DIMENSION AUCL(NA,NU,*)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC*NL,NCHO)
*      Logical Incore
      DIMENSION IOFF1(8),IOFF2(8)
Case G:
C   WP(u,l,ac)=  ((aucl)+cual))/SQRT(2+2*Kron(ab))
C   WM(u,l,ac)=  ((aucl)-cual))*SQRT(1.5D0)
*                                                                      *
************************************************************************
*                                                                      *
*

      ISYA=MUL(JSYM,ISYU)
      ISYC=MUL(JSYM,ISYL)
      ISYM=ISYU
      ISYAC=MUL(ISYA,ISYC)
C Set up offset table:
      IO1=0
      IO2=0
      DO ISI=1,NSYM
        IOFF1(ISI)=IO1
        IOFF2(ISI)=IO2
        ISAB=MUL(ISI,ISYM)
        IO1=IO1+NISH(ISI)*NAGEB(ISAB)
        IO2=IO2+NISH(ISI)*NAGTB(ISAB)
      END DO

C   Allocate W with parts WP,WM
      NAS=NASH(ISYM)
      NISP=NISUP(ISYM,10)
      NISM=NISUP(ISYM,11)
      NWGP=NAS*NISP
      NWGM=NAS*NISM
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option?
*
*     Incore=nBuff.ge.Max(NWGP,NWGM)+NA*NU*NC*NL
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSG'
*        Call Abend()
*     End If
*                                                                      *
************************************************************************
*                                                                      *
C     CALL DGEMM_('N','T',NA*NU,NC*NL,NCHO,
C    &            1.0D0,Cho_Bra,NA*NU,
C    &                  Cho_Ket,NC*NL,
C    &            0.0D0,AUCL,NA*NU)
*
*     LWG=1+NA*NU*NC*NL
*     LWGP=LWG
*     LWGM=LWG
      LDGP=NAS
      LDGM=NAS
*
* The plus combination:
      IF (NWGP.GT.0) THEN
       ICASE=10
C Read WP:
       CALL RHS_ALLO (NAS,NISP,lg_GP)
       CALL RHS_READ (NAS,NISP,lg_GP,iCASE,iSYM,iVEC)

C to keep memory advantage, scale NL so that NC*NL fits in KCL
C (scaling NL does not affect ordering of integrals = safe)
       NBXSZC=NSECBX
C      NBXSZJ=NINABX
       KCL=NAUCL/(NA*NU)
       NBXSZL=KCL/NC
       IF (NBXSZL.LE.0) THEN
         Write (6,*) 'Not enough memory in ADDRHSG, I give up'
         CALL Abend()
       ENDIF


      DO ICSTA=1,NC,NBXSZC
        ICEND=MIN(ICSTA-1+NBXSZC,NC)
        NCSZ=ICEND-ICSTA+1
        DO ILSTA=1,NL,NBXSZL
          ILEND=MIN(ILSTA-1+NBXSZL,NL)
          NLSZ=ILEND-ILSTA+1

          ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
          CALL DGEMM_('N','T',NA*NU,NCSZ*NLSZ,NCHO,
     &         1.0D0,Cho_Bra,NA*NU,
     &         Cho_Ket(ICLSTA,1),NC*NL,
     &         0.0D0,AUCL,NA*NU)

      ICL=0
      IBUF=0
      DO IL=ILSTA,ILEND
       DO IC=ICSTA,ICEND
        ICABS=IC+NSES(ISYC)
        ICL=ICL+1

        DO IA=1,NA
         IAABS=IA+NSES(ISYA)
         SCL=SQRT(0.5D0)
         IF(IAABS.GE.ICABS) THEN
          IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
          IF(IAABS.EQ.ICABS) SCL=1.0D0
         ELSE
          IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
         END IF
          DO IU=1,NU
           IW1=IU
           IW2=IL+NL*(IAGEC-1)+IOFF1(ISYL)
           IW=IW1+NAS*(IW2-1)
           IBUF=IBUF+1
*          Buff(LWGP-1+IW)=Buff(LWGP-1+IW)+SCL*AUCL(IA,IU,ICL)
           idxBuf(IBUF)=IW
           Buff(IBUF)=SCL*AUCL(IA,IU,ICL)
           IF (IBUF.EQ.NBUFF) THEN
             CALL RHS_SCATTER(LDGP,lg_GP,Buff,idxBuf,IBUF)
             IBUF=0
           END IF
          END DO
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDGP,lg_GP,Buff,idxBuf,IBUF)
       END IF

         ENDDO
       ENDDO

       CALL RHS_SAVE (NAS,NISP,lg_GP,iCASE,iSYM,iVEC)
       CALL RHS_FREE (NAS,NISP,lg_GP)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
* The minus combination:
      IF (NWGM.GT.0) THEN
       ICASE=11
C Read WGM:
       CALL RHS_ALLO (NAS,NISM,lg_GM)
       CALL RHS_READ (NAS,NISM,lg_GM,iCASE,iSYM,iVEC)

C to keep memory advantage, scale NL so that NC*NL fits in KCL
C (scaling NL does not affect ordering of integrals = safe)
       NBXSZC=NSECBX
C      NBXSZJ=NINABX
       KCL=NAUCL/(NA*NU)
       NBXSZL=KCL/NC
       IF (NBXSZL.LE.0) THEN
         Write (6,*) 'Not enough memory in ADDRHSG, I give up'
         CALL Abend()
       ENDIF


      DO ICSTA=1,NC,NBXSZC
        ICEND=MIN(ICSTA-1+NBXSZC,NC)
        NCSZ=ICEND-ICSTA+1
        DO ILSTA=1,NL,NBXSZL
          ILEND=MIN(ILSTA-1+NBXSZL,NL)
          NLSZ=ILEND-ILSTA+1

          ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
          CALL DGEMM_('N','T',NA*NU,NCSZ*NLSZ,NCHO,
     &         1.0D0,Cho_Bra,NA*NU,
     &         Cho_Ket(ICLSTA,1),NC*NL,
     &         0.0D0,AUCL,NA*NU)

       ICL=0
       IBUF=0
       DO IL=ILSTA,ILEND
        DO IC=ICSTA,ICEND
         ICABS=IC+NSES(ISYC)
         ICL=ICL+1

         DO IA=1,NA
          IAABS=IA+NSES(ISYA)
          IF(IAABS.GT.ICABS) THEN
           IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
           SCL=SQRT(1.5D0)
           DO IU=1,NU
            IW1=IU
            IW2=IL+NL*(IAGTC-1)+IOFF2(ISYL)
            IW=IW1+NAS*(IW2-1)
*           Buff(LWGM-1+IW)=Buff(LWGM-1+IW)+SCL*AUCL(IA,IU,ICL)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=SCL*AUCL(IA,IU,ICL)
            IF (IBUF.EQ.NBUFF) THEN
              CALL RHS_SCATTER(LDGM,lg_GM,Buff,idxBuf,IBUF)
              IBUF=0
            END IF
           END DO
          ELSE IF(IAABS.LT.ICABS) THEN
           IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
           SCL=-SQRT(1.5D0)
           DO IU=1,NU
            IW1=IU
            IW2=IL+NL*(IAGTC-1)+IOFF2(ISYL)
            IW=IW1+NAS*(IW2-1)
*           Buff(LWGM-1+IW)=Buff(LWGM-1+IW)+SCL*AUCL(IA,IU,ICL)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=SCL*AUCL(IA,IU,ICL)
            IF (IBUF.EQ.NBUFF) THEN
              CALL RHS_SCATTER(LDGM,lg_GM,Buff,idxBuf,IBUF)
              IBUF=0
            END IF
           END DO
          END IF
         END DO
        END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER(LDGM,lg_GM,Buff,idxBuf,IBUF)
       END IF

         ENDDO
       ENDDO

       CALL RHS_SAVE (NAS,NISM,lg_GM,iCASE,iSYM,iVEC)
       CALL RHS_FREE (NAS,NISM,lg_GM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSH(IVEC,JSYM,ISYJ,ISYL,NA,NJ,NC,NL,AJCL,NAJCL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION AJCL(NC*NL,*)
      DIMENSION Buff(nBuff)
      DIMENSION idxBuf(nBuff)
      DIMENSION Cho_Bra(NA*NJ,NCHO), Cho_Ket(NC*NL,NCHO)
*      Logical Incore
* Case H:
C   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
C   WM(jl,ac)=((ajcl)-(alcj))*SQRT(3.0D0)
*                                                                      *
************************************************************************
*                                                                      *
      IF(ISYJ.LT.ISYL) GOTO 900
      ISYA=MUL(JSYM,ISYJ)
      ISYC=MUL(JSYM,ISYL)
      ISYM=MUL(ISYA,ISYC)
      ISYAC=ISYM
      ISYJL=ISYM

C   Allocate WHP,WHM
      NASP=NAGEB(ISYM)
      NISP=NIGEJ(ISYM)
      NWHP=NASP*NISP
      IF(NWHP.EQ.0) GOTO 900
      NASM=NAGTB(ISYM)
      NISM=NIGTJ(ISYM)
      NWHM=NASM*NISM
*                                                                      *
************************************************************************
*                                                                      *
*     Incore or partitioned option?
*
*     Incore=nBuff.ge.Max(NWHP,NWHM)+NC*NL
*     If (.NOT.Incore) Then
*        Write (6,*) 'Sort out of memory in ADDRHSH'
*        Call Abend()
*     End If

*     nBatch=MIN((nBuff-MAX(NWHP,NWHM))/(NC*NL),NA*NJ)
*     IF((iPrGlb.GE.DEBUG).AND.(nBatch.lt.NA*NJ)) THEN
*       WRITE(6,'(1X,A)') 'less memory than ideal for ADDRHSH:'
*       WRITE(6,'(1X,A12,I12)') 'needed    = ', NA*NJ*NC*NL
*       WRITE(6,'(1X,A12,I12)') 'available = ', nBatch*NC*NL
*     ENDIF

*                                                                      *
************************************************************************
*                                                                      *
*     CALL DGEMM_('N','T',NA*NJ,NC*NL,NCHO,
*    &            1.0D0,Cho_Bra,NA*NJ,
*    &                  Cho_Ket,NC*NL,
*    &            0.0D0,AJCL,NA*NJ)
*
*     LWHP=1+nBatch*NC*NL
*     LWHM=LWHP
      LDHP=NASP
      LDHM=NASM
*
* The plus combination:
      IF (NWHP.GT.0) THEN
       ICASE=12
C Read WP:
       CALL RHS_ALLO (NASP,NISP,lg_HP)
       CALL RHS_READ (NASP,NISP,lg_HP,iCASE,iSYM,iVEC)

C WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
C-SVC20100313: use only part of aj but all of cl, this introduces
C a reordered loop over parts of the RHS array.

C to keep memory advantage, scale NJ so that NA*NJ fits in KAJ
C (scaling NJ or NL does not affect ordering of integrals = safe)
       NBXSZA=NSECBX
C      NBXSZJ=NINABX
       KAJ=NAJCL/(NC*NL)
       NBXSZJ=KAJ/NA
       IF (NBXSZJ.LE.0) THEN
         Write (6,*) 'Not enough memory in ADDRHSH, I give up'
         CALL Abend()
       ENDIF

       NBXSZC=NSECBX
       NBXSZL=NINABX

       DO IASTA=1,NA,NBXSZA
         IAEND=MIN(IASTA-1+NBXSZA,NA)
         NASZ=IAEND-IASTA+1
         DO IJSTA=1,NJ,NBXSZJ
           IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
           NJSZ=IJEND-IJSTA+1

           IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
           CALL DGEMM_('N','T',NC*NL,NASZ*NJSZ,NCHO,
     &          1.0D0,Cho_Ket,NC*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          0.0D0,AJCL,NC*NL)

           DO ICSTA=1,NC,NBXSZC
             ICEND=MIN(ICSTA-1+NBXSZC,NC)
             NCSZ=ICEND-ICSTA+1
             DO ILSTA=1,NL,NBXSZL
               ILEND=MIN(ILSTA-1+NBXSZL,NL)

               ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

       IAJ=0
       IBUF=0
       DO IJ=IJSTA,IJEND
         IJABS=IJ+NIES(ISYJ)
         ILMAX=NL
         IF(ISYJ.EQ.ISYL) ILMAX=IJ
         DO IA=IASTA,IAEND
           IAABS=IA+NSES(ISYA)
           IAJ=IAJ+1

           ICL=0
           DO IL=ILSTA,MIN(ILEND,ILMAX)
             ILABS=IL+NIES(ISYL)
             SCL1=1.0D0
             IJGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
             IF(IJABS.EQ.ILABS) SCL1=SQRT(0.5D0)
             DO IC=ICSTA,ICEND
               ICABS=IC+NSES(ISYC)
               ICL=ICL+1

               SCL=SCL1
               IF(IAABS.GE.ICABS) THEN
                 IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                 IF(IAABS.EQ.ICABS) SCL=SQRT(2.0D0)*SCL1
               ELSE
                 IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
               END IF
               IW=IAGEC+NAGEB(ISYM)*(IJGEL-1)
*              Buff(LWHP-1+IW)=Buff(LWHP-1+IW)+
*    &                          SCL*AJCL(ICLSTA+ICL-1,IAJ)
               IBUF=IBUF+1
               idxBuf(IBUF)=IW
               Buff(IBUF)=SCL*AJCL(ICLSTA+ICL-1,IAJ)
               IF (IBUF.EQ.NBUFF) THEN
                 CALL RHS_SCATTER (LDHP,lg_HP,Buff,idxBuf,IBUF)
                 IBUF=0
               END IF
             END DO
           END DO
         END DO
       END DO
       IF (IBUF.NE.0) THEN
         CALL RHS_SCATTER (LDHP,lg_HP,Buff,idxBuf,IBUF)
       END IF

             ENDDO
           ENDDO
         ENDDO
       ENDDO

       CALL RHS_SAVE (NASP,NISP,lg_HP,iCASE,iSYM,iVEC)
       CALL RHS_FREE (NASP,NISP,lg_HP)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
* The minus combination:
      IF (NWHM.EQ.0) GO TO 900
      ICASE=13
C Read WM:
      CALL RHS_ALLO (NASM,NISM,lg_HM)
      CALL RHS_READ (NASM,NISM,lg_HM,iCASE,iSYM,iVEC)
*
C VM(jl,ac)=((ajcl)-(alcj))*SQRT(3.0D0)
       NBXSZA=NSECBX
C      NBXSZJ=NINABX
       KAJ=NAJCL/(NC*NL)
       NBXSZJ=KAJ/NA
       IF (NBXSZJ.LE.0) THEN
         Write (6,*) 'Not enough memory in ADDRHSH, I give up'
         CALL Abend()
       ENDIF
       NBXSZC=NSECBX
       NBXSZL=NINABX

       DO IASTA=1,NA,NBXSZA
         IAEND=MIN(IASTA-1+NBXSZA,NA)
         NASZ=IAEND-IASTA+1
         DO IJSTA=1,NJ,NBXSZJ
           IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
           NJSZ=IJEND-IJSTA+1

           IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
           CALL DGEMM_('N','T',NC*NL,NASZ*NJSZ,NCHO,
     &          1.0D0,Cho_Ket,NC*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          0.0D0,AJCL,NC*NL)

           DO ICSTA=1,NC,NBXSZC
             ICEND=MIN(ICSTA-1+NBXSZC,NC)
             NCSZ=ICEND-ICSTA+1
             DO ILSTA=1,NL,NBXSZL
               ILEND=MIN(ILSTA-1+NBXSZL,NL)

               ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

      IAJ=0
      IBUF=0
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ.EQ.ISYL) ILMAX=IJ-1
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
          IAJ=IAJ+1

          ICL=0
          DO IL=ILSTA,MIN(ILMAX,ILEND)
            ILABS=IL+NIES(ISYL)
            IJGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
              ICL=ICL+1

              IF (IAABS.GT.ICABS) THEN
                IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL= SQRT(3.0D0)
              ELSE IF(IAABS.LT.ICABS) THEN
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQRT(3.0D0)
              ELSE
                GO TO 700
              ENDIF
              IW=IAGTC+NAGTB(ISYM)*(IJGTL-1)
*             Buff(LWHM-1+IW)=Buff(LWHM-1+IW)+
*    &                         SCL*AJCL(ICLSTA+ICL-1,IAJ)
              IBUF=IBUF+1
              idxBuf(IBUF)=IW
              Buff(IBUF)=SCL*AJCL(ICLSTA+ICL-1,IAJ)
              IF (IBUF.EQ.NBUFF) THEN
                CALL RHS_SCATTER (LDHM,lg_HM,Buff,idxBuf,IBUF)
                IBUF=0
              END IF
 700        CONTINUE
            END DO
          END DO
        END DO
      END DO
      IF (IBUF.NE.0) THEN
        CALL RHS_SCATTER (LDHM,lg_HM,Buff,idxBuf,IBUF)
      END IF

             ENDDO
           ENDDO
         ENDDO
       ENDDO

      CALL RHS_SAVE (NASM,NISM,lg_HM,iCASE,iSYM,iVEC)
      CALL RHS_FREE (NASM,NISM,lg_HM)
*                                                                      *
************************************************************************
*                                                                      *
 900  CONTINUE
      RETURN
      END
