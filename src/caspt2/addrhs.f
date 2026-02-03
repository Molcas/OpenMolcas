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
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYJ,ISYX,NT,NJ,NV,NX,
     &                                nBuff, NCHO
      real(kind=wp) TJVX(NT,NJ,NV,NX)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
*      Logical Incore
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYJ
      IF(NINDEP(ISYM,1).EQ.0) RETURN
      NAS=NTUV(ISYM)
      NIS=NISH(ISYM)
      NWA=NAS*NIS
      IF(NWA.EQ.0) RETURN
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
     &            One,Cho_Bra,NT*NJ,
     &                  Cho_Ket,NV*NX,
     &            Zero,TJVX,NT*NJ)
*
C Compute W(tvx,j)=(tj,vx) + FIMO(t,j)*delta(v,x)/NACTEL
      ICASE=1
*     LWA=1+NT*NJ*NV*NX
      LDA=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_A)
      CALL RHS_READ (NAS,NIS,lg_A,iCASE,iSYM,iVEC)

      if (iParRHS == 1) then
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
*SVC       Buff(LWA-1+IW)=Buff(LWA-1+IW)+TJVX(IT,IJ,IV,IX)
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
#ifdef _MOLCAS_MPP_
      ! Note that iParRHS = 2 only when is_real_par() is true and
      ! PRHS = 2 is specified in the input file (see procinp_caspt2.F90)
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(TJVX,NT*NJ*NV*NX)
       myRank = GA_NodeID()
       CALL GA_Distribution (lg_A,myRank,ILOV,IHIV,JLOV,JHIV)
       if (JLOV > 0) then
        CALL GA_Access(lg_A,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        DO IT=1,NT
         ITABS=IT+NAES(ISYT)
         DO IJ=JLOV,JHIV
          IW2=IJ
          DO IV=1,NV
           IVABS=IV+NAES(ISYV)
           DO IX=1,NX
            IXABS=IX+NAES(ISYX)
            IW1=KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
            if (IW1 >= ILOV .and. IW1 <= IHIV) then
             DBL_MB(MV+IW1-ILOV+LDA*(IW2-JLOV))
     *         = DBL_MB(MV+IW1-ILOV+LDA*(IW2-JLOV))
     *         + TJVX(IT,IJ,IV,IX)
            end if
           END DO
          END DO
         END DO
        END DO
        CALL GA_Release_Update(lg_A,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
      end if

C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_A,iCASE,iSYM,iVEC)
      CALL RHS_FREE(lg_A)
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSA

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSB(IVEC,JSYM,ISYJ,ISYL,NT,NJ,NV,NL,TJVL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Half, Two, Quart
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
Case B:
* t>v j>l WP(tv,jl)=add ((tj,vl))*(1/2)
* t>v j=l WP(tv,jl)=add ((tj,vl))*(1/2)*(SQRT(2))
* t>v j<l WP(tv,lj)=add ((tj,vl))*(1/2)

* t=v j>l WP(tv,jl)=add ((tj,vl))*(1/4)
* t=v j=l WP(tv,jl)=add ((tj,vl))*(1/4)*(SQRT(2))
* t=v j<l WP(tv,lj)=add ((tj,vl))*(1/4)

      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYJ,ISYL,NT,NJ,NV,NL,
     &                                nBuff, NCHO
      real(kind=wp) TJVL(NT,NJ,NV,NL)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      IF(ISYT.LT.ISYV) RETURN
      SQ2=SQRT(Two)
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
      If (Max(NWBP,NWBM).le.0) RETURN
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
     &            One,Cho_Bra,NT*NJ,
     &                  Cho_Ket,NV*NL,
     &            Zero,TJVL,NT*NJ)
#ifdef _MOLCAS_MPP_
      if (iParRHS == 2) call GADSUM_ADDRHS(TJVL,NT*NJ*NV*NL)
#endif
      If (NWBP>0) THEN
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

       if (iParRHS == 1) then
        IBUF=0
         DO IT=1,NT
         ITABS=IT+NAES(ISYT)
         IVMAX=NV
         IF(ISYV.EQ.ISYT) IVMAX=IT
         DO IV=1,IVMAX
          IVABS=IV+NAES(ISYV)
          SCL1=Half
          IW1=KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
          IF(ITABS.EQ.IVABS) SCL1=Quart
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
*           Buff(LWBP-1+IW)=Buff(LWBP-1+IW)+SCL*TJVL(IT,IJ,IV,IL)
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
#ifdef _MOLCAS_MPP_
       else if (iParRHS == 2) then
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_BP,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) then
         CALL GA_Access(lg_BP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
         DO IT=1,NT
          ITABS=IT+NAES(ISYT)
          IVMAX=NV
          IF(ISYV.EQ.ISYT) IVMAX=IT
          DO IV=1,IVMAX
           IVABS=IV+NAES(ISYV)
           SCL1=Half
           IW1=KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
           if (IW1 < ILOV .or. IW1 > IHIV) cycle
           IF(ITABS.EQ.IVABS) SCL1=Quart
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
             if (IW2 >= JLOV .and. IW2 <= JHIV) then
              DBL_MB(MV+IW1-ILOV+LDBP*(IW2-JLOV))
     *          = DBL_MB(MV+IW1-ILOV+LDBP*(IW2-JLOV))
     *          + SCL*TJVL(IT,IJ,IV,IL)
             end if
            END DO
           END DO
          END DO
         END DO
         CALL GA_Release_Update(lg_BP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        end if
#endif
       end if

C Put WBP on disk:
       CALL RHS_SAVE (NASP,NISP,lg_BP,iCASE,iSYM,iVEC)
       CALL RHS_FREE(lg_BP)
      END IF
      END IF
*                                                                      *
************************************************************************
*                                                                      *
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

      if (iParRHS == 1) then
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
*           Buff(LWBM-1+IW)=Buff(LWBM-1+IW)+Half*TJVL(IT,IJ,IV,IL)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=Half*TJVL(IT,IJ,IV,IL)
           ELSE IF (IJABS.LT.ILABS) THEN
            IW2=KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
            IW=IW1+NASM*(IW2-1)
*           Buff(LWBM-1+IW)=Buff(LWBM-1+IW)-Half*TJVL(IT,IJ,IV,IL)
            IBUF=IBUF+1
            idxBuf(IBUF)=IW
            Buff(IBUF)=-Half*TJVL(IT,IJ,IV,IL)
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       myRank = GA_NodeID()
       CALL GA_Distribution (lg_BM,myRank,ILOV,IHIV,JLOV,JHIV)
       if (JLOV > 0) then
        CALL GA_Access(lg_BM,ILOV,IHIV,JLOV,JHIV,MV,LDV)

        DO IT=1,NT
         ITABS=IT+NAES(ISYT)
         IVMAX=NV
         IF(ISYV.EQ.ISYT) IVMAX=IT-1
         DO IV=1,IVMAX
          IVABS=IV+NAES(ISYV)
          IW1=KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
          if (IW1 < ILOV .or. IW1 > IHIV) cycle
          DO IJ=1,NJ
           IJABS=IJ+NIES(ISYJ)
           DO IL=1,NL
            ILABS=IL+NIES(ISYL)
            IF(IJABS.GT.ILABS) THEN
             IW2=KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
             if (IW2 >= JLOV .and. IW2 <= JHIV) then
               DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV))
     *           = DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV))
     *           + Half*TJVL(IT,IJ,IV,IL)
             end if
            ELSE IF (IJABS.LT.ILABS) THEN
             IW2=KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
             if (IW2 >= JLOV .and. IW2 <= JHIV) then
               DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV))
     *           = DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV))
     *           - Half*TJVL(IT,IJ,IV,IL)
             end if
            END IF
            IF (IBUF.EQ.NBUFF) THEN
              CALL RHS_SCATTER(LDBM,lg_BM,Buff,idxBuf,IBUF)
              IBUF=0
            END IF
           END DO
          END DO
         END DO
        END DO
        CALL GA_Release_Update(lg_BM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
      end if

C Put WBM on disk:
       CALL RHS_SAVE (NASM,NISM,lg_BM,iCASE,iSYM,iVEC)
       CALL RHS_FREE (lg_BM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSB

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSC(IVEC,JSYM,ISYU,ISYX,NA,NU,NV,NX,AUVX,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYU,ISYX,NA,NU,NV,NX,
     &                                nBuff, NCHO
      real(kind=wp) AUVX(NA,NU,NV,NX)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NX,NCHO)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
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
      IF(NINDEP(ISYM,4).EQ.0) RETURN
      NAS=NTUV(ISYM)
      NIS=NSSH(ISYM)
      NWC=NAS*NIS
      IF(NWC.EQ.0) RETURN
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
     &            One,Cho_Bra,NA*NU,
     &                  Cho_Ket,NV*NX,
     &            Zero,AUVX,NA*NU)
*
      ICASE=4
*     LWC=1+NA*NU*NV*NX
      LDC=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_C)
      CALL RHS_READ (NAS,NIS,lg_C,iCASE,iSYM,iVEC)
      if (iParRHS == 1) then
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
*          Buff(LWC-1+IW)=Buff(LWC-1+IW)+AUVX(IA,IU,IV,IX)
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AUVX,NA*NU*NV*NX)
       myRank = GA_NodeID()
       CALL GA_Distribution (lg_C,myRank,ILOV,IHIV,JLOV,JHIV)
       if (JLOV > 0) then
        CALL GA_Access(lg_C,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        DO IA=JLOV,JHIV
         IW2=IA
         DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          DO IV=1,NV
           IVABS=IV+NAES(ISYV)
           DO IX=1,NX
            IXABS=IX+NAES(ISYX)
            IW1=KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
            if (IW1 >= ILOV .and. IW1 <= IHIV) then
             DBL_MB(MV+IW1-ILOV+LDC*(IW2-JLOV))
     *         = DBL_MB(MV+IW1-ILOV+LDC*(IW2-JLOV))
     *         + AUVX(IA,IU,IV,IX)
            end if
           END DO
          END DO
         END DO
        END DO
        CALL GA_Release_Update(lg_C,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
      end if
C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_C,iCASE,iSYM,iVEC)
      CALL RHS_FREE (lg_C)

*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSC

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSD1(IVEC,JSYM,ISYJ,ISYX,NA,NJ,NV,NX,AJVX,
     &                    nBuff,Buff,idxBuf,
     &                    Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYJ,ISYX,NA,NJ,NV,NX,
     &                                nBuff, NCHO
      real(kind=wp) AJVX(NV,NX,*)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV,NX,NCHO)
*      Logical Incore
      integer(kind=iwp) IOFFD(8,8)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
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
      IF(NINDEP(ISYM,5).EQ.0) RETURN
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) RETURN
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
C    &            One,Cho_Bra,NA*NJ,
C    &                  Cho_Ket,NV*NX,
C    &            Zero,AJVX,NA*NJ)
*

C Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
      ICASE=5
*     LWD=1+NA*NJ*NV*NX
      LDD=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_D)
      CALL RHS_READ (NAS,NIS,lg_D,iCASE,iSYM,iVEC)
#ifdef _MOLCAS_MPP_
      if (iParRHS == 2) then
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_D,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) CALL GA_Access(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif

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
     &         One,Cho_Ket,NV*NX,
     &         Cho_Bra(IAJSTA,1),NA*NJ,
     &         Zero,AJVX,NV*NX)

      if (iParRHS == 1) then
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
*          WD=Buff(LWD-1+IW)+AJVX(IV,IX,IAJ)
*          Buff(LWD-1+IW)=WD
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AJVX,NV*NX*NASZ*NJSZ)
       if (JLOV > 0) then
        IAJ=0
        DO IJ=IJSTA,IJEND
         DO IA=IASTA,IAEND
         IAJ=IAJ+1
          IW2=IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
          if (IW2 < JLOV .or. IW2 > JHIV) cycle

          DO IX=1,NX
           IXABS=IX+NAES(ISYX)
            DO IV=1,NV
            IVABS=IV+NAES(ISYV)

            IW1=KTU(IVABS,IXABS)-NTUES(ISYM)
            if (IW1 >= ILOV .and. IW1 <= IHIV) then
             DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV))
     *         = DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV))
     *         + AJVX(IV,IX,IAJ)
            end if
           END DO
          END DO
         END DO
        END DO
       end if
#endif
      end if

        ENDDO
      ENDDO

#ifdef _MOLCAS_MPP_
      if (iParRHS == 2 .and. JLOV > 0) then
        CALL GA_Release_Update(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif
C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_D,iCASE,iSYM,iVEC)
      CALL RHS_FREE (lg_D)
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSD1

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSD2(IVEC,JSYM,ISYU,ISYL,NA,NU,NV,NL,AUVL,
     &                    nBuff,Buff,idxBuf,
     &                    Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYU,ISYL,NA,NU,NV,NL,
     &                                nBuff, NCHO
      real(kind=wp) AUVL(NA,NU,NV,NL)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
      integer(kind=iwp) IOFFD(8,8)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
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
      IF(NINDEP(ISYM,5).EQ.0) RETURN
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) RETURN
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
     &            One,Cho_Bra,NA*NU,
     &                  Cho_Ket,NV*NL,
     &            Zero,AUVL,NA*NU)
*
C Compute W2(vu,al)=(au,vl)
      ICASE=5
*     LWD=1+NA*NU*NV*NL
      LDD=NAS
C Read W:
      CALL RHS_ALLO (NAS,NIS,lg_D)
      CALL RHS_READ (NAS,NIS,lg_D,iCASE,iSYM,iVEC)

      if (iParRHS == 1) then
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
*          Buff(LWD-1+IW)=Buff(LWD-1+IW)+AUVL(IA,IU,IV,IL)
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AUVL,NA*NU*NV*NL)
       myRank = GA_NodeID()
       CALL GA_Distribution (lg_D,myRank,ILOV,IHIV,JLOV,JHIV)
       if (JLOV > 0) then
        CALL GA_Access(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        DO IA=1,NA
         DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          DO IV=1,NV
           IVABS=IV+NAES(ISYV)
           IW1=NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
           if (IW1 < ILOV .and. IW1 > IHIV) cycle
           DO IL=1,NL
            IW2=IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
            if (IW2 >= JLOV .and. IW2 <= JHIV) then
              DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV))
     *          = DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV))
     *          + AUVL(IA,IU,IV,IL)
            end if
           END DO
          END DO
         END DO
        END DO
        CALL GA_Release_Update(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
      end if
*
C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_D,iCASE,iSYM,iVEC)
      CALL RHS_FREE (lg_D)
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSD2

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSE(IVEC,JSYM,ISYJ,ISYL,NA,NJ,NV,NL,AJVL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Half, OneHalf
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYJ,ISYL,NA,NJ,NV,NL,
     &                                nBuff, NCHO
      real(kind=wp) AJVL(NV,NL,*)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
      integer(kind=iwp) IOFF1(8),IOFF2(8)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
Case E:
*                                                                      *
************************************************************************
*                                                                      *
*
      SQ32=SQRT(OneHalf)
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
      If (NW.eq.0) RETURN
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
C    &            One,Cho_Bra,NA*NJ,
C    &                  Cho_Ket,NV*NL,
C    &            Zero,AJVL,NA*NJ)
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
#ifdef _MOLCAS_MPP_
      if (iParRHS == 2) then
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_EP,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) CALL GA_Access(lg_EP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif

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
     &          One,Cho_Ket,NV*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          Zero,AJVL,NV*NL)
      if (iParRHS == 1) then
       IAJ=0
       IBUF=0
       DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        DO IA=IASTA,IAEND
        IAJ=IAJ+1

         DO IV=1,NV
          DO IL=1,NL
           ILABS=IL+NIES(ISYL)
           SCL=SQRT(Half)
           IF(IJABS.GE.ILABS) THEN
            JGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
            IF(IJABS.EQ.ILABS) SCL=One
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
#ifdef _MOLCAS_MPP_
      else
       call GADSUM_ADDRHS(AJVL,NV*NL*NASZ*NJSZ)
       if (JLOV > 0) then
        IAJ=0
        DO IJ=IJSTA,IJEND
         IJABS=IJ+NIES(ISYJ)
         DO IA=IASTA,IAEND
         IAJ=IAJ+1

          DO IV=ILOV,IHIV
           IW1=IV
           DO IL=1,NL
            ILABS=IL+NIES(ISYL)
            SCL=SQRT(Half)
            IF(IJABS.GE.ILABS) THEN
             JGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
             IF(IJABS.EQ.ILABS) SCL=One
            ELSE
             JGEL=KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
            END IF
            IW2=IA+NA*(JGEL-1)+IOFF1(ISYA)
            if (IW2 >= JLOV .and. IW2 <= JHIV) then
              DBL_MB(MV+IW1-ILOV+LDEP*(IW2-JLOV))
     *          = DBL_MB(MV+IW1-ILOV+LDEP*(IW2-JLOV))
     *          + SCL*AJVL(IV,IL,IAJ)
            end if
           END DO
          END DO
         END DO
        END DO
       end if
#endif
      end if

         ENDDO
       ENDDO

#ifdef _MOLCAS_MPP_
      if (iParRHS == 2 .and. JLOV > 0) then
        CALL GA_Release_Update(lg_EP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif
      CALL RHS_SAVE (NAS,NISP,lg_EP,iCASE,iSYM,iVEC)
      CALL RHS_FREE (lg_EP)
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
#ifdef _MOLCAS_MPP_
      if (iParRHS == 2) then
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_EM,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) CALL GA_Access(lg_EM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif

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
     &          One,Cho_Ket,NV*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          Zero,AJVL,NV*NL)

      if (iParRHS == 1) then
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AJVL,NV*NL*NASZ*NJSZ)
       if (JLOV > 0) then
        IAJ=0
        DO IJ=IJSTA,IJEND
         IJABS=IJ+NIES(ISYJ)
         DO IA=IASTA,IAEND
         IAJ=IAJ+1

          DO IV=ILOV,IHIV
           IW1=IV
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
             IW2=IA+NA*(JGTL-1)+IOFF2(ISYA)
             if (IW2 >= JLOV .and. IW2 <= JHIV) then
               DBL_MB(MV+IW1-ILOV+LDEM*(IW2-JLOV))
     *           = DBL_MB(MV+IW1-ILOV+LDEM*(IW2-JLOV))
     *           + SCL*AJVL(IV,IL,IAJ)
             end if
            END IF
           END DO
          END DO
         END DO
        END DO
       end if
#endif
      end if

         ENDDO
       ENDDO

#ifdef _MOLCAS_MPP_
      if (iParRHS == 2 .and. JLOV > 0) then
        CALL GA_Release_Update(lg_EM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif
      CALL RHS_SAVE (NAS,NISM,lg_EM,iCASE,iSYM,iVEC)
      CALL RHS_FREE (lg_EM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSE

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSF(IVEC,JSYM,ISYU,ISYX,NA,NU,NC,NX,AUCX,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Half, Two, Quart
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYU,ISYX,NA,NU,NC,NX,
     &                                nBuff, NCHO
      real(kind=wp) AUCX(NA,NU,NC,NX)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
*      Logical Incore
Case F:
C   WP(ux,ac)=((aucx)+(axcu))*(1-Kron(x,u)/2) /2
C With new normalisation, replace /2 with /(2*SQRT(1+Kron(ac))
C   WM(ux,ac)= -((aucx)-(axcu))/2
*                                                                      *
************************************************************************
*                                                                      *
*
      IF(ISYU.LT.ISYX) RETURN

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
      If (NWFP+NWFM.le.0) RETURN
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
     &            One,Cho_Bra,NA*NU,
     &                  Cho_Ket,NC*NX,
     &            Zero,AUCX,NA*NU)
#ifdef _MOLCAS_MPP_
      IF (iParRHS == 2) call GADSUM_ADDRHS(AUCX,NA*NU*NC*NX)
#endif
*
      IF (NWFP>0) THEN
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

       if (iParRHS == 1) then
        IBUF=0
        DO IU=1,NU
         IUABS=IU+NAES(ISYU)
         IXMAX=NX
         IF(ISYU.EQ.ISYX) IXMAX=IU
         DO IX=1,IXMAX
          IXABS=IX+NAES(ISYX)
          SCL1=Half
          IF(IUABS.EQ.IXABS) SCL1=Quart
          IW1=KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
          DO IA=1,NA
           IAABS=IA+NSES(ISYA)
           DO IC=1,NC
            ICABS=IC+NSES(ISYC)
            SCL=SCL1
            IF(IAABS.GE.ICABS) THEN
             IW2=KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
             IF(IAABS.EQ.ICABS) SCL=SQRT(Two)*SCL1
            ELSE
             IW2=KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
            END IF
            IW=IW1+NASP*(IW2-1)
*           WFP=Buff(LWFP-1+IW)+SCL*AUCX(IA,IU,IC,IX)
*           Buff(LWFP-1+IW)=WFP
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
#ifdef _MOLCAS_MPP_
       else if (iParRHS == 2) then
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_FP,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) then
         CALL GA_Access(lg_FP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
         IBUF=0
         DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          IXMAX=NX
          IF(ISYU.EQ.ISYX) IXMAX=IU
          DO IX=1,IXMAX
           IXABS=IX+NAES(ISYX)
           SCL1=Half
           IF(IUABS.EQ.IXABS) SCL1=Quart
           IW1=KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
           if (IW1 < ILOV .or. IW1 > IHIV) cycle
           DO IA=1,NA
            IAABS=IA+NSES(ISYA)
            DO IC=1,NC
             ICABS=IC+NSES(ISYC)
             SCL=SCL1
             IF(IAABS.GE.ICABS) THEN
              IW2=KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
              IF(IAABS.EQ.ICABS) SCL=SQRT(Two)*SCL1
             ELSE
              IW2=KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
             END IF
             if (IW2 >= JLOV .and. IW2 <= JHIV) then
              DBL_MB(MV+IW1-ILOV+LDFP*(IW2-JLOV))
     *        = DBL_MB(MV+IW1-ILOV+LDFP*(IW2-JLOV))
     *        + SCL*AUCX(IA,IU,IC,IX)
             end if
            END DO
           END DO
          END DO
         END DO
         CALL GA_Release_Update(lg_FP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        end if
#endif
       end if

       CALL RHS_SAVE (NASP,NISP,lg_FP,iCASE,iSYM,iVEC)
       CALL RHS_FREE (lg_FP)
      END IF
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      IF (NWFM.le.0) RETURN
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
#ifdef _MOLCAS_MPP_
       if (Is_Real_Par()) then
         myRank = GA_NodeID()
         CALL GA_Distribution (lg_FM,myRank,ILOV,IHIV,JLOV,JHIV)
         if (JLOV > 0) CALL GA_Access(lg_FM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif

       if (iParRHS == 1) then
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
*            Buff(LWFM-1+IW)=Buff(LWFM-1+IW)-Half*AUCX(IA,IU,IC,IX)
             IBUF=IBUF+1
             idxBuf(IBUF)=IW
             Buff(IBUF)=-Half*AUCX(IA,IU,IC,IX)
            ELSE IF(IAABS.LT.ICABS) THEN
             IW2=KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
             IW=IW1+NASM*(IW2-1)
*            Buff(LWFM-1+IW)=Buff(LWFM-1+IW)+Half*AUCX(IA,IU,IC,IX)
             IBUF=IBUF+1
             idxBuf(IBUF)=IW
             Buff(IBUF)=Half*AUCX(IA,IU,IC,IX)
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
#ifdef _MOLCAS_MPP_
       else if (iParRHS == 2) then
        if (JLOV > 0) then
         DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          IXMAX=NX
          IF(ISYU.EQ.ISYX) IXMAX=IU-1
          DO IX=1,IXMAX
           IXABS=IX+NAES(ISYX)
           IW1=KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
           if (IW1 < ILOV .or. IW1 > IHIV) cycle
           DO IA=1,NA
            IAABS=IA+NSES(ISYA)
            DO IC=1,NC
             ICABS=IC+NSES(ISYC)
             IF(IAABS.GT.ICABS) THEN
              IW2=KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
              if (IW2 >= JLOV .and. IW2 <= JHIV) then
               DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV))
     *           = DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV))
     *           - Half*AUCX(IA,IU,IC,IX)
              end if
             ELSE IF(IAABS.LT.ICABS) THEN
              IW2=KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
              if (IW2 >= JLOV .and. IW2 <= JHIV) then
               DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV))
     *           = DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV))
     *           + Half*AUCX(IA,IU,IC,IX)
              end if
             END IF
            END DO
           END DO
          END DO
         END DO
         CALL GA_Release_Update(lg_FM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        end if
#endif
       end if

C Put WFM on disk:
       CALL RHS_SAVE (NASM,NISM,lg_FM,iCASE,iSYM,iVEC)
       CALL RHS_FREE (lg_FM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSF

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSG(IVEC,JSYM,ISYU,ISYL,NA,NU,NC,NL,AUCL,NAUCL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Half, OneHalf
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYU,ISYL,NA,NU,NC,NL,
     &                                NAUCL, nBuff, NCHO
      real(kind=wp) AUCL(NA,NU,*)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA,NU,NCHO), Cho_Ket(NC*NL,NCHO)
*      Logical Incore
      integer(kind=iwp) IOFF1(8),IOFF2(8)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV
#endif
Case G:
C   WP(u,l,ac)=  ((aucl)+cual))/SQRT(2+2*Kron(ab))
C   WM(u,l,ac)=  ((aucl)-cual))*SQRT(OneHalf)
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
C    &            One,Cho_Bra,NA*NU,
C    &                  Cho_Ket,NC*NL,
C    &            Zero,AUCL,NA*NU)
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
#ifdef _MOLCAS_MPP_
       if (iParRHS == 2) then
         myRank = GA_NodeID()
         CALL GA_Distribution (lg_GP,myRank,ILOV,IHIV,JLOV,JHIV)
         if (JLOV > 0) CALL GA_Access(lg_GP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif

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
     &         One,Cho_Bra,NA*NU,
     &         Cho_Ket(ICLSTA,1),NC*NL,
     &         Zero,AUCL,NA*NU)

      if (iParRHS == 1) then
       ICL=0
       IBUF=0
       DO IL=ILSTA,ILEND
        DO IC=ICSTA,ICEND
         ICABS=IC+NSES(ISYC)
         ICL=ICL+1

         DO IA=1,NA
          IAABS=IA+NSES(ISYA)
          SCL=SQRT(Half)
          IF(IAABS.GE.ICABS) THEN
           IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
           IF(IAABS.EQ.ICABS) SCL=One
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AUCL,NA*NU*NCSZ*NLSZ)
       if (JLOV > 0) then
        ICL=0
        DO IL=ILSTA,ILEND
         DO IC=ICSTA,ICEND
          ICABS=IC+NSES(ISYC)
          ICL=ICL+1

          DO IA=1,NA
           IAABS=IA+NSES(ISYA)
           SCL=SQRT(Half)
           IF(IAABS.GE.ICABS) THEN
            IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
            IF(IAABS.EQ.ICABS) SCL=One
           ELSE
            IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
           END IF
           IW2 = IL+NL*(IAGEC-1)+IOFF1(ISYL)
           if (IW2 < JLOV .or. IW2 > JHIV) cycle
           DO IU=1,NU
            IW1 = IU
            if (IW1 >= ILOV .and. IW1 <= IHIV) then
              DBL_MB(MV+IW1-ILOV+LDGP*(IW2-JLOV))
     *          = DBL_MB(MV+IW1-ILOV+LDGP*(IW2-JLOV))
     *          + SCL*AUCL(IA,IU,ICL)
            end if
           END DO
          END DO
         END DO
        END DO
       end if
#endif
      end if

         ENDDO
       ENDDO

#ifdef _MOLCAS_MPP_
       if (iParRHS == 2 .and. JLOV > 0) then
         CALL GA_Release_Update(lg_GP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
       CALL RHS_SAVE (NAS,NISP,lg_GP,iCASE,iSYM,iVEC)
       CALL RHS_FREE (lg_GP)
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
#ifdef _MOLCAS_MPP_
       if (iParRHS == 2) then
         myRank = GA_NodeID()
         CALL GA_Distribution (lg_GM,myRank,ILOV,IHIV,JLOV,JHIV)
         if (JLOV > 0) CALL GA_Access(lg_GM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif

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
     &         One,Cho_Bra,NA*NU,
     &         Cho_Ket(ICLSTA,1),NC*NL,
     &         Zero,AUCL,NA*NU)

      if (iParRHS == 1) then
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
           SCL=SQRT(OneHalf)
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
           SCL=-SQRT(OneHalf)
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AUCL,NA*NU*NCSZ*NLSZ)
       if (JLOV > 0) then
        ICL=0
        DO IL=ILSTA,ILEND
         DO IC=ICSTA,ICEND
          ICABS=IC+NSES(ISYC)
          ICL=ICL+1

          DO IA=1,NA
           IAABS=IA+NSES(ISYA)
           IF(IAABS.GT.ICABS) THEN
            IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
            ITMP2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
            if (ITMP2 < JLOV .or. ITMP2 > JHIV) cycle
            SCL=SQRT(OneHalf)
            DO IU=1,NU
              ITMP1 = IU
              if (ITMP1 >= ILOV .and. ITMP1 <= IHIV) then
                DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV))
     *            = DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV))
     *            + SCL*AUCL(IA,IU,ICL)
              end if
            END DO
           ELSE IF(IAABS.LT.ICABS) THEN
            IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
            ITMP2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
            if (ITMP2 < JLOV .or. ITMP2 > JHIV) cycle
            SCL=-SQRT(OneHalf)
            DO IU=1,NU
              ITMP1 = IU
              if (ITMP1 >= ILOV .and. ITMP1 <= IHIV) then
                DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV))
     *            = DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV))
     *            + SCL*AUCL(IA,IU,ICL)
              end if
            END DO
           END IF
          END DO
         END DO
        END DO
       end if
#endif
      end if

         ENDDO
       ENDDO

#ifdef _MOLCAS_MPP_
       if (iParRHS == 2 .and. JLOV > 0) then
         CALL GA_Release_Update(lg_GM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
       CALL RHS_SAVE (NAS,NISM,lg_GM,iCASE,iSYM,iVEC)
       CALL RHS_FREE (lg_GM)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSG

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE ADDRHSH(IVEC,JSYM,ISYJ,ISYL,NA,NJ,NC,NL,AJCL,NAJCL,
     &                   nBuff,Buff,idxBuf,
     &                   Cho_Bra,Cho_Ket,NCHO)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Half, Two, Three
      use caspt2_global, only: iParRHS
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), Intent(in):: IVEC,JSYM,ISYJ,ISYL,NA,NJ,NC,NL,
     &                                NAJCL, nBuff, NCHO
      real(kind=wp) AJCL(NC*NL,*)
      real(kind=wp) Buff(nBuff)
      integer(kind=iwp) idxBuf(nBuff)
      real(kind=wp) Cho_Bra(NA*NJ,NCHO), Cho_Ket(NC*NL,NCHO)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank,ILOV,IHIV,JLOV,JHIV,MV,LDV,ITMP1,ITMP2
#endif
*      Logical Incore
* Case H:
C   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
C   WM(jl,ac)=((ajcl)-(alcj))*SQRT(Three)
*                                                                      *
************************************************************************
*                                                                      *
      IF(ISYJ.LT.ISYL) RETURN
      ISYA=MUL(JSYM,ISYJ)
      ISYC=MUL(JSYM,ISYL)
      ISYM=MUL(ISYA,ISYC)
      ISYAC=ISYM
      ISYJL=ISYM

C   Allocate WHP,WHM
      NASP=NAGEB(ISYM)
      NISP=NIGEJ(ISYM)
      NWHP=NASP*NISP
      IF(NWHP.EQ.0) RETURN
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
*    &            One,Cho_Bra,NA*NJ,
*    &                  Cho_Ket,NC*NL,
*    &            Zero,AJCL,NA*NJ)
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

#ifdef _MOLCAS_MPP_
       if (iParRHS == 2) then
         myRank = GA_NodeID()
         CALL GA_Distribution (lg_HP,myRank,ILOV,IHIV,JLOV,JHIV)
         if (JLOV > 0) CALL GA_Access(lg_HP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif

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
     &          One,Cho_Ket,NC*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          Zero,AJCL,NC*NL)

      if (iParRHS == 1) then
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
             SCL1=One
             IJGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
             IF(IJABS.EQ.ILABS) SCL1=SQRT(Half)
             DO IC=ICSTA,ICEND
               ICABS=IC+NSES(ISYC)
               ICL=ICL+1

               SCL=SCL1
               IF(IAABS.GE.ICABS) THEN
                 IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                 IF(IAABS.EQ.ICABS) SCL=SQRT(Two)*SCL1
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
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
       call GADSUM_ADDRHS(AJCL,NC*NL*NASZ*NJSZ)
        if (JLOV > 0) then

           DO ICSTA=1,NC,NBXSZC
             ICEND=MIN(ICSTA-1+NBXSZC,NC)
             NCSZ=ICEND-ICSTA+1
             DO ILSTA=1,NL,NBXSZL
               ILEND=MIN(ILSTA-1+NBXSZL,NL)

               ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

       IAJ=0
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
             SCL1=One
             IJGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
             if (IJGEL < JLOV .or. IJGEL > JHIV) then
               ICL = ICL + ICEND-ICSTA+1
               cycle
             end if
             IF(IJABS.EQ.ILABS) SCL1=SQRT(Half)
             DO IC=ICSTA,ICEND
               ICABS=IC+NSES(ISYC)
               ICL=ICL+1

               SCL=SCL1
               IF(IAABS.GE.ICABS) THEN
                 IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                 IF(IAABS.EQ.ICABS) SCL=SQRT(Two)*SCL1
               ELSE
                 IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
               END IF
               if (IAGEC >= ILOV .and. IAGEC <= IHIV) then
                 ITMP1 = IAGEC
                 ITMP2 = IJGEL
                 DBL_MB(MV+ITMP1-ILOV+LDHP*(ITMP2-JLOV))
     *             = DBL_MB(MV+ITMP1-ILOV+LDHP*(ITMP2-JLOV))
     *             + SCL*AJCL(ICLSTA+ICL-1,IAJ)
               end if
             END DO
           END DO
         END DO
       END DO
             ENDDO
           ENDDO
         end if
#endif
       end if
         ENDDO
       ENDDO

#ifdef _MOLCAS_MPP_
       if (iParRHS == 2 .and. JLOV > 0) then
         CALL GA_Release_Update(lg_HP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
       end if
#endif
       CALL RHS_SAVE (NASP,NISP,lg_HP,iCASE,iSYM,iVEC)
       CALL RHS_FREE (lg_HP)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
* The minus combination:
      IF (NWHM.EQ.0) RETURN
      ICASE=13
C Read WM:
      CALL RHS_ALLO (NASM,NISM,lg_HM)
      CALL RHS_READ (NASM,NISM,lg_HM,iCASE,iSYM,iVEC)

#ifdef _MOLCAS_MPP_
      if (iParRHS == 2) then
        myRank = GA_NodeID()
        CALL GA_Distribution (lg_HM,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) CALL GA_Access(lg_HM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif
*
C VM(jl,ac)=((ajcl)-(alcj))*SQRT(Three)
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
     &          One,Cho_Ket,NC*NL,
     &          Cho_Bra(IAJSTA,1),NA*NJ,
     &          Zero,AJCL,NC*NL)

      if (iParRHS == 1) then
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
                SCL= SQRT(Three)
              ELSE IF(IAABS.LT.ICABS) THEN
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQRT(Three)
              ELSE
                CYCLE
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
            END DO
          END DO
        END DO
      END DO
      IF (IBUF.NE.0) THEN
        CALL RHS_SCATTER (LDHM,lg_HM,Buff,idxBuf,IBUF)
      END IF

             ENDDO
           ENDDO
#ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
        call GADSUM_ADDRHS(AJCL,NC*NL*NASZ*NJSZ)
        if (JLOV > 0) then

           DO ICSTA=1,NC,NBXSZC
             ICEND=MIN(ICSTA-1+NBXSZC,NC)
             NCSZ=ICEND-ICSTA+1
             DO ILSTA=1,NL,NBXSZL
               ILEND=MIN(ILSTA-1+NBXSZL,NL)

               ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

      IAJ=0
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
            if (IJGTL < JLOV .or. IJGTL > JHIV) then
              ICL = ICL + ICEND-ICSTA+1
              cycle
            end if
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
              ICL=ICL+1

              IF (IAABS.GT.ICABS) THEN
                IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL= SQRT(Three)
              ELSE IF(IAABS.LT.ICABS) THEN
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQRT(Three)
              ELSE
                cycle
              ENDIF
              if (IAGTC >= ILOV .and. IAGTC <= IHIV) then
                ITMP1 = IAGTC
                ITMP2 = IJGTL
                DBL_MB(MV+ITMP1-ILOV+LDHM*(ITMP2-JLOV))
     *            = DBL_MB(MV+ITMP1-ILOV+LDHM*(ITMP2-JLOV))
     *            + SCL*AJCL(ICLSTA+ICL-1,IAJ)
              end if
            END DO
          END DO
        END DO
      END DO
             ENDDO
           ENDDO
         end if
#endif
      end if
         ENDDO
       ENDDO

#ifdef _MOLCAS_MPP_
      if (iParRHS == 2 .and. JLOV > 0) then
        CALL GA_Release_Update(lg_HM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#endif
      CALL RHS_SAVE (NASM,NISM,lg_HM,iCASE,iSYM,iVEC)
      CALL RHS_FREE (lg_HM)
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE ADDRHSH

      subroutine GADSUM_ADDRHS(buff,nbuff)
      use definitions, only: iwp, wp
      use caspt2_global, only: MAXBUF
      use definitions, only: iwp,wp
      implicit none
      integer(kind=iwp), intent(in) :: nbuff
      real(kind=wp), intent(inout) :: buff(nbuff)
      integer(kind=iwp) :: istart
      ! GADSUM wrapper: avoid the 2 GB limit of 32-bit MPI
      do istart = 1, nbuff, MAXBUF
       CALL GADSUM(buff(istart),MIN(nbuff-istart+1,MAXBUF))
      end do
      end subroutine GADSUM_ADDRHS
