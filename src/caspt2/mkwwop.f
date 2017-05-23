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
      SUBROUTINE MKWWOP(IVEC,JVEC,OP0,OP1,NOP2,OP2,NOP3,OP3)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

C Presently symmetry blocking is disregarded for OP2, OP3, but
C index pair C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      DIMENSION OP1(NASHT,NASHT),OP2(NOP2),OP3(NOP3)

C Given the coefficients for two excitation operators in the
C vectors numbered IVEC and C JVEC on file, construct the
C zero-, one-, two-, and three-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.

      OP0=0.0D0
      CALL DCOPY_(NASHT**2,0.0D0,0,OP1,1)
      CALL DCOPY_(NOP2,0.0D0,0,OP2,1)
      CALL DCOPY_(NOP3,0.0D0,0,OP3,1)
      CALL MKWWOPA(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      CALL MKWWOPB(IVEC,JVEC,OP0,OP1,NOP2,OP2)
      CALL MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      CALL MKWWOPD(IVEC,JVEC,OP1,NOP2,OP2)
      CALL MKWWOPE(IVEC,JVEC,OP0,OP1)
      CALL MKWWOPF(IVEC,JVEC,NOP2,OP2)
      CALL MKWWOPG(IVEC,JVEC,OP1)
      CALL MKWWOPH(IVEC,JVEC,OP0)

      RETURN
      END
      SUBROUTINE MKWWOPA(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      DIMENSION OP1(NASHT,NASHT),OP2(NOP2),OP3(NOP3)

C Given the coefficients for two excitation operators of the
C type VJTU = Case A, available in vectors numbered IVEC and
C JVEC on file, construct the zero-, one-, two-, and three-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.
C Formula used:
C  W1(tuv,i)(conj)*W2(xyz,j) = dij * (  -Evuxtyz -dyu Evzxt
C                     - dyt Evuxz - dxu Evtyz - dxu dyt Evz
C                     + 2 dtx Evuyz + 2 dtx dyu Evz )
* ------------------------------------------------------------
* PAM 2008: Sectioning over non-active superindices added
* at Krapperup Labour Camp, May 2008. Some comments of changes
* only at this routine; similar changes in MKWWOPB--MKWWOPH.
* ------------------------------------------------------------

      ICASE=1
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
* PAM2008: Added sectioning over non-active superindex
* but this will obviously hardly affect this case.
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for this block of excitation amplitudes:
* Sectioning sizes instead. Replaced code:
*        CALL GETMEM('WWW1','ALLO','REAL',LW1,NW)
*        CALL GETMEM('WWW2','ALLO','REAL',LW2,NW)
* replace with:
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
C Pick up a symmetry block of W1 and W2
*        CALL RDBLKC(ISYM,ICASE,IVEC,WORK(LW1))
*        CALL RDBLKC(ISYM,ICASE,JVEC,WORK(LW2))
C Allocate space for the contraction:
        NWSCT=MIN(NAS,1000)
        NWPROD=NWSCT**2
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
* End of addition
C Loop over sections of WW1 and WW2:
        DO ITUVSTA=1,NAS,NWSCT
          LW1A=LW1-1+ITUVSTA
          ITUVEND=MIN(ITUVSTA-1+NWSCT,NAS)
          MWS1=ITUVEND+1-ITUVSTA
          DO IXYZSTA=1,NAS,NWSCT
            LW2A=LW2-1+IXYZSTA
            IXYZEND=MIN(IXYZSTA-1+NWSCT,NAS)
            MWS2=IXYZEND+1-IXYZSTA
C Multiply WProd = (W1 sect )*(W2 sect transpose)
*            CALL DGEMM_('N','T',
*     &                  MWS1,MWS2,NIS,
*     &                  1.0d0,WORK(LW1A),NAS,
*     &                  WORK(LW2A),NAS,
*     &                  0.0d0,WORK(LWPROD),NWSCT)
* Replaced, due to sectioning over inactives:
            CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
            CALL DGEMM_('N','T',
     &                  MWS1,MWS2,NCOL,
     &                  1.0d0,WORK(LW1A),NAS,
     &                  WORK(LW2A),NAS,
     &                  1.0d0,WORK(LWPROD),NWSCT)
* End of replacement

C Loop over (TUV) in its section
          DO ITUV=ITUVSTA,ITUVEND
            IW1=ITUV+1-ITUVSTA
            ITUVABS=ITUV+NTUVES(ISYM)
            ITABS=MTUV(1,ITUVABS)
            IUABS=MTUV(2,ITUVABS)
            IVABS=MTUV(3,ITUVABS)
            IVU=IVABS+NASHT*(IUABS-1)
C Loop over (XYZ) in its section
          DO IXYZ=IXYZSTA,IXYZEND
            IW2=IXYZ+1-IXYZSTA
            IXYZABS=IXYZ+NTUVES(ISYM)
            IXABS=MTUV(1,IXYZABS)
            IYABS=MTUV(2,IXYZABS)
            IZABS=MTUV(3,IXYZABS)
            IXT=IXABS+NASHT*(ITABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IWPROD=IW1+NWSCT*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
C Remember:
C  W1(tuv,i)(conj)*W2(xyz,j) = dij * (  -Evuxtyz -dyu Evzxt
C                     - dyt Evuxz - dxu Evtyz - dxu dyt Evz
C                     + 2 dtx Evuyz + 2 dtx dyu Evz )
C Contrib to 3-particle operator:
            IF(IVU.LT.IXT) THEN
              IF(IVU.GE.IYZ) THEN
                JVU=IXT
                JXT=IVU
                JYZ=IYZ
              ELSE IF(IXT.LT.IYZ) THEN
                  JVU=IYZ
                  JXT=IXT
                  JYZ=IVU
              ELSE
                  JVU=IXT
                  JXT=IYZ
                  JYZ=IVU
              END IF
            ELSE
              IF(IVU.LT.IYZ) THEN
                JVU=IYZ
                JXT=IVU
                JYZ=IXT
              ELSE IF (IXT.GE.IYZ) THEN
                JVU=IVU
                JXT=IXT
                JYZ=IYZ
              ELSE
                JVU=IVU
                JXT=IYZ
                JYZ=IXT
              END IF
            END IF
            JVUXTYZ=((JVU+1)*JVU*(JVU-1))/6+(JXT*(JXT-1))/2+JYZ
            OP3(JVUXTYZ)=OP3(JVUXTYZ)-WPROD
C Contrib to 2-particle operator, from -dyu Evzxt:
            IF(IYABS.EQ.IUABS) THEN
              IVZ=IVABS+NASHT*(IZABS-1)
              IXT=IXABS+NASHT*(ITABS-1)
              IF(IVZ.GE.IXT) THEN
                JVZXT=(IVZ*(IVZ-1))/2+IXT
              ELSE
                JVZXT=(IXT*(IXT-1))/2+IVZ
              END IF
              OP2(JVZXT)=OP2(JVZXT)-WPROD
            END IF
C Contrib to 2-particle operator, from -dyt Evuxz:
            IF(IYABS.EQ.ITABS) THEN
              IVU=IVABS+NASHT*(IUABS-1)
              IXZ=IXABS+NASHT*(IZABS-1)
              IF(IVU.GE.IXZ) THEN
                JVUXZ=(IVU*(IVU-1))/2+IXZ
              ELSE
                JVUXZ=(IXZ*(IXZ-1))/2+IVU
              END IF
              OP2(JVUXZ)=OP2(JVUXZ)-WPROD
C Contrib to 1-particle operator, from -dxu dyt Evz:
              IF(IXABS.EQ.IUABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)-WPROD
              END IF
            END IF
C Contrib to 2-particle operator, from -dxu Evtyz:
            IF(IXABS.EQ.IUABS) THEN
              IVT=IVABS+NASHT*(ITABS-1)
              IYZ=IYABS+NASHT*(IZABS-1)
              IF(IVT.GE.IYZ) THEN
                JVTYZ=(IVT*(IVT-1))/2+IYZ
              ELSE
                JVTYZ=(IYZ*(IYZ-1))/2+IVT
              END IF
              OP2(JVTYZ)=OP2(JVTYZ)-WPROD
            END IF
C Contrib to 2-particle operator, from +2 dtx Evuyz:
            IF(ITABS.EQ.IXABS) THEN
              IVU=IVABS+NASHT*(IUABS-1)
              IYZ=IYABS+NASHT*(IZABS-1)
              IF(IVU.GE.IYZ) THEN
                JVUYZ=(IVU*(IVU-1))/2+IYZ
              ELSE
                JVUYZ=(IYZ*(IYZ-1))/2+IVU
              END IF
              OP2(JVUYZ)=OP2(JVUYZ)+2.0D0*WPROD
C Contrib to 1-particle operator, from +2 dtx dyu Evz:
              IF(IYABS.EQ.IUABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)+2.0D0*WPROD
              END IF
            END IF
           END DO
          END DO
         END DO
        END DO
* PAM2008, an added sectioning loop
        END DO
C Deallocate temporary space:
*        CALL GETMEM('WWW1','FREE','REAL',LW1,NW)
*        CALL GETMEM('WWW2','FREE','REAL',LW2,NW)
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPB(IVEC,JVEC,OP0,OP1,NOP2,OP2)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      DIMENSION OP1(NASHT,NASHT),OP2(NOP2)

C Given the coefficients for two excitation operators, available in
C vectors numbered IVEC and JVEC on file, use the blocks for
C excitation cases VJTI(+) and VJTI(-), i.e. cases 2 and 3, to
C construct the zero-, one-, and two-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.
C Formulae used:
C For the B+ case (i.e. case 2)
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu + 2 Eytxu -2dxt Eyu
C           -2dyu Ext -2dyt Exu -2dxu Eyt + 4 dxt dyu + 4 dxu dyt)
C For the B- case (i.e. case 3)
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu - 2 Eytxu -6dxt Eyu
C           -6dyu Ext +6dyt Exu +6dxu Eyt +12 dxt dyu -12 dxu dyt)

C FIRST THE B+ i.e. VJTI+ i.e. CASE 2 -----------------------------
      ICASE=2
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 888
C Allocate space for one section of excitation amplitudes:
C Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        IF(IVEC.EQ.JVEC) THEN
          LW2=LW1
        ELSE
          CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        END IF
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Loop over sections:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA+MDVEC-1,MDVEC)
         NCOL=IIEND-IISTA+1
         NSCT=NAS*NCOL
        CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
        IF (IVEC.NE.JVEC) CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1 sect )*(W2 sect transpose)
        CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
        END DO
C Deallocate W1 and W2
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        IF(IVEC.NE.JVEC) CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (TU)
        DO ITU=1,NAS
          IW1=ITU
          ITUABS=ITU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITUABS)
          IUABS=MTGEU(2,ITUABS)
C Loop over (XY)
          DO IXY=1,NAS
            IW2=IXY
            IXYABS=IXY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXYABS)
            IYABS=MTGEU(2,IXYABS)
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IYT=IYABS+NASHT*(ITABS-1)
            IXU=IXABS+NASHT*(IUABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
C Remember:
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu + 2 Eytxu -2dxt Eyu
C           -2dyu Ext -2dyt Exu -2dxu Eyt + 4 dxt dyu + 4 dxu dyt)
C Contrib to 2-particle operator, from 2 Extyu:
            IF(IXT.GE.IYU) THEN
              JXTYU=(IXT*(IXT-1))/2+IYU
            ELSE
              JXTYU=(IYU*(IYU-1))/2+IXT
            END IF
            OP2(JXTYU)=OP2(JXTYU)+2.0D0*WPROD
C Contrib to 1-particle operator, from -2dxt Eyu
            IF(IXABS.EQ.ITABS) THEN
              OP1(IYABS,IUABS)=OP1(IYABS,IUABS)-2.0D0*WPROD
            END IF
C Contrib to 1-particle operator, from -2dyu Ext
            IF(IYABS.EQ.IUABS) THEN
              OP1(IXABS,ITABS)=OP1(IXABS,ITABS)-2.0D0*WPROD
C Contrib to 0-particle operator, from +4 dxt dyu
              IF(IXABS.EQ.ITABS) OP0=OP0 + 4.0D0*WPROD
            END IF
C Contrib to 2-particle operator, from 2 Eytxu:
            IF(IYT.GT.IXU) THEN
              JYTXU=(IYT*(IYT-1))/2+IXU
            ELSE
              JYTXU=(IXU*(IXU-1))/2+IYT
            END IF
            OP2(JYTXU)=OP2(JYTXU)+2.0D0*WPROD
C Contrib to 1-particle operator, from -2dyt Exu
            IF(IYABS.EQ.ITABS) THEN
              OP1(IXABS,IUABS)=OP1(IXABS,IUABS)-2.0D0*WPROD
            END IF
C Contrib to 1-particle operator, from -2dxu Eyt
            IF(IXABS.EQ.IUABS) THEN
              OP1(IYABS,ITABS)=OP1(IYABS,ITABS)-2.0D0*WPROD
C Contrib to 0-particle operator, from +4 dyt dxu
              IF(IYABS.EQ.ITABS) OP0=OP0 + 4.0D0*WPROD
            END IF
          END DO
        END DO
C Deallocate matrix product:
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 888    CONTINUE
      END DO
C Then THE B- i.e. VJTI- i.e. CASE 3 -----------------------------
      ICASE=3
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1 sect )*(W2 sect transpose)
        CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
        END DO
C Deallocate W1, W2
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (TU)
          DO ITU=1,NAS
            IW1=ITU
            ITUABS=ITU+NTGTUES(ISYM)
            ITABS=MTGTU(1,ITUABS)
            IUABS=MTGTU(2,ITUABS)
C Loop over (XY)
          DO IXY=1,NAS
            IW2=IXY
            IXYABS=IXY+NTGTUES(ISYM)
            IXABS=MTGTU(1,IXYABS)
            IYABS=MTGTU(2,IXYABS)
            IXT=IXABS+NASHT*(ITABS-1)
            IYU=IYABS+NASHT*(IUABS-1)
            IYT=IYABS+NASHT*(ITABS-1)
            IXU=IXABS+NASHT*(IUABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
C Remember:
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu - 2 Eytxu -6dxt Eyu
C           -6dyu Ext +6dyt Exu +6dxu Eyt +12 dxt dyu -12 dxu dyt)
C Contrib to 2-particle operator, from 2 Extyu:
            IF(IXT.GE.IYU) THEN
              JXTYU=(IXT*(IXT-1))/2+IYU
            ELSE
              JXTYU=(IYU*(IYU-1))/2+IXT
            END IF
            OP2(JXTYU)=OP2(JXTYU)+2.0D0*WPROD
C Contrib to 1-particle operator, from -6dxt Eyu
            IF(IXABS.EQ.ITABS) THEN
              OP1(IYABS,IUABS)=OP1(IYABS,IUABS)-6.0D0*WPROD
            END IF
C Contrib to 1-particle operator, from -6dyu Ext
            IF(IYABS.EQ.IUABS) THEN
              OP1(IXABS,ITABS)=OP1(IXABS,ITABS)-6.0D0*WPROD
C Contrib to 0-particle operator, from +12 dxt dyu
              IF(IXABS.EQ.ITABS) OP0=OP0 + 12.0D0*WPROD
            END IF
C Contrib to 2-particle operator, from -2 Eytxu:
            IF(IYT.GE.IXU) THEN
              JYTXU=(IYT*(IYT-1))/2+IXU
            ELSE
              JYTXU=(IXU*(IXU-1))/2+IYT
            END IF
            OP2(JYTXU)=OP2(JYTXU)-2.0D0*WPROD
C Contrib to 1-particle operator, from +6dyt Exu
            IF(IYABS.EQ.ITABS) THEN
              OP1(IXABS,IUABS)=OP1(IXABS,IUABS)+6.0D0*WPROD
            END IF
C Contrib to 1-particle operator, from +6dxu Eyt
            IF(IXABS.EQ.IUABS) THEN
              OP1(IYABS,ITABS)=OP1(IYABS,ITABS)+6.0D0*WPROD
C Contrib to 0-particle operator, from -12 dyt dxu
              IF(IYABS.EQ.ITABS) OP0=OP0 -12.0D0*WPROD
            END IF
          END DO
        END DO
C Deallocate matrix product
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      DIMENSION OP1(NASHT,NASHT),OP2(NOP2),OP3(NOP3)

C Given the coefficients for two excitation operators of the
C type ATVX = Case C, available in vectors numbered IVEC and
C JVEC on file, construct the zero-, one-, two-, and three-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.
C Formula used:
C  W1(tuv,a)(conj)*W2(xyz,b) = dab * ( Evutxyz +dyu Evztx
C                       + dyx Evutz + dtu Evxyz + dtu dyx Evz )

      ICASE=4
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        NWSCT=MIN(NAS,1000)
        NWPROD=NWSCT**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Loop over sections of WW1 and WW2:
        DO ITUVSTA=1,NAS,NWSCT
          LW1A=LW1-1+ITUVSTA
          ITUVEND=MIN(ITUVSTA-1+NWSCT,NAS)
          MWS1=ITUVEND+1-ITUVSTA
          DO IXYZSTA=1,NAS,NWSCT
            IXYZEND=MIN(IXYZSTA-1+NWSCT,NAS)
            LW2A=LW2-1+IXYZSTA
            MWS2=IXYZEND+1-IXYZSTA
C Multiply WProd = (W1 sect )*(W2 sect transpose)
            CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
            CALL DGEMM_('N','T',
     &                  MWS1,MWS2,NCOL,
     &                  1.0d0,WORK(LW1A),NAS,
     &                  WORK(LW2A),NAS,
     &                  1.0d0,WORK(LWPROD),NWSCT)

C Loop over (TUV) in its section
          DO ITUV=ITUVSTA,ITUVEND
            IW1=ITUV+1-ITUVSTA
            ITUVABS=ITUV+NTUVES(ISYM)
            ITABS=MTUV(1,ITUVABS)
            IUABS=MTUV(2,ITUVABS)
            IVABS=MTUV(3,ITUVABS)
            IVU=IVABS+NASHT*(IUABS-1)
C Loop over (XYZ) in its section
          DO IXYZ=IXYZSTA,IXYZEND
            IW2=IXYZ+1-IXYZSTA
            IXYZABS=IXYZ+NTUVES(ISYM)
            IXABS=MTUV(1,IXYZABS)
            IYABS=MTUV(2,IXYZABS)
            IZABS=MTUV(3,IXYZABS)
            ITX=ITABS+NASHT*(IXABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IWPROD=IW1+NWSCT*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
C Remember:
C  W1(tuv,a)(conj)*W2(xyz,b) = dab * ( Evutxyz +dyu Evztx
C                       + dyx Evutz + dtu Evxyz + dtu dyx Evz )
C Contrib to 3-particle operator:
            IF(IVU.LT.ITX) THEN
              IF(IVU.GE.IYZ) THEN
                JVU=ITX
                JTX=IVU
                JYZ=IYZ
              ELSE IF(ITX.LT.IYZ) THEN
                  JVU=IYZ
                  JTX=ITX
                  JYZ=IVU
              ELSE
                  JVU=ITX
                  JTX=IYZ
                  JYZ=IVU
              END IF
            ELSE
              IF(IVU.LT.IYZ) THEN
                JVU=IYZ
                JTX=IVU
                JYZ=ITX
              ELSE IF (ITX.GE.IYZ) THEN
                JVU=IVU
                JTX=ITX
                JYZ=IYZ
              ELSE
                JVU=IVU
                JTX=IYZ
                JYZ=ITX
              END IF
            END IF
            JVUTXYZ=((JVU+1)*JVU*(JVU-1))/6+(JTX*(JTX-1))/2+JYZ
            OP3(JVUTXYZ)=OP3(JVUTXYZ)+WORK(LWPROD-1+IWPROD)
C Contrib to 2-particle operator, from  dyu Evztx:
            IF(IYABS.EQ.IUABS) THEN
              IVZ=IVABS+NASHT*(IZABS-1)
              ITX=ITABS+NASHT*(IXABS-1)
              IF(IVZ.GE.ITX) THEN
                JVZTX=(IVZ*(IVZ-1))/2+ITX
              ELSE
                JVZTX=(ITX*(ITX-1))/2+IVZ
              END IF
              OP2(JVZTX)=OP2(JVZTX)+WPROD
            END IF
C Contrib to 2-particle operator, from  dyx Evutz:
            IF(IYABS.EQ.IXABS) THEN
              IVU=IVABS+NASHT*(IUABS-1)
              ITZ=ITABS+NASHT*(IZABS-1)
              IF(IVU.GE.ITZ) THEN
                JVUTZ=(IVU*(IVU-1))/2+ITZ
              ELSE
                JVUTZ=(ITZ*(ITZ-1))/2+IVU
              END IF
              OP2(JVUTZ)=OP2(JVUTZ)+WPROD
            END IF
C Contrib to 2-particle operator, from  dtu Evxyz:
            IF(ITABS.EQ.IUABS) THEN
              IVX=IVABS+NASHT*(IXABS-1)
              IYZ=IYABS+NASHT*(IZABS-1)
              IF(IVX.GE.IYZ) THEN
                JVXYZ=(IVX*(IVX-1))/2+IYZ
              ELSE
                JVXYZ=(IYZ*(IYZ-1))/2+IVX
              END IF
              OP2(JVXYZ)=OP2(JVXYZ)+WPROD
C Contrib to 1-particle operator, from  dtu dyx Evz:
              IF(IYABS.EQ.IXABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)+WPROD
              END IF
            END IF
           END DO
          END DO
         END DO
        END DO
* Extra sectioning loop added...
        END DO
C Deallocate temporary space:
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPD(IVEC,JVEC,OP1,NOP2,OP2)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      DIMENSION OP1(NASHT,NASHT),OP2(NOP2)

C Given the coefficients for two excitation operators, available in
C vectors numbered IVEC and JVEC on file, use the blocks for
C excitation case AIVX, i.e. case 5, to
C construct the zero-, one-, and two-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.
C Formulae used:
C  (W1A(tu,ai) conj)*(W2A(tu,ai)) = 2*(Eutxy + dtx Euy)
C  (W1A(tu,ai) conj)*(W2B(tu,ai)) =  -(Eutxy + dtx Euy)
C  (W1B(tu,ai) conj)*(W2A(tu,ai)) =  -(Eutxy + dtx Euy)
C  (W1B(tu,ai) conj)*(W2B(tu,ai)) =  -Extuy + 2dtx Euy

      ICASE=5
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NAS1=NAS/2
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
        END DO
C Deallocate space for this block of excitation amplitudes:
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (TU)
          DO JTU=1,NAS1
            IW1A=JTU
            IW1B=JTU+NAS1
            JTUABS=JTU+NTUES(ISYM)
            ITABS=MTU(1,JTUABS)
            IUABS=MTU(2,JTUABS)
C Loop over (XY)
          DO JXY=1,NAS1
            IW2A=JXY
            IW2B=JXY+NAS1
            JXYABS=JXY+NTUES(ISYM)
            IXABS=MTU(1,JXYABS)
            IYABS=MTU(2,JXYABS)
            IUT=IUABS+NASHT*(ITABS-1)
            IXY=IXABS+NASHT*(IYABS-1)
            IXT=IXABS+NASHT*(ITABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            IWPRAA=IW1A+NAS*(IW2A-1)
            IWPRAB=IW1A+NAS*(IW2B-1)
            IWPRBA=IW1B+NAS*(IW2A-1)
            IWPRBB=IW1B+NAS*(IW2B-1)
            WPRAA=WORK(LWPROD-1+IWPRAA)
            WPRAB=WORK(LWPROD-1+IWPRAB)
            WPRBA=WORK(LWPROD-1+IWPRBA)
            WPRBB=WORK(LWPROD-1+IWPRBB)
C Remember:
C (W1A(tu,ai) conj)*(W2A(tu,ai)) = 2*(Eutxy + dtx Euy)
C (W1A(tu,ai) conj)*(W2B(tu,ai)) =  -(Eutxy + dtx Euy)
C (W1B(tu,ai) conj)*(W2A(tu,ai)) =  -(Eutxy + dtx Euy)
C (W1B(tu,ai) conj)*(W2B(tu,ai)) =  -Extuy + 2dtx Euy
C Contrib to 2-particle operator, from Eutxy:
            IF(IUT.GE.IXY) THEN
              JUTXY=(IUT*(IUT-1))/2+IXY
            ELSE
              JUTXY=(IXY*(IXY-1))/2+IUT
            END IF
            OP2(JUTXY)=OP2(JUTXY)+(2.0D0*WPRAA-WPRAB-WPRBA)
C Contrib to 1-particle operator, from Euy:
            IF(ITABS.EQ.IXABS) THEN
              OP1(IUABS,IYABS)= OP1(IUABS,IYABS)
     &               +(2.0D0*WPRAA-WPRAB-WPRBA+2.0D0*WPRBB)
            END IF
C Contrib to 2-particle operator, from Extuy:
            IF(IXT.GE.IUY) THEN
              JXTUY=(IXT*(IXT-1))/2+IUY
            ELSE
              JXTUY=(IUY*(IUY-1))/2+IXT
            END IF
            OP2(JXTUY)=OP2(JXTUY)-WPRBB
          END DO
        END DO
C Deallocate matrix product:
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPE(IVEC,JVEC,OP0,OP1)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Presently symmetry blocking is disregarded.
      DIMENSION OP1(NASHT,NASHT)

C Given the coefficients for two excitation operators, available in
C vectors numbered IVEC and JVEC on file, use the blocks for
C excitation case VJAI, i.e. cases 6 and 7, to express the
C product (Op in IVEC conjugated)(Op in JVEC) as a one-body
C operator on the CASSCF space.
C Formula used:
C  (W1(t,aij) conj)*(W2(x,bkl)) = dik*djl*dab*(2*dtx - Etx)
C the same for both cases 6 and 7.

C Loop over cases
      DO ICASE=6,7
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
         END DO
C Deallocate space for this block of excitation amplitudes:
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (T)
          DO IT=1,NAS
            IW1=IT
            ITABS=IT+NAES(ISYM)
C Loop over (X)
          DO IX=1,NAS
            IW2=IX
            IXABS=IX+NAES(ISYM)
            IWPROD=IW1+NAS*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
            OP1(ITABS,IXABS)=OP1(ITABS,IXABS)-WPROD
            IF(ITABS.EQ.IXABS) OP0=OP0+2.0D0*WPROD
          END DO
        END DO
C Deallocate matrix product
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
C End of loop over cases.
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPF(IVEC,JVEC,NOP2,OP2)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      DIMENSION OP2(NOP2)
C Given the coefficients for two excitation operators, available in
C vectors numbered IVEC and JVEC on file, use the blocks for
C excitation cases BVAT(+) and BVAT(-), i.e. cases 8 and 9, to
C construct the zero-, one-, and two-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.
C Formulae used:
C For the F+ case (i.e. case 8)
C W1(tu,ab)(conj)*W2(xy,cd) = (dac*dbd)*(2 Etxuy + 2 Etyux)
C For the F- case (i.e. case 9)
C W1(tu,ab)(conj)*W2(xy,cd) = (dac*dbd)*(2 Etxuy - 2 Etyux)

C FIRST THE F+ i.e. BVAT+ i.e. CASE 8 -----------------------------
      ICASE=8
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 888
C Allocate space for one section of excitation amplitudes:
C Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        IF(JVEC.EQ.IVEC) THEN
          LW2=LW1
        ELSE
          CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        END IF
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1 sect )*(W2 sect transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
         END DO
C Deallocate W1 and W2
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        IF(JVEC.NE.IVEC) CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (TU)
        DO ITU=1,NAS
          IW1=ITU
          ITUABS=ITU+NTGEUES(ISYM)
          ITABS=MTGEU(1,ITUABS)
          IUABS=MTGEU(2,ITUABS)
C Loop over (XY)
          DO IXY=1,NAS
            IW2=IXY
            IXYABS=IXY+NTGEUES(ISYM)
            IXABS=MTGEU(1,IXYABS)
            IYABS=MTGEU(2,IXYABS)
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            ITY=ITABS+NASHT*(IYABS-1)
            IUX=IUABS+NASHT*(IXABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
C Remember: C For the F+ case (i.e. case 8)
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy + 2 Etyux)
C Contrib to 2-particle operator, from 2 Etxuy:
            IF(ITX.GE.IUY) THEN
              JTXUY=(ITX*(ITX-1))/2+IUY
            ELSE
              JTXUY=(IUY*(IUY-1))/2+ITX
            END IF
            OP2(JTXUY)=OP2(JTXUY)+2.0D0*WPROD
C Contrib to 2-particle operator, from 2 Etyux:
            IF(ITY.GE.IUX) THEN
              JTYUX=(ITY*(ITY-1))/2+IUX
            ELSE
              JTYUX=(IUX*(IUX-1))/2+ITY
            END IF
            OP2(JTYUX)=OP2(JTYUX)+2.0D0*WPROD
          END DO
        END DO
C Deallocate matrix product:
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
888     CONTINUE
      END DO
C THEN THE F- i.e. BVAT- i.e. CASE 9 -----------------------------
      ICASE=9
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
C Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        IF(JVEC.EQ.IVEC) THEN
          LW2=LW1
        ELSE
          CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        END IF
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1 sect )*(W2 sect transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
        END DO
C Deallocate W1 and W2
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        IF(JVEC.NE.IVEC) CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (TU)
        DO ITU=1,NAS
          IW1=ITU
          ITUABS=ITU+NTGTUES(ISYM)
          ITABS=MTGTU(1,ITUABS)
          IUABS=MTGTU(2,ITUABS)
C Loop over (XY)
          DO IXY=1,NAS
            IW2=IXY
            IXYABS=IXY+NTGTUES(ISYM)
            IXABS=MTGTU(1,IXYABS)
            IYABS=MTGTU(2,IXYABS)
            ITX=ITABS+NASHT*(IXABS-1)
            IUY=IUABS+NASHT*(IYABS-1)
            ITY=ITABS+NASHT*(IYABS-1)
            IUX=IUABS+NASHT*(IXABS-1)
            IWPROD=IW1+NAS*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
C Remember: C For the F- case (i.e. case 9)
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy - 2 Etyux)
C Contrib to 2-particle operator, from 2 Etxuy:
            IF(ITX.GE.IUY) THEN
              JTXUY=(ITX*(ITX-1))/2+IUY
            ELSE
              JTXUY=(IUY*(IUY-1))/2+ITX
            END IF
            OP2(JTXUY)=OP2(JTXUY)+2.0D0*WPROD
C Contrib to 2-particle operator, from -2 Etyux:
            IF(ITY.GE.IUX) THEN
              JTYUX=(ITY*(ITY-1))/2+IUX
            ELSE
              JTYUX=(IUX*(IUX-1))/2+ITY
            END IF
            OP2(JTYUX)=OP2(JTYUX)-2.0D0*WPROD
          END DO
        END DO
C Deallocate matrix product:
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPG(IVEC,JVEC,OP1)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Presently symmetry blocking is disregarded.
      DIMENSION OP1(NASHT,NASHT)

C Given the coefficients for two excitation operators, available in
C vectors numbered IVEC and JVEC on file, use the blocks for
C excitation case BJAT, i.e. cases 10 and 11, to express the
C product (Op in IVEC conjugated)(Op in JVEC) as a one-body
C operator on the CASSCF space.
C Formula used:
C  (W1(t,aij) conj)*(W2(x,bkl)) = dik*djl*dab* Etx
C the same for both cases 10 and 11.

C Loop over cases
      DO ICASE=10,11
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL GETMEM('WWPROD','ALLO','REAL',LWPROD,NWPROD)
        CALL DCOPY_(NWPROD,0.0D0,0,WORK(LWPROD),1)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              1.0d0,WORK(LW1),NAS,
     &              WORK(LW2),NAS,
     &              1.0d0,WORK(LWPROD),NAS)
        END DO
C Deallocate space for this block of excitation amplitudes:
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)

C Loop over (T)
          DO IT=1,NAS
            IW1=IT
            ITABS=IT+NAES(ISYM)
C Loop over (X)
          DO IX=1,NAS
            IW2=IX
            IXABS=IX+NAES(ISYM)
            IWPROD=IW1+NAS*(IW2-1)
            WPROD=WORK(LWPROD-1+IWPROD)
            OP1(ITABS,IXABS)=OP1(ITABS,IXABS)+WPROD
          END DO
        END DO
C Deallocate matrix product
        CALL GETMEM('WWPROD','FREE','REAL',LWPROD,NWPROD)
 999    CONTINUE
      END DO
C End of loop over cases.
      END DO
      RETURN
      END
      SUBROUTINE MKWWOPH(IVEC,JVEC,OP0)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"

C Given the coefficients for two excitation operators, available in
C vectors numbered IVEC and JVEC on file, use the blocks for
C excitation case BJAI, i.e. cases 12 and 13, to express the
C product (Op in IVEC conjugated)(Op in JVEC) as a zero-body
C operator, i.e. a scalar factor, in the CASSCF space.
C Formula used:
C  (W1(ij,ab) conj)*(W2(kl,cd)) = dik*djl*dac*dbd
C the same for both cases 10 and 11.

C Loop over cases
      DO ICASE=12,13
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) GOTO 999
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL GETMEM('WWW1','ALLO','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','ALLO','REAL',LW2,NAS*MDVEC)
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         NSCT=NAS*NCOL
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,WORK(LW1))
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,WORK(LW2))
C Pick up a symmetry block of W1 and W2
         OP0=OP0+DDOT_(NSCT,WORK(LW1),1,WORK(LW2),1)
        END DO
        CALL GETMEM('WWW1','FREE','REAL',LW1,NAS*MDVEC)
        CALL GETMEM('WWW2','FREE','REAL',LW2,NAS*MDVEC)
 999    CONTINUE
      END DO
C End of loop over cases.
      END DO
      RETURN
      END
