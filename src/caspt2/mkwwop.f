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
      use definitions, only: iwp, wp
      use constants, only: Zero
      use caspt2_module, only: NASHT
      IMPLICIT None

C Presently symmetry blocking is disregarded for OP2, OP3, but
C index pair C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2, NOP3
      real(kind=wp), intent(out):: OP0, OP1(NASHT,NASHT), OP2(NOP2),
     &                             OP3(NOP3)

C Given the coefficients for two excitation operators in the
C vectors numbered IVEC and C JVEC on file, construct the
C zero-, one-, two-, and three-body
C expansions of the product (Op in IVEC conjugated)(Op in JVEC)
C as operating on the CASSCF space.

      OP0=Zero
      OP1(:,:)=Zero
      OP2(:)=Zero
      OP3(:)=Zero
      CALL MKWWOPA(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      CALL MKWWOPB(IVEC,JVEC,OP0,OP1,NOP2,OP2)
      CALL MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      CALL MKWWOPD(IVEC,JVEC,OP1,NOP2,OP2)
      CALL MKWWOPE(IVEC,JVEC,OP0,OP1)
      CALL MKWWOPF(IVEC,JVEC,NOP2,OP2)
      CALL MKWWOPG(IVEC,JVEC,OP1)
      CALL MKWWOPH(IVEC,JVEC,OP0)

      END SUBROUTINE MKWWOP

      SUBROUTINE MKWWOPA(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      use definitions, only: iwp, wp
      USE SUPERINDEX, only: MTUV
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,
     &                         NTUVES
      IMPLICIT None

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2, NOP3
      real(kind=wp), Intent(inout) :: OP1(NASHT,NASHT),OP2(NOP2),
     &                                OP3(NOP3)

      real(kind=wp), Allocatable:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS, ITUV,
     &                  ITUVABS, ITUVEND, ITUVSTA, IUABS, IVABS, IVT,
     &                  IVU, IVZ, IW1, IW2, IWPROD, IXABS, IXT, IXYZ,
     &                  IXYZABS, IXYZEND, IXYZSTA, IXZ, IYABS, IYZ,
     &                  IZABS, JVTYZ, JVU, JVUXTYZ, JVUXZ, JVUYZ, JVZXT,
     &                  JXT, JYZ, LW1A, LW2A, MDVEC, MWS1, MWS2, NAS,
     &                  NCOL, NIS, NWPROD, NWSCT
      real(kind=wp) W_PROD
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
*        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for this block of excitation amplitudes:
* Sectioning sizes instead. Replaced code:
*        CALL mma_allocate(W1,NW,Label='W1')
*        CALL mma_allocate(W2,NW,Label='W2')
* replace with:
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,LABEL='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
C Pick up a symmetry block of W1 and W2
*        CALL RDBLKC(ISYM,ICASE,IVEC,W1)
*        CALL RDBLKC(ISYM,ICASE,JVEC,W2)
C Allocate space for the contraction:
        NWSCT=MIN(NAS,1000)
        NWPROD=NWSCT**2
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
* End of addition
C Loop over sections of WW1 and WW2:
        DO ITUVSTA=1,NAS,NWSCT
          LW1A=ITUVSTA
          ITUVEND=MIN(ITUVSTA-1+NWSCT,NAS)
          MWS1=ITUVEND+1-ITUVSTA
          DO IXYZSTA=1,NAS,NWSCT
            LW2A=IXYZSTA
            IXYZEND=MIN(IXYZSTA-1+NWSCT,NAS)
            MWS2=IXYZEND+1-IXYZSTA
C Multiply WProd = (W1 sect )*(W2 sect transpose)
*            CALL DGEMM_('N','T',
*     &                  MWS1,MWS2,NIS,
*     &                  1.0d0,W1(LW1A),NAS,
*     &                  W2(LW2A),NAS,
*     &                  0.0d0,WPROD,NWSCT)
* Replaced, due to sectioning over inactives:
            WPROD(:)=0.0D0
            CALL DGEMM_('N','T',
     &                  MWS1,MWS2,NCOL,
     &                  1.0d0,W1(LW1A),NAS,
     &                  W2(LW2A),NAS,
     &                  1.0d0,WPROD,NWSCT)
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
            W_PROD=WPROD(IWPROD)
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
            OP3(JVUXTYZ)=OP3(JVUXTYZ)-W_PROD
C Contrib to 2-particle operator, from -dyu Evzxt:
            IF(IYABS.EQ.IUABS) THEN
              IVZ=IVABS+NASHT*(IZABS-1)
              IXT=IXABS+NASHT*(ITABS-1)
              IF(IVZ.GE.IXT) THEN
                JVZXT=(IVZ*(IVZ-1))/2+IXT
              ELSE
                JVZXT=(IXT*(IXT-1))/2+IVZ
              END IF
              OP2(JVZXT)=OP2(JVZXT)-W_PROD
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
              OP2(JVUXZ)=OP2(JVUXZ)-W_PROD
C Contrib to 1-particle operator, from -dxu dyt Evz:
              IF(IXABS.EQ.IUABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)-W_PROD
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
              OP2(JVTYZ)=OP2(JVTYZ)-W_PROD
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
              OP2(JVUYZ)=OP2(JVUYZ)+2.0D0*W_PROD
C Contrib to 1-particle operator, from +2 dtx dyu Evz:
              IF(IYABS.EQ.IUABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)+2.0D0*W_PROD
              END IF
            END IF
           END DO
          END DO
         END DO
        END DO
* PAM2008, an added sectioning loop
        END DO
C Deallocate temporary space:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)
        CALL mma_deallocate(WPROD)
      END DO

      END SUBROUTINE MKWWOPA

      SUBROUTINE MKWWOPB(IVEC,JVEC,OP0,OP1,NOP2,OP2)
      use definitions, only: iwp, wp
      use Constants, only: Zero, One, Two, Four, Six, Twelve
      USE SUPERINDEX, only: MTGEU, MTGTU
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,
     &                         NTGEUES, NTGTUES
      IMPLICIT None

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      integer(kind=iwp), Intent(in):: IVEC, JVEC, NOP2
      real(kind=wp), intent(inout):: OP0, OP1(NASHT,NASHT),OP2(NOP2)

      real(kind=wp), Allocatable, TARGET:: W1(:), W2_H(:), WPROD(:)
      real(kind=wp), POINTER::  W2(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS, ITU,
     &                  ITUABS, IUABS, IW1, IW2, IWPROD, IXABS, IXT,
     &                  IXU, IXY, IXYABS, IYABS, IYT, IYU, JXTYU,
     &                  JYTXU, MDVEC, NAS, NCOL, NIS, NWPROD
      real(kind=wp) W_PROD
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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
C Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        IF(IVEC.EQ.JVEC) THEN
          W2=>W1
        ELSE
          CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
          W2=>W2_H
        END IF
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Loop over sections:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA+MDVEC-1,MDVEC)
         NCOL=IIEND-IISTA+1
        CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
        IF (IVEC.NE.JVEC) CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1 sect )*(W2 sect transpose)
        CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
        END DO
C Deallocate W1 and W2
        CALL mma_deallocate(W1)
        IF(IVEC.NE.JVEC) CALL mma_deallocate(W2_H)
        W2=>Null()

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
            W_PROD=WPROD(IWPROD)
C Remember:
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu + 2 Eytxu -2dxt Eyu
C           -2dyu Ext -2dyt Exu -2dxu Eyt + 4 dxt dyu + 4 dxu dyt)
C Contrib to 2-particle operator, from 2 Extyu:
            IF(IXT.GE.IYU) THEN
              JXTYU=(IXT*(IXT-1))/2+IYU
            ELSE
              JXTYU=(IYU*(IYU-1))/2+IXT
            END IF
            OP2(JXTYU)=OP2(JXTYU)+Two*W_PROD
C Contrib to 1-particle operator, from -2dxt Eyu
            IF(IXABS.EQ.ITABS) THEN
              OP1(IYABS,IUABS)=OP1(IYABS,IUABS)-Two*W_PROD
            END IF
C Contrib to 1-particle operator, from -2dyu Ext
            IF(IYABS.EQ.IUABS) THEN
              OP1(IXABS,ITABS)=OP1(IXABS,ITABS)-Two*W_PROD
C Contrib to 0-particle operator, from +4 dxt dyu
              IF(IXABS.EQ.ITABS) OP0=OP0 + Four*W_PROD
            END IF
C Contrib to 2-particle operator, from 2 Eytxu:
            IF(IYT.GT.IXU) THEN
              JYTXU=(IYT*(IYT-1))/2+IXU
            ELSE
              JYTXU=(IXU*(IXU-1))/2+IYT
            END IF
            OP2(JYTXU)=OP2(JYTXU)+Two*W_PROD
C Contrib to 1-particle operator, from -2dyt Exu
            IF(IYABS.EQ.ITABS) THEN
              OP1(IXABS,IUABS)=OP1(IXABS,IUABS)-Two*W_PROD
            END IF
C Contrib to 1-particle operator, from -2dxu Eyt
            IF(IXABS.EQ.IUABS) THEN
              OP1(IYABS,ITABS)=OP1(IYABS,ITABS)-Two*W_PROD
C Contrib to 0-particle operator, from +4 dyt dxu
              IF(IYABS.EQ.ITABS) OP0=OP0 + Four*W_PROD
            END IF
          END DO
        END DO
C Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
C Then THE B- i.e. VJTI- i.e. CASE 3 -----------------------------
      ICASE=3
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
        W2=>W2_H
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1 sect )*(W2 sect transpose)
        CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
        END DO
C Deallocate W1, W2
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2_H)
        W2=>Null()

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
            W_PROD=WPROD(IWPROD)
C Remember:
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Extyu - 2 Eytxu -6dxt Eyu
C           -6dyu Ext +6dyt Exu +6dxu Eyt +12 dxt dyu -12 dxu dyt)
C Contrib to 2-particle operator, from 2 Extyu:
            IF(IXT.GE.IYU) THEN
              JXTYU=(IXT*(IXT-1))/2+IYU
            ELSE
              JXTYU=(IYU*(IYU-1))/2+IXT
            END IF
            OP2(JXTYU)=OP2(JXTYU)+Two*W_PROD
C Contrib to 1-particle operator, from -6dxt Eyu
            IF(IXABS.EQ.ITABS) THEN
              OP1(IYABS,IUABS)=OP1(IYABS,IUABS)-Six*W_PROD
            END IF
C Contrib to 1-particle operator, from -6dyu Ext
            IF(IYABS.EQ.IUABS) THEN
              OP1(IXABS,ITABS)=OP1(IXABS,ITABS)-Six*W_PROD
C Contrib to 0-particle operator, from +12 dxt dyu
              IF(IXABS.EQ.ITABS) OP0=OP0 + Twelve*W_PROD
            END IF
C Contrib to 2-particle operator, from -2 Eytxu:
            IF(IYT.GE.IXU) THEN
              JYTXU=(IYT*(IYT-1))/2+IXU
            ELSE
              JYTXU=(IXU*(IXU-1))/2+IYT
            END IF
            OP2(JYTXU)=OP2(JYTXU)-Two*W_PROD
C Contrib to 1-particle operator, from +6dyt Exu
            IF(IYABS.EQ.ITABS) THEN
              OP1(IXABS,IUABS)=OP1(IXABS,IUABS)+Six*W_PROD
            END IF
C Contrib to 1-particle operator, from +6dxu Eyt
            IF(IXABS.EQ.IUABS) THEN
              OP1(IYABS,ITABS)=OP1(IYABS,ITABS)+Six*W_PROD
C Contrib to 0-particle operator, from -12 dyt dxu
              IF(IYABS.EQ.ITABS) OP0=OP0 -Twelve*W_PROD
            END IF
          END DO
        END DO
C Deallocate matrix product
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPB

      SUBROUTINE MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      use definitions, only: iwp, wp
      USE SUPERINDEX, only: MTUV
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP, NTUVES
      IMPLICIT None

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      integer(kind=iwp), intent(in)::  IVEC, JVEC, NOP2, NOP3
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT),OP2(NOP2),
     &                               OP3(NOP3)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS,
     &                  IUABS, IW1, IW2, IWPROD, IXABS, IYABS,
     &                  MDVEC, NAS, NCOL, NIS, NWPROD,
     &                  ITUV, ITUVABS, ITUVEND, ITUVSTA, ITX, ITZ,
     &                  IVABS, IVU, IVX, IVZ, IXYZ, IXYZABS, IXYZEND,
     &                  IXYZSTA, IYZ, IZABS, JTX, JVU, JVUTXYZ, JVUTZ,
     &                  JVXYZ, JVZTX, JYZ, LW1A, LW2A, MWS1, MWS2, NWSCT
      real(kind=wp) W_PROD
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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
        NWSCT=MIN(NAS,1000)
        NWPROD=NWSCT**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Loop over sections of WW1 and WW2:
        DO ITUVSTA=1,NAS,NWSCT
          LW1A=ITUVSTA
          ITUVEND=MIN(ITUVSTA-1+NWSCT,NAS)
          MWS1=ITUVEND+1-ITUVSTA
          DO IXYZSTA=1,NAS,NWSCT
            IXYZEND=MIN(IXYZSTA-1+NWSCT,NAS)
            LW2A=IXYZSTA
            MWS2=IXYZEND+1-IXYZSTA
C Multiply WProd = (W1 sect )*(W2 sect transpose)
            WPROD(:)=0.0D0
            CALL DGEMM_('N','T',
     &                  MWS1,MWS2,NCOL,
     &                  1.0d0,W1(LW1A),NAS,
     &                  W2(LW2A),NAS,
     &                  1.0d0,WPROD,NWSCT)

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
            W_PROD=WPROD(IWPROD)
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
            OP3(JVUTXYZ)=OP3(JVUTXYZ)+WPROD(IWPROD)
C Contrib to 2-particle operator, from  dyu Evztx:
            IF(IYABS.EQ.IUABS) THEN
              IVZ=IVABS+NASHT*(IZABS-1)
              ITX=ITABS+NASHT*(IXABS-1)
              IF(IVZ.GE.ITX) THEN
                JVZTX=(IVZ*(IVZ-1))/2+ITX
              ELSE
                JVZTX=(ITX*(ITX-1))/2+IVZ
              END IF
              OP2(JVZTX)=OP2(JVZTX)+W_PROD
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
              OP2(JVUTZ)=OP2(JVUTZ)+W_PROD
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
              OP2(JVXYZ)=OP2(JVXYZ)+W_PROD
C Contrib to 1-particle operator, from  dtu dyx Evz:
              IF(IYABS.EQ.IXABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)+W_PROD
              END IF
            END IF
           END DO
          END DO
         END DO
        END DO
* Extra sectioning loop added...
        END DO
C Deallocate temporary space:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPC

      SUBROUTINE MKWWOPD(IVEC,JVEC,OP1,NOP2,OP2)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      USE SUPERINDEX, only: MTU
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,
     &                         NTUES
      IMPLICIT None

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT),OP2(NOP2)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS,
     &                  IUABS, IXABS, IYABS,
     &                  MDVEC, NAS, NCOL, NIS, NWPROD,
     &                  IUT, IUY, IW1A, IW1B, IW2A, IW2B, IWPRAA,
     &                  IWPRAB, IWPRBA, IWPRBB, IXT, IXY, JTU, JTUABS,
     &                  JUTXY, JXTUY, JXY, JXYABS, NAS1
      real(kind=wp) WPRAA, WPRAB, WPRBA, WPRBB
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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
        END DO
C Deallocate space for this block of excitation amplitudes:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)

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
            WPRAA=WPROD(IWPRAA)
            WPRAB=WPROD(IWPRAB)
            WPRBA=WPROD(IWPRBA)
            WPRBB=WPROD(IWPRBB)
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
            OP2(JUTXY)=OP2(JUTXY)+(Two*WPRAA-WPRAB-WPRBA)
C Contrib to 1-particle operator, from Euy:
            IF(ITABS.EQ.IXABS) THEN
              OP1(IUABS,IYABS)= OP1(IUABS,IYABS)
     &               +(Two*WPRAA-WPRAB-WPRBA+Two*WPRBB)
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
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPD

      SUBROUTINE MKWWOPE(IVEC,JVEC,OP0,OP1)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,
     &                         NAES
      IMPLICIT None

C Presently symmetry blocking is disregarded.
      integer(kind=iwp), intent(in):: IVEC, JVEC
      real(kind=wp), intent(inout):: OP0, OP1(NASHT,NASHT)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, NWPROD, ISCT,
     &                  IISTA, IIEND, NCOL, IT, IW1, ITABS, IX, IW2,
     &                  IWPROD, IXABS
      real(kind=wp) W_PROD

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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
         END DO
C Deallocate space for this block of excitation amplitudes:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)

C Loop over (T)
          DO IT=1,NAS
            IW1=IT
            ITABS=IT+NAES(ISYM)
C Loop over (X)
          DO IX=1,NAS
            IW2=IX
            IXABS=IX+NAES(ISYM)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
            OP1(ITABS,IXABS)=OP1(ITABS,IXABS)-W_PROD
            IF(ITABS.EQ.IXABS) OP0=OP0+Two*W_PROD
          END DO
        END DO
C Deallocate matrix product
        CALL mma_deallocate(WPROD)
      END DO
C End of loop over cases.
      END DO
      END SUBROUTINE MKWWOPE

      SUBROUTINE MKWWOPF(IVEC,JVEC,NOP2,OP2)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      USE SUPERINDEX, only: MTGEU, MTGTU
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,
     &                         NTGEUES, NTGTUES
      IMPLICIT None

C Presently symmetry blocking is disregarded, but index pair
C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2
      real(kind=wp), intent(inout):: OP2(NOP2)

      real(kind=wp), ALLOCATABLE, TARGET:: W1(:), W2_H(:)
      real(kind=wp), ALLOCATABLE:: WPROD(:)
      real(kind=wp), POINTER:: W2(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, NWPROD, ISCT,
     &                  IISTA, IIEND, NCOL, IW1, ITABS, IW2,
     &                  IWPROD, IXABS, ITU, ITUABS, ITX, ITY, IUABS,
     &                  IUX, IUY, IXY, IXYABS, IYABS, JTXUY, JTYUX
      real(kind=wp) W_PROD

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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
C Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        IF(JVEC.EQ.IVEC) THEN
          W2=>W1
        ELSE
          CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
          W2=>W2_H
        END IF
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1 sect )*(W2 sect transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
         END DO
C Deallocate W1 and W2
        CALL mma_deallocate(W1)
        IF(JVEC.NE.IVEC) CALL mma_deallocate(W2_H)
        W2=>Null()

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
            W_PROD=WPROD(IWPROD)
C Remember: C For the F+ case (i.e. case 8)
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy + 2 Etyux)
C Contrib to 2-particle operator, from 2 Etxuy:
            IF(ITX.GE.IUY) THEN
              JTXUY=(ITX*(ITX-1))/2+IUY
            ELSE
              JTXUY=(IUY*(IUY-1))/2+ITX
            END IF
            OP2(JTXUY)=OP2(JTXUY)+Two*W_PROD
C Contrib to 2-particle operator, from 2 Etyux:
            IF(ITY.GE.IUX) THEN
              JTYUX=(ITY*(ITY-1))/2+IUX
            ELSE
              JTYUX=(IUX*(IUX-1))/2+ITY
            END IF
            OP2(JTYUX)=OP2(JTYUX)+Two*W_PROD
          END DO
        END DO
C Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
C THEN THE F- i.e. BVAT- i.e. CASE 9 -----------------------------
      ICASE=9
C Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
C Pick up a symmetry block of W1 and W2
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        IF(JVEC.EQ.IVEC) THEN
          W2=>W1
        ELSE
          CALL mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
          W2=>W2_H
        END IF
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1 sect )*(W2 sect transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
        END DO
C Deallocate W1 and W2
        CALL mma_deallocate(W1)
        IF(JVEC.NE.IVEC) CALL mma_deallocate(W2_H)
        W2=>Null()

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
            W_PROD=WPROD(IWPROD)
C Remember: C For the F- case (i.e. case 9)
C W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy - 2 Etyux)
C Contrib to 2-particle operator, from 2 Etxuy:
            IF(ITX.GE.IUY) THEN
              JTXUY=(ITX*(ITX-1))/2+IUY
            ELSE
              JTXUY=(IUY*(IUY-1))/2+ITX
            END IF
            OP2(JTXUY)=OP2(JTXUY)+Two*W_PROD
C Contrib to 2-particle operator, from -2 Etyux:
            IF(ITY.GE.IUX) THEN
              JTYUX=(ITY*(ITY-1))/2+IUX
            ELSE
              JTYUX=(IUX*(IUX-1))/2+ITY
            END IF
            OP2(JTYUX)=OP2(JTYUX)-Two*W_PROD
          END DO
        END DO
C Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPF

      SUBROUTINE MKWWOPG(IVEC,JVEC,OP1)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NASUP,NISUP,NINDEP,NASHT,NAES
      IMPLICIT None

C Presently symmetry blocking is disregarded.
      integer(kind=iwp), intent(in):: IVEC, JVEC
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, NWPROD, ISCT,
     &                  IISTA, IIEND, NCOL, IW1, ITABS, IW2, IT, IX,
     &                  IWPROD, IXABS
      real(kind=wp) W_PROD

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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,LABEL='W1')
        CALL mma_allocate(W2,NAS*MDVEC,LABEL='W2')
        NWPROD=NAS**2
C Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',
     &              NAS,NAS,NCOL,
     &              One,W1,NAS,
     &              W2,NAS,
     &              One,WPROD,NAS)
        END DO
C Deallocate space for this block of excitation amplitudes:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)

C Loop over (T)
          DO IT=1,NAS
            IW1=IT
            ITABS=IT+NAES(ISYM)
C Loop over (X)
          DO IX=1,NAS
            IW2=IX
            IXABS=IX+NAES(ISYM)
            IWPROD=IW1+NAS*(IW2-1)
            W_PROD=WPROD(IWPROD)
            OP1(ITABS,IXABS)=OP1(ITABS,IXABS)+W_PROD
          END DO
        END DO
C Deallocate matrix product
        CALL mma_deallocate(WPROD)
      END DO
C End of loop over cases.
      END DO
      END SUBROUTINE MKWWOPG

      SUBROUTINE MKWWOPH(IVEC,JVEC,OP0)
      use definitions, only: iwp, wp
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, NASUP, NISUP, NINDEP
      IMPLICIT None

      integer(kind=iwp), intent(in):: IVEC, JVEC
      real(kind=wp), intent(inout):: OP0

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:)
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, ISCT, IISTA,
     &                  IIEND, NCOL, NSCT
      real(kind=wp), external:: DDot_

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
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
C Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
* Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         NSCT=NAS*NCOL
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2)
C Pick up a symmetry block of W1 and W2
         OP0=OP0+DDOT_(NSCT,W1,1,W2,1)
        END DO
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)
      END DO
C End of loop over cases.
      END DO
      END SUBROUTINE MKWWOPH
