!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE MKWWOPA(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      USE SUPERINDEX, only: MTUV
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,       &
     &                         NTUVES
      use constants, only: Zero, One, Two
      use definitions, only: iwp, wp
      IMPLICIT None

! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
! NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2, NOP3
      real(kind=wp), Intent(inout) :: OP1(NASHT,NASHT),OP2(NOP2),       &
     &                                OP3(NOP3)

      real(kind=wp), Allocatable:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS, ITUV,   &
     &                  ITUVABS, ITUVEND, ITUVSTA, IUABS, IVABS, IVT,   &
     &                  IVU, IVZ, IW1, IW2, IWPROD, IXABS, IXT, IXYZ,   &
     &                  IXYZABS, IXYZEND, IXYZSTA, IXZ, IYABS, IYZ,     &
     &                  IZABS, JVTYZ, JVU, JVUXTYZ, JVUXZ, JVUYZ, JVZXT,&
     &                  JXT, JYZ, LW1A, LW2A, MDVEC, MWS1, MWS2, NAS,   &
     &                  NCOL, NIS, NWPROD, NWSCT
      real(kind=wp) W_PROD
! Given the coefficients for two excitation operators of the
! type VJTU = Case A, available in vectors numbered IVEC and
! JVEC on file, construct the zero-, one-, two-, and three-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formula used:
!  W1(tuv,i)(conj)*W2(xyz,j) = dij * (  -Evuxtyz -dyu Evzxt
!                     - dyt Evuxz - dxu Evtyz - dxu dyt Evz
!                     + 2 dtx Evuyz + 2 dtx dyu Evz )
! ------------------------------------------------------------
! PAM 2008: Sectioning over non-active superindices added
! at Krapperup Labour Camp, May 2008. Some comments of changes
! only at this routine; similar changes in MKWWOPB--MKWWOPH.
! ------------------------------------------------------------

      ICASE=1
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
! PAM2008: Added sectioning over non-active superindex
! but this will obviously hardly affect this case.
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
!        NW=NAS*NIS
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for this block of excitation amplitudes:
! Sectioning sizes instead. Replaced code:
!        CALL mma_allocate(W1,NW,Label='W1')
!        CALL mma_allocate(W2,NW,Label='W2')
! replace with:
! Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,LABEL='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
! Pick up a symmetry block of W1 and W2
!        CALL RDBLKC(ISYM,ICASE,IVEC,W1,NAS*MDVEC)
!        CALL RDBLKC(ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! Allocate space for the contraction:
        NWSCT=MIN(NAS,1000)
        NWPROD=NWSCT**2
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
! Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! End of addition
! Loop over sections of WW1 and WW2:
        DO ITUVSTA=1,NAS,NWSCT
          LW1A=ITUVSTA
          ITUVEND=MIN(ITUVSTA-1+NWSCT,NAS)
          MWS1=ITUVEND+1-ITUVSTA
          DO IXYZSTA=1,NAS,NWSCT
            LW2A=IXYZSTA
            IXYZEND=MIN(IXYZSTA-1+NWSCT,NAS)
            MWS2=IXYZEND+1-IXYZSTA
! Multiply WProd = (W1 sect )*(W2 sect transpose)
!            CALL DGEMM_('N','T',
!     &                  MWS1,MWS2,NIS,
!     &                  One,W1(LW1A),NAS,
!     &                  W2(LW2A),NAS,
!     &                  Zero,WPROD,NWSCT)
! Replaced, due to sectioning over inactives:
            WPROD(:)=Zero
            CALL DGEMM_('N','T',                                        &
     &                  MWS1,MWS2,NCOL,                                 &
     &                  One,W1(LW1A),NAS,                               &
     &                  W2(LW2A),NAS,                                   &
     &                  One,WPROD,NWSCT)
! End of replacement

! Loop over (TUV) in its section
          DO ITUV=ITUVSTA,ITUVEND
            IW1=ITUV+1-ITUVSTA
            ITUVABS=ITUV+NTUVES(ISYM)
            ITABS=MTUV(1,ITUVABS)
            IUABS=MTUV(2,ITUVABS)
            IVABS=MTUV(3,ITUVABS)
            IVU=IVABS+NASHT*(IUABS-1)
! Loop over (XYZ) in its section
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
! Remember:
!  W1(tuv,i)(conj)*W2(xyz,j) = dij * (  -Evuxtyz -dyu Evzxt
!                     - dyt Evuxz - dxu Evtyz - dxu dyt Evz
!                     + 2 dtx Evuyz + 2 dtx dyu Evz )
! Contrib to 3-particle operator:
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
! Contrib to 2-particle operator, from -dyu Evzxt:
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
! Contrib to 2-particle operator, from -dyt Evuxz:
            IF(IYABS.EQ.ITABS) THEN
              IVU=IVABS+NASHT*(IUABS-1)
              IXZ=IXABS+NASHT*(IZABS-1)
              IF(IVU.GE.IXZ) THEN
                JVUXZ=(IVU*(IVU-1))/2+IXZ
              ELSE
                JVUXZ=(IXZ*(IXZ-1))/2+IVU
              END IF
              OP2(JVUXZ)=OP2(JVUXZ)-W_PROD
! Contrib to 1-particle operator, from -dxu dyt Evz:
              IF(IXABS.EQ.IUABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)-W_PROD
              END IF
            END IF
! Contrib to 2-particle operator, from -dxu Evtyz:
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
! Contrib to 2-particle operator, from +2 dtx Evuyz:
            IF(ITABS.EQ.IXABS) THEN
              IVU=IVABS+NASHT*(IUABS-1)
              IYZ=IYABS+NASHT*(IZABS-1)
              IF(IVU.GE.IYZ) THEN
                JVUYZ=(IVU*(IVU-1))/2+IYZ
              ELSE
                JVUYZ=(IYZ*(IYZ-1))/2+IVU
              END IF
              OP2(JVUYZ)=OP2(JVUYZ)+Two*W_PROD
! Contrib to 1-particle operator, from +2 dtx dyu Evz:
              IF(IYABS.EQ.IUABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)+Two*W_PROD
              END IF
            END IF
           END DO
          END DO
         END DO
        END DO
! PAM2008, an added sectioning loop
        END DO
! Deallocate temporary space:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)
        CALL mma_deallocate(WPROD)
      END DO

      END SUBROUTINE MKWWOPA
