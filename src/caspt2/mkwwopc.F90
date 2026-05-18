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
      SUBROUTINE MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
      USE SUPERINDEX, only: MTUV
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP, NTUVES
      use constants, only: Zero, One
      use definitions, only: iwp, wp
      IMPLICIT None

! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
! NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)
      integer(kind=iwp), intent(in)::  IVEC, JVEC, NOP2, NOP3
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT),OP2(NOP2),        &
     &                               OP3(NOP3)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS,         &
     &                  IUABS, IW1, IW2, IWPROD, IXABS, IYABS,          &
     &                  MDVEC, NAS, NCOL, NIS, NWPROD,                  &
     &                  ITUV, ITUVABS, ITUVEND, ITUVSTA, ITX, ITZ,      &
     &                  IVABS, IVU, IVX, IVZ, IXYZ, IXYZABS, IXYZEND,   &
     &                  IXYZSTA, IYZ, IZABS, JTX, JVU, JVUTXYZ, JVUTZ,  &
     &                  JVXYZ, JVZTX, JYZ, LW1A, LW2A, MWS1, MWS2, NWSCT
      real(kind=wp) W_PROD
! Given the coefficients for two excitation operators of the
! type ATVX = Case C, available in vectors numbered IVEC and
! JVEC on file, construct the zero-, one-, two-, and three-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formula used:
!  W1(tuv,a)(conj)*W2(xyz,b) = dab * ( Evutxyz +dyu Evztx
!                       + dyx Evutz + dtu Evxyz + dtu dyx Evz )

      ICASE=4
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
        NWSCT=MIN(NAS,1000)
        NWPROD=NWSCT**2
! Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
! Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! Loop over sections of WW1 and WW2:
        DO ITUVSTA=1,NAS,NWSCT
          LW1A=ITUVSTA
          ITUVEND=MIN(ITUVSTA-1+NWSCT,NAS)
          MWS1=ITUVEND+1-ITUVSTA
          DO IXYZSTA=1,NAS,NWSCT
            IXYZEND=MIN(IXYZSTA-1+NWSCT,NAS)
            LW2A=IXYZSTA
            MWS2=IXYZEND+1-IXYZSTA
! Multiply WProd = (W1 sect )*(W2 sect transpose)
            WPROD(:)=Zero
            CALL DGEMM_('N','T',                                        &
     &                  MWS1,MWS2,NCOL,                                 &
     &                  One,W1(LW1A),NAS,                               &
     &                  W2(LW2A),NAS,                                   &
     &                  One,WPROD,NWSCT)

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
            ITX=ITABS+NASHT*(IXABS-1)
            IYZ=IYABS+NASHT*(IZABS-1)
            IWPROD=IW1+NWSCT*(IW2-1)
            W_PROD=WPROD(IWPROD)
! Remember:
!  W1(tuv,a)(conj)*W2(xyz,b) = dab * ( Evutxyz +dyu Evztx
!                       + dyx Evutz + dtu Evxyz + dtu dyx Evz )
! Contrib to 3-particle operator:
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
! Contrib to 2-particle operator, from  dyu Evztx:
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
! Contrib to 2-particle operator, from  dyx Evutz:
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
! Contrib to 2-particle operator, from  dtu Evxyz:
            IF(ITABS.EQ.IUABS) THEN
              IVX=IVABS+NASHT*(IXABS-1)
              IYZ=IYABS+NASHT*(IZABS-1)
              IF(IVX.GE.IYZ) THEN
                JVXYZ=(IVX*(IVX-1))/2+IYZ
              ELSE
                JVXYZ=(IYZ*(IYZ-1))/2+IVX
              END IF
              OP2(JVXYZ)=OP2(JVXYZ)+W_PROD
! Contrib to 1-particle operator, from  dtu dyx Evz:
              IF(IYABS.EQ.IXABS) THEN
                OP1(IVABS,IZABS)=OP1(IVABS,IZABS)+W_PROD
              END IF
            END IF
           END DO
          END DO
         END DO
        END DO
! Extra sectioning loop added...
        END DO
! Deallocate temporary space:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPC
