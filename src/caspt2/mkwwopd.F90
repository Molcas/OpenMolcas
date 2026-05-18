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
      SUBROUTINE MKWWOPD(IVEC,JVEC,OP1,NOP2,OP2)
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      USE SUPERINDEX, only: MTU
      use EQSOLV, only: MODVEC
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP,       &
     &                         NTUES
      IMPLICIT None

! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
      integer(kind=iwp), intent(in):: IVEC, JVEC, NOP2
      real(kind=wp), intent(inout):: OP1(NASHT,NASHT),OP2(NOP2)

      real(kind=wp), ALLOCATABLE:: W1(:), W2(:), WPROD(:)
      integer(kind=iwp) ICASE, IIEND, IISTA, ISCT, ISYM, ITABS,         &
     &                  IUABS, IXABS, IYABS,                            &
     &                  MDVEC, NAS, NCOL, NIS, NWPROD,                  &
     &                  IUT, IUY, IW1A, IW1B, IW2A, IW2B, IWPRAA,       &
     &                  IWPRAB, IWPRBA, IWPRBB, IXT, IXY, JTU, JTUABS,  &
     &                  JUTXY, JXTUY, JXY, JXYABS, NAS1
      real(kind=wp) WPRAA, WPRAB, WPRBA, WPRBB
! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation case AIVX, i.e. case 5, to
! construct the zero-, one-, and two-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formulae used:
!  (W1A(tu,ai) conj)*(W2A(tu,ai)) = 2*(Eutxy + dtx Euy)
!  (W1A(tu,ai) conj)*(W2B(tu,ai)) =  -(Eutxy + dtx Euy)
!  (W1B(tu,ai) conj)*(W2A(tu,ai)) =  -(Eutxy + dtx Euy)
!  (W1B(tu,ai) conj)*(W2B(tu,ai)) =  -Extuy + 2dtx Euy

      ICASE=5
! Loop over symmetry ISYM
      DO ISYM=1,NSYM
        NAS=NASUP(ISYM,ICASE)
        NAS1=NAS/2
        NIS=NISUP(ISYM,ICASE)
        IF(NINDEP(ISYM,ICASE).EQ.0) CYCLE
! Allocate space for one section of excitation amplitudes:
        MDVEC=MODVEC(ISYM,ICASE)
        CALL mma_allocate(W1,NAS*MDVEC,Label='W1')
        CALL mma_allocate(W2,NAS*MDVEC,Label='W2')
        NWPROD=NAS**2
! Allocate space for the contraction:
        CALL mma_allocate(WPROD,NWPROD,Label='WPROD')
        WPROD(:)=Zero
! Sectioning loop added:
        ISCT=0
        DO IISTA=1,NIS,MDVEC
         ISCT=ISCT+1
         IIEND=MIN(IISTA-1+MDVEC,NIS)
         NCOL=1+IIEND-IISTA
         CALL RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
         CALL RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
! Multiply WProd = (W1)*(W2 transpose)
         CALL DGEMM_('N','T',                                           &
     &              NAS,NAS,NCOL,                                       &
     &              One,W1,NAS,                                         &
     &              W2,NAS,                                             &
     &              One,WPROD,NAS)
        END DO
! Deallocate space for this block of excitation amplitudes:
        CALL mma_deallocate(W1)
        CALL mma_deallocate(W2)

! Loop over (TU)
          DO JTU=1,NAS1
            IW1A=JTU
            IW1B=JTU+NAS1
            JTUABS=JTU+NTUES(ISYM)
            ITABS=MTU(1,JTUABS)
            IUABS=MTU(2,JTUABS)
! Loop over (XY)
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
! Remember:
! (W1A(tu,ai) conj)*(W2A(tu,ai)) = 2*(Eutxy + dtx Euy)
! (W1A(tu,ai) conj)*(W2B(tu,ai)) =  -(Eutxy + dtx Euy)
! (W1B(tu,ai) conj)*(W2A(tu,ai)) =  -(Eutxy + dtx Euy)
! (W1B(tu,ai) conj)*(W2B(tu,ai)) =  -Extuy + 2dtx Euy
! Contrib to 2-particle operator, from Eutxy:
            IF(IUT.GE.IXY) THEN
              JUTXY=(IUT*(IUT-1))/2+IXY
            ELSE
              JUTXY=(IXY*(IXY-1))/2+IUT
            END IF
            OP2(JUTXY)=OP2(JUTXY)+(Two*WPRAA-WPRAB-WPRBA)
! Contrib to 1-particle operator, from Euy:
            IF(ITABS.EQ.IXABS) THEN
              OP1(IUABS,IYABS)= OP1(IUABS,IYABS)                        &
     &               +(Two*WPRAA-WPRAB-WPRBA+Two*WPRBB)
            END IF
! Contrib to 2-particle operator, from Extuy:
            IF(IXT.GE.IUY) THEN
              JXTUY=(IXT*(IXT-1))/2+IUY
            ELSE
              JXTUY=(IUY*(IUY-1))/2+IXT
            END IF
            OP2(JXTUY)=OP2(JXTUY)-WPRBB
          END DO
        END DO
! Deallocate matrix product:
        CALL mma_deallocate(WPROD)
      END DO
      END SUBROUTINE MKWWOPD
