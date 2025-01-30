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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************
      SUBROUTINE DIATERM2_GAS( FACTOR,  ITASK,    VEC, NBLOCK, IBLOCK,  &
     &                           IOFF,  JPERT,    J12,    JDC)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas, only: NSTSO
      use lucia_data, only: ECORE_ORIG,ECORE
      use lucia_data, only: IPRDIA
      use lucia_data, only: MXNSTR
      use lucia_data, only: NTOOB,IREOST,IREOTS,NACOB
#ifdef _DEBUGPRINT_
      use lucia_data, only: IDC,IPERTOP,IBSPGPFTP
#endif
      use lucia_data, only: NOCTYP
      use lucia_data, only: NELEC
      use csm_data, only: NSMST
! = DIATERM_GAS, just J12 added !
!
! Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
! Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
!
! For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
!
! If JPERT.NE.0, the perturbation operator as defined by IPART is used.
!
! Jeppe Olsen, August 1995
!
      IMPLICIT NONE
!
      REAL*8 FACTOR
      INTEGER ITASK,IBLOCK(8,*)
      REAL*8 VEC(*)
      INTEGER IOFF,JPERT,J12,JDC

      Integer, Allocatable:: LASTR(:), LBSTR(:)
      Real*8, Allocatable:: LSCR2(:)
      Real*8, Allocatable:: LJ(:), LK(:), LXB(:), LH1D(:), LRJKA(:)
      INTEGER, EXTERNAL:: IMNMX
      INTEGER NTEST,IATP,IBTP,NAEL,NBEL,NOCTPA,MAXA,NBLOCK
      REAL*8 ECOREP,SHIFT,FACTORX
#ifdef _DEBUGPRINT_
      INTEGER NOCTPB,IOCTPA,IOCTPB
#endif
!
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
!
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
!
!
!     IF(JPERT.EQ.0) THEN
!. Use full Hamiltonian
!       I12 = 2
!       IPERTOP = 0
!     ELSE
!. Use perturbation operator
!       IF(IPART.EQ.1) THEN
!. Moller-Plesset partitioning
!         I12 = 1
!         IPERTOP = 1
!       ELSE IF(IPART.EQ.2) THEN
!. Epstein-Nesbet Partitioning
!         I12 = 2
!         IPERTOP = 0
!       END IF
!     END IF

#ifdef _DEBUGPRINT_
!. Offsets for alpha and beta supergroups
      NOCTPB = NOCTYP(IBTP)
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   DIATERM2_GAS speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,J12,IPERTOP',JPERT,J12,IPERTOP
      END IF
#else
      Call Unused_INteger(JPERT)
#endif
!. A bit of scracth
      CALL mma_allocate(LJ   ,NTOOB**2,Label='LJ')
      CALL mma_allocate(LK   ,NTOOB**2,Label='LK')
      Call mma_allocate(LSCR2,2*NTOOB**2,Label='LSCR2')
      CALL mma_allocate(LXB  ,NACOB,Label='LX')
      CALL mma_allocate(LH1D ,NACOB,Label='LH1D')
!. Space for blocks of strings
      Call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
      Call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
      MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
      CALL mma_allocate(LRJKA,MAXA,Label='LRJKA')
!. Diagonal of one-body integrals and coulomb and exchange integrals
!. Integrals assumed in place so :
      CALL GT1DIA(LH1D)
      IF(J12.EQ.2)                                                      &
     &CALL GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
!. Core energy not included
      ECOREP = 0.0D0
      SHIFT = ECORE_ORIG-ECORE
      FACTORX = FACTOR + SHIFT
      CALL DIATERMS_GAS(NAEL,LASTR,NBEL,LBSTR,                          &
     &                  NACOB,VEC,NSMST,                                &
     &                  LH1D,JDC,LXB,LJ,LK,                             &
     &                  NSTSO(IATP)%I,NSTSO(IBTP)%I,                    &
     &                  ECOREP,0,0,IPRDIA,NTOOB,LRJKA,J12,              &
     &                  IBLOCK(1,IOFF),NBLOCK,ITASK, FACTORX,0,[0])
!
!    &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
!.Flush local memory
      CALL mma_deallocate(LJ)
      CALL mma_deallocate(LK)
      Call mma_deallocate(LSCR2)
      CALL mma_deallocate(LXB)
      CALL mma_deallocate(LH1D)
      Call mma_deallocate(LASTR)
      Call mma_deallocate(LBSTR)
      CALL mma_deallocate(LRJKA)
!
#ifdef _DEBUGPRINT_
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTS(      VEC,IBLOCK(1,IOFF),NBLOCK,  NSMST,            &
     &                 NSTSO(IATP)%I,NSTSO(IBTP)%I,IDC)
      END IF
#endif
!
      END SUBROUTINE DIATERM2_GAS
