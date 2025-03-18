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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************
      SUBROUTINE CNTOST(ICONF,ICTSDT,NAEL,NBEL,                         &
     &                 IPRODT,IREFSM,                                   &
     &                 NORB,NEL,                                        &
     &                 IGENSG,ISGNA,ISGNB,ICNSTR,IAGRP,IBGRP,IOOS,      &
     &                 PSSIGN,IPRNT)

!
! Obtain pointer abs(ICTSDT(I)) giving address of determinant I in
! STRING ordering for determinant I in CSF ordering.
! Going between the two formats can involve a sign change . this is
! stored in the sign of ICTSDT)
! SGNCTS is thus to be multiplied with vector ordered in CSF ordering.
!
! December 1990 : NCNFCN,ICNFOK added
! January 1991  : IGENSG,ISGNA,ISGNB added
! April   1991  : LUCIA version
! September 1993 > Sign and address stored together
!
! ICNSTR .ne. 0 indicates that additional constraints on configurations
! should be checked  (IS = 0 )
! by calling CICNCH.ICNFOK(ICNF) is 1 of tests are passed, ICNFOK(ICNF)
! is zero if test fails
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: NTYP,NCNATS,NDPCNT,MINOP
      IMPLICIT NONE
      INTEGER ICONF(*),ICTSDT(*)
      INTEGER NAEL,NBEL
      INTEGER IPRODT(*)
      INTEGER IREFSM,NORB,NEL,IGENSG
      INTEGER ISGNA(*),ISGNB(*)
      INTEGER ICNSTR,IAGRP,IBGRP
      INTEGER IOOS(*)
      REAL*8 PSSIGN
      INTEGER IPRNT
!
! IWORK should at least be of length (MXDT+2)*NEL,
! where MXDT is the largest number of prototype determinants occuring
! in a single block.
!
!./SPINFO/
      Integer, Allocatable:: LDTBL(:), LIA(:), LIB(:), SCR23(:)
      INTEGER NTEST,MXDT,ITYP,ICNF,JDTABS,IPSFAC,ISGNAB,ICNBS0,IPBAS,   &
     &        IJKL_NUM,IDET,IOPEN,ICL,IOCC,IC,ICNBS,JDET,ISIGN,IABNUM
!
       NEL = NAEL + NBEL
       NTEST=0000

!.. Local memory
!
! Largest number of dets for a given type
      MXDT = 0
      DO 10 ITYP = 1, NTYP
        MXDT   = MAX(MXDT,NDPCNT(ITYP) )
   10 CONTINUE
      CALL mma_allocate(LDTBL,MXDT*NEL,Label='LDTBL')
      CALL mma_allocate(LIA,NAEL,Label='LIA')
      CALL mma_allocate(LIB,NBEL,Label='LIB')
      CALL mma_allocate(SCR23,NEL,Label='SCR23')
!
!.. Loop over configurations and generate determinants in compact form
!
      ICNF = 0
      JDTABS = 0
      IPSFAC=0 ! Removes a compiler error
      ISGNAB=0 ! Removes a compiler error
      ICNBS0=0 ! dummy initialize
      IPBAS =0 ! dummy initialize
      ijkl_num =0 !yma counter
      DO 1000 ITYP = 1, NTYP
        IDET = NDPCNT(ITYP)
        IOPEN = ITYP + MINOP - 1
        ICL = (NEL - IOPEN) / 2
        IOCC = IOPEN + ICL
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
        ELSE
          ICNBS0 = ICNBS0 + NCNATS(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
        END IF
! Base for prototype determinants
        IF( ITYP .EQ. 1 ) THEN
          IPBAS = 1
        ELSE
          IPBAS = IPBAS + NDPCNT(ITYP-1)*(IOPEN-1)
        END IF
! Determinants for this configuration
        DO 900  IC = 1, NCNATS(ITYP,IREFSM)
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
!. Check orbital occupancy with additional constraints
          IF( NTEST .GE. 10 ) WRITE(6,*) ' IC ICNF ICNBS',IC,ICNF,ICNBS
          CALL CNDET(ICONF(ICNBS),IPRODT(IPBAS),IDET,                   &
     &               NEL,IOCC,IOPEN,ICL, LDTBL,IPRNT)
! Separate determinants into strings and determine string number .
          DO 800 JDET = 1,IDET
!            write(117,"(1X,I8,1X,A,1X)",advance='no')ITYP,"ITYP"  ! yma
            JDTABS = JDTABS + 1
            CALL DETSTR_MCLR(LDTBL(1+(JDET-1)*NEL),LIA,                 &
     &             LIB,NEL,NAEL,NBEL,NORB,ISIGN,SCR23,IPRNT)
            ijkl_num=ijkl_num+1
! Find number (and sign)of this determinant in string ordering
            ICTSDT(JDTABS) =IABNUM(LIA,LIB,IAGRP,IBGRP,IGENSG,          &
     &             ISGNA,ISGNB,ISGNAB,IOOS,NORB,IPSFAC,PSSIGN,          &
     &             IPRNT)
             IF(  DBLE(ISIGN*ISGNAB*IPSFAC) .eq. -1.0d0)then
               ICTSDT(JDTABS) = - ICTSDT(JDTABS)
             END IF
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
      CALL mma_deallocate(SCR23)
      CALL mma_deallocate(LIB)
      CALL mma_deallocate(LIA)
      CALL mma_deallocate(LDTBL)

!
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(ICNSTR)
      END SUBROUTINE CNTOST
