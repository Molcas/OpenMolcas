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
      SUBROUTINE INTCSF(NACTOB,NACTEL,MULTP,MS2,NORB1,NORB2,NORB3,      &
     &                  NEL1MN,NEL3MX,                                  &
     &                  LLCSF,NCNSM,ICNSTR,PSSIGN,                      &
     &                  IPRNT,lconf,lldet)

!
! Initializing routine for CSF-DET expansions of internal space
!
! Set up common block /CSFDIM/
! This gives information about the number of dets,csf's for
! each symmetry
!
! find local memory requirements for CSF routines
! Largest local memory requirements in CNFORD,CSFDET_MCLR is returned in
! LLCSF
!
!             DFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
!             CFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
!             DTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
!             CNSM(I)%ICONF : SPACE FOR STORING  NCNSM
!                        CONFIGURATION EXPANSIONS
!
! If PSSIGN .ne. 0, spin combinations are used !!
      Use Str_Info, only: DFTP, CFTP, DTOC, CNSM
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: MULTSP,MS2P,MINOP,MAXOP,NTYP,NCPCNT,NDPCNT,  &
     &                       NCNASM,NCNATS,NCSASM,NDTASM
      use DetDim, only: MXPCTP,MXPCSM,MXCNSM
      IMPLICIT NONE
      INTEGER NACTOB,NACTEL,MULTP,MS2,NORB1,NORB2,NORB3,NEL1MN,NEL3MX,  &
     &        LLCSF,NCNSM,ICNSTR
      REAL*8 PSSIGN
      INTEGER IPRNT,lconf,lldet
!     local variables
      Integer, Allocatable:: IICL(:), IIOP(:), IIOC(:)
      Integer NTEST,IMSCMB,MULTS,NEL,IEL1,IEL2,IEL3,IOP1,IOP2,IOP3,IOP, &
     &        ITP,IOPEN,IAEL,IBEL,ITYPE,LIDT,LICS,LDTOC,MXPTBL,MXDT,    &
     &        LCSFDT,LCNFOR,LDET,ILCNF,ISYM,ILLCNF,LLCONF,ITYP,ICL,     &
     &        ICNSM,IBION,IWEYLF
!
! Last modification : Sept 20 : sign and address of dets goes together
!                      in CNSM(:)%ICTS
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
!
      IF(PSSIGN.NE.0.0D0) THEN
       IMSCMB = 1
      ELSE
       IMSCMB = 0
      END IF
!
!*. Define parameters in SPINFO
!
      MULTSP = MULTP
      MS2P = MS2
      MULTS = MULTSP

      NEL = NACTEL
! ================================
!. Allowed number of open orbitals
! ================================
      MINOP = ABS(MS2)
      MAXOP = 0
      DO 5 IEL1 = NEL1MN,2*NORB1
        DO 4 IEL3 = 0,NEL3MX
          IEL2 = NACTEL-IEL1-IEL3
          IF(IEL2 .LT. 0 ) GOTO 4
          IOP1 = MIN(NORB1,2*NORB1-IEL1,IEL1)
          IOP2 = MIN(NORB2,2*NORB2-IEL2,IEL2)
          IOP3 = MIN(NORB3,2*NORB3-IEL3,IEL3)
          IOP = IOP1 + IOP2 + IOP3
          MAXOP = MAX(MAXOP,IOP)
    4   CONTINUE
    5 CONTINUE
!?    WRITE(6,*) ' MAXOP with RAS constraints :' ,MAXOP
      NTYP = MAXOP-MINOP + 1
!
      IF( NTYP .GT. MXPCTP ) THEN
        WRITE(6,*) '  NUMBER OF CONFIGURATION TYPES TO LARGE '
        WRITE(6,*) '  CHANGE PARAMETER MXPCTP TO AT LEAST ',NTYP
        WRITE(6,*) '  CURRENT VALUE OF MXPCTP ',MXPCTP
        write(6,*) ' MTYP IN LUSPIN TO SMALL '
        call Abend()
      END IF
!
      IF( NTEST .GE. 5 )                                                &
     &WRITE(6,*) ' MINOP MAXOP NTYP ',MINOP,MAXOP,NTYP
! ================================================
!. Number of sd's and csf's per configuration type
! ================================================
      DO 10 ITP = 1,NTYP
        IOPEN = MINOP+ITP - 1
        IAEL = (IOPEN + MS2 ) / 2
        IBEL = (IOPEN - MS2 ) / 2
        IF(IAEL+IBEL .EQ. IOPEN ) THEN
          NDPCNT(ITP) = IBION(IOPEN,IAEL)
          IF(IMSCMB.NE.0.AND.IOPEN.NE.0)                                &
     &    NDPCNT(ITP) = NDPCNT(ITP)/2
          IF(IOPEN .GE. MULTS-1) THEN
            NCPCNT(ITP) = IWEYLF(IOPEN,MULTS)
          ELSE
            NCPCNT(ITP) = 0
          END IF
        ELSE
          NDPCNT(ITP) = 0
          NCPCNT(ITP) = 0
        END IF
   10 CONTINUE
      IF(NTEST.GE.2) THEN
      WRITE(6,'(/A)') ' Information about prototype configurations '
      WRITE(6,'( A)') ' ========================================== '
      WRITE(6,'(/A)')
      IF(IMSCMB.EQ.0) THEN
        WRITE(6,'(/A)') ' Combinations = Slater determinants'
      ELSE
        WRITE(6,'(/A)') ' Combinations = Spin combinations '
      END IF
      WRITE(6,'(/A)')                                                   &
     &'  Open orbitals   Combinations    CSFs '
      DO 580 IOPEN = MINOP,MAXOP,2
        ITYPE = IOPEN - MINOP + 1
        WRITE(6,'(5X,I3,10X,I6,7X,I6)')                                 &
     &  IOPEN,NDPCNT(ITYPE),NCPCNT(ITYPE)
  580 CONTINUE
      END IF
! =================================================
!*. Number of Combinations and CSF's per  symmetry
! =================================================
      CALL mma_allocate(IICL,NACTOB,Label='IICL')
      CALL mma_allocate(IIOP,NACTOB,Label='IIOP')
      CALL mma_allocate(IIOC,NORB1+NORB2+NORB3,Label='IIOC')

      CALL CISIZE(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,NACTEL,               &
     &            MINOP,MAXOP,MXPCTP,MXPCSM,NCNATS,NCNASM,NDTASM,       &
     &            NCSASM,                                               &
     &            NDPCNT,NCPCNT,                                        &
     &            IICL,IIOP,IIOC,IPRNT)

      CALL mma_deallocate(IIOC)
      CALL mma_deallocate(IIOP)
      CALL mma_deallocate(IICL)
! ==============================================
!   Permanent and local memory for csf routines
! ==============================================
!
!    memory for CSDTMT arrays.
!    Largest block of proto type determinants .
!    Largest number of prototype determinants
!    All configurations( of any specific symmetry )
!

      LIDT = 0
      LICS = 0
      LDTOC = 0
      MXPTBL = 0
      MXDT = 0
      LCONF = 0
      DO 11 ITP = 1,NTYP
        IOPEN = MINOP+ITP - 1
        LIDT = LIDT + NDPCNT(ITP) * IOPEN
        LICS = LICS + NCPCNT(ITP) * IOPEN
        LDTOC= LDTOC + NCPCNT(ITP)*NDPCNT(ITP)
        MXDT =   MAX(MXDT,NDPCNT(ITP) )
        MXPTBL = MAX(NDPCNT(ITP)*IOPEN,MXPTBL)
   11 CONTINUE
!. local memory for CSFDET_MCLR
      LCSFDT = MXPTBL + MAXOP
!. local memory for CNFORD
      LCNFOR = MAX(2*NTYP+NACTOB,(MXDT+2)*NACTEL)
!. local memory for any routine used in construction of csf basis
      LLCSF = MAX(LCSFDT,LCNFOR)
!. Memory needed to store ICONF array
      LCONF = 0
      LDET = 0
      ILCNF = 0
      DO 30 ISYM = 1, MXPCSM
        ILLCNF = 0
        LLCONF = 0
        LDET = MAX(LDET,NDTASM(ISYM))
        DO 25 ITYP = 1, NTYP
          IOPEN = ITYP+MINOP-1
          ICL = ( NEL-IOPEN)/2
          LLCONF = LLCONF + NCNATS(ITYP,ISYM)*(IOPEN+ICL)
          ILLCNF = ILLCNF + NCNATS(ITYP,ISYM)
   25   CONTINUE
!?      WRITE(6,*) ' MEMORY FOR HOLDING CONFS OF SYM... ',ISYM,LLCONF
        LCONF = MAX(LCONF,LLCONF)
        ILCNF = MAX(ILCNF,ILLCNF)
   30 CONTINUE

      ! notice the ILCNF number ! yma

       IF(NTEST.GE.5) THEN
       WRITE(6,'(/A,I8)')                                               &
     & '  Memory for holding largest list of configurations ',LCONF
       WRITE(6,'(/A,I8)')                                               &
     & '  Size of largest CI expansion (combinations)',LDET
       WRITE(6,'(/A,I8)')                                               &
     & '  Size of largest CI expansion (confs)',ILCNF
       END IF
       call xflush(6) !yma

!. permanent memory for csf proto type arrays

      Call mma_allocate(DFTP,LIDT,Label='DFTP')
      Call mma_allocate(CFTP,LICS,Label='CFTP')
      CALL mma_allocate(DTOC,LDTOC,Label='DTOC')

!. Permanent arrays for reordering and phases
      IF(NCNSM .GT. MXCNSM ) THEN
        WRITE(6,'(A,2I2)')                                              &
     &  '  TROUBLE IN CSFDIM NCNSM > MXCNSM : NCNSM,MXCNSM',            &
     &  NCNSM,MXCNSM
        write(6,*) ' CSFDIM : NCNSM  IS GREATER THAN MXCNSM '
        Call Abend()
      END IF
      DO 60 ICNSM = 1, NCNSM
        CALL mma_allocate(CNSM(ICNSM)%ICONF,LCONF,Label='ICONF')
        CALL mma_allocate(CNSM(ICNSM)%ICTS,LDET,Label='ICTS')
   60 CONTINUE
!
!
!
      lldet=ldet

! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(ICNSTR)
      END SUBROUTINE INTCSF
