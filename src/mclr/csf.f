************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE CNFORD(ICTSDT,ICONF,
     *           IREFSM,NORB,IPRODT,NCNFTP,
     *           NEL,ICNSTR,IGENSG,ISGNA,ISGNB,IAGRP,IBGRP,IOOS,
     *           NORB1,NORB2,NORB3,NEL1MN,NEL3MX,NAEL,NBEL,MINOP,MAXOP,
     *           PSSIGN,IPRNT)
*
* Generate configurations in ICONF
*
* Generate determinants in configuration order and obtain
* sign array for switching between the two formats .
*
* Jeppe Olsen January 1989
*
* December 1990 : ICNFOK added
* September 1993 : Combinations added
*
* NCNFCN .ne. 0 indicates that additional constraints on configurations
* should be checked
* by calling CICNCH.ICNFOK(ICNF) is 1 of tests are passed, ICNFOK(ICNF)
* is zero if test fails
* IGENSG .ne. 0 assumes general signs of strings given in ISGNA,ISGNB
      IMPLICIT REAL*8 (A-H,O-Z)
*.Specific input
      DIMENSION ISGNA(*),ISGNB(*)
      DIMENSION ICONF(*),NCNFTP(*)
      DIMENSION IOOS(*)
*.General input
      INTEGER IPRODT(*)
*.Output
      DIMENSION ICTSDT(*)
*.Scratch
#include "WrkSpc.fh"
#include "dmrginfo_mclr.fh"
* NOTE : NCNFTP IS COLUMN FOR SYMMETRY GIVEN , NOT COMPLETE MATRIX.
* Dim of IWORK : MAX(3*NORB,(MXDT+2)*NEL),
* where MXDT is the largest number of prototype determinants of
* a given type.
*
* ================================================================
*. Construct list of configurations,offset for each configuration
*   and type for each configuration
* ================================================================
*

      CALL GETMEM('KSCR11','ALLO','INTEGER',KL1,NORB1+NORB2+NORB3)
      CALL GETMEM('KSCR12','ALLO','INTEGER',KL2,NORB1+NORB2+NORB3)
      CALL GETMEM('KSCR13','ALLO','INTEGER',KL3,NORB1+NORB2+NORB3)
      CALL CONFG2(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,
     *            MINOP,MAXOP,IREFSM,NEL,ICONF,
     *            NCNFTP,IWORK(KL1),IWORK(KL2),IWORK(KL3),IPRNT)
      CALL GETMEM('KSCR13','FREE','INTEGER',KL3,NORB1+NORB2+NORB3)
      CALL GETMEM('KSCR12','FREE','INTEGER',KL2,NORB1+NORB2+NORB3)
      CALL GETMEM('KSCR11','FREE','INTEGER',KL1,NORB1+NORB2+NORB3)
*
* =========================================================
* Obtain determinants for each configuration and determine
* the corresponding adress and phaseshift to reform into
* string form and ordering.
* ==========================================================
*
      CALL CNTOST(ICONF,ICTSDT,NAEL,NBEL,
     &            IPRODT,IREFSM,
     &            NORB,NEL,
     &            IGENSG,ISGNA,ISGNB,ICNSTR,IAGRP,IBGRP,IOOS,PSSIGN,
     &            IPRNT)
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE CNTOST(ICONF,ICTSDT,NAEL,NBEL,
     &                 IPRODT,IREFSM,
     &                 NORB,NEL,
     &                 IGENSG,ISGNA,ISGNB,ICNSTR,IAGRP,IBGRP,IOOS,
     &                 PSSIGN,IPRNT)

*
* Obtain pointer abs(ICTSDT(I)) giving adress of determinant I in
* STRING ordering for determinant I in CSF ordering.
* Going between the two formats can involve a sign change . this is
* stored in the sign of ICTSDT)
* SGNCTS is thus to be multiplied with vector ordered in CSF ordering.
*
* December 1990 : NCNFCN,ICNFOK added
* January 1991  : IGENSG,ISGNA,ISGNB added
* April   1991  : LUCIA version
* September 1993 > Sign and address stored together
*
* ICNSTR .ne. 0 indicates that additional constraints on configurations
* should be checked  (IS = 0 )
* by calling CICNCH.ICNFOK(ICNF) is 1 of tests are passed, ICNFOK(ICNF)
* is zero if test fails
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION ICONF(*),IPRODT(*)
      DIMENSION ICTSDT(*)
      DIMENSION ISGNA(*),ISGNB(*)
      DIMENSION IOOS(*)
C
C IWORK should at least be of length (MXDT+2)*NEL,
C where MXDT is the largest number of prototype determinants occuring
C in a single block.
C
*./SPINFO/
#include "detdim.fh"
#include "spinfo_mclr.fh"
#include "WrkSpc.fh"
#include "dmrginfo_mclr.fh"
*
       NEL = NAEL + NBEL
       NTEST=0000

C.. Local memory
C
C Largest number of dets for a given type
      MXDT = 0
      DO 10 ITYP = 1, NTYP
        MXDT   = MAX(MXDT,NDPCNT(ITYP) )
   10 CONTINUE
      CALL GETMEM('KSCR20','ALLO','INTEGER',KLDTBL,MXDT*NEL)
      CALL GETMEM('KSCR21','ALLO','INTEGER',KLIA,NAEL)
      CALL GETMEM('KSCR22','ALLO','INTEGER',KLIB,NBEL)
      CALL GETMEM('KSCR23','ALLO','INTEGER',KSCR23,NEL)
C
C.. Loop over configurations and generate determinants in compact form
C
      ICNF = 0
      IDT  = 0
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
C Base for prototype determinants
        IF( ITYP .EQ. 1 ) THEN
          IPBAS = 1
        ELSE
          IPBAS = IPBAS + NDPCNT(ITYP-1)*(IOPEN-1)
        END IF
* Determinants for this configuration
        DO 900  IC = 1, NCNATS(ITYP,IREFSM)
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
*. Check orbital occupancy with additional constraints
          IF( NTEST .GE. 10 ) WRITE(6,*) ' IC ICNF ICNBS',IC,ICNF,ICNBS
          CALL CNDET(ICONF(ICNBS),IPRODT(IPBAS),IDET,
     &               NEL,IOCC,IOPEN,ICL, IWORK(KLDTBL),IPRNT)
C Separate determinants into strings and determine string number .
          DO 800 JDET = 1,IDET
!            write(117,"(1X,I8,1X,A,1X)",advance='no')ITYP,"ITYP"  ! yma
            JDTABS = JDTABS + 1
            CALL DETSTR_MCLR(IWORK(KLDTBL+(JDET-1)*NEL),IWORK(KLIA),
     &      IWORK(KLIB),NEL,NAEL,NBEL,NORB,ISIGN,IWORK(KSCR23),IPRNT)
            ijkl_num=ijkl_num+1
* Find number (and sign)of this determinant in string ordering
            ICTSDT(JDTABS) =
     &      IABNUM(IWORK(KLIA),IWORK(KLIB),IAGRP,IBGRP,IGENSG,
     &             ISGNA,ISGNB,ISGNAB,IOOS,NORB,IPSFAC,PSSIGN,
     &             IPRNT)
             IF(  DBLE(ISIGN*ISGNAB*IPSFAC) .eq. -1.0d0)then
               ICTSDT(JDTABS) = - ICTSDT(JDTABS)
             END IF
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
      CALL GETMEM('KSCR23','FREE','INTEGER',KSCR23,NEL)
      CALL GETMEM('KSCR22','FREE','INTEGER',KLIB,NBEL)
      CALL GETMEM('KSCR21','FREE','INTEGER',KLIA,NAEL)
      CALL GETMEM('KSCR20','FREE','INTEGER',KLDTBL,MXDT*NEL)

C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(ICNSTR)
      END
      FUNCTION IABNUM(IASTR,IBSTR,IAGRP,IBGRP,IGENSG,
     &                ISGNA,ISGNB,ISGNAB,IOOS,NORB,IPSFAC,PSSIGN,
     &                IPRNT)
*
* Encapsulation routine for IABNUS
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IASTR(*),IBSTR(*)
      DIMENSION ISGNA(*),ISGNB(*)
#include "detdim.fh"
#include "strinp_mclr.fh"
#include "stinf_mclr.fh"
#include "strbas_mclr.fh"
#include "WrkSpc.fh"
      DIMENSION IOOS(NOCTYP(IAGRP),NOCTYP(IBGRP),*)
*
      IABNUM = IABNUS(IASTR,NELEC(IAGRP),iWORK(KSTREO(IAGRP)),
     &         iWORK(KSTCL(IAGRP)),iWORK(KSTSM(IAGRP)),
     &             NOCTYP(IAGRP),
     &         iWORK(KZ(IAGRP)),iWORK(KISTSO(IAGRP)),
     &         iWORK(KNSTSO(IAGRP)),
     &                IBSTR,NELEC(IBGRP),iWORK(KSTREO(IBGRP)),
     &         iWORK(KSTCL(IBGRP)),iWORK(KSTSM(IBGRP)),NOCTYP(IBGRP),
     &         iWORK(KZ(IBGRP)),iWORK(KISTSO(IBGRP)),
     &     iWORK(KNSTSO(IBGRP)),
     &         IOOS,NORB,IGENSG,ISGNA,ISGNB,ISGNAB,PSSIGN,IPSFAC,
     &         IPRNT)
      RETURN
      END
      FUNCTION IABNUS(IASTR,NAEL,IAORD,ITPFSA,ISMFSA,NOCTPA,ZA,
     &                ISSOA,NSSOA,
     &                IBSTR,NBEL,IBORD,ITPFSB,ISMFSB,NOCTPB,ZB,
     &                ISSOB,NSSOB,
     &                IOOS,NORB,IGENSG,ISGNA,ISGNB,ISGNAB,
     &                PSSIGN,IPSFAC,IPRNT)
*
* A determinant is given by strings IASTR,IBSTR .
* Find number of this determinant
*
* If PSSIGN .ne. 0, the determinant with higher alpha number is picked
* and phase factor IPSFAC calculated. This corresponds to
* configuration order
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IASTR(NAEL),IBSTR(NBEL)
      DIMENSION IAORD(*),IBORD(*)
      INTEGER ZA(*),ZB(*)
      INTEGER NSSOA(NOCTPA,*),NSSOB(NOCTPB,*)
      INTEGER ISSOA(NOCTPA,*),ISSOB(NOCTPB,*)
      INTEGER IOOS(NOCTPA,NOCTPB,*)
      INTEGER ISGNA(*),ISGNB(*)
      INTEGER ITPFSA(*),ITPFSB(*)
      INTEGER ISMFSA(*),ISMFSB(*)
*
* Jeppe Olsen
*
      NTEST =  000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 300) THEN
       WRITE(6,*) ' >>> IABNUS SPEAKING <<< '
       WRITE(6,*) ' NOCTPA,NOCTPB ', NOCTPA,NOCTPB
       WRITE(6,*) ' ALPHA AND BETA STRING '
       CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
       CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
      END IF
*.Number of alpha- and beta-string
C             ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
      IANUM = ISTRNM(IASTR,NORB,NAEL,ZA,IAORD,1)
      IBNUM = ISTRNM(IBSTR,NORB,NBEL,ZB,IBORD,1)
      IF( NTEST .GE. 10 ) WRITE(6,*) ' IANUM AND IBNUM ',IANUM,IBNUM
*
      IF(IGENSG.NE.0) THEN
        ISGNAB = ISGNA(IANUM)*ISGNB(IBNUM)
      ELSE
        ISGNAB = 1
      END IF
*. Symmetries and types
      IASYM = ISMFSA(IANUM)
      IBSYM = ISMFSB(IBNUM)
C?    IF( NTEST .GE.10) WRITE(6,*) ' IASYM IBSYM ',IASYM,IBSYM
      IATP = ITPFSA(IANUM)
      IBTP = ITPFSB(IBNUM)
C?    IF(NTEST.GE.10) WRITE(6,*) ' IATP,IBTP ', IATP,IBTP
      IAREL = IANUM - ISSOA(IATP,IASYM)+1
      IBREL = IBNUM - ISSOB(IBTP,IBSYM)+1
C?    IF(NTEST .GE.10) WRITE(6,*) ' IAREL IBREL ', IAREL,IBREL
*
      IF(PSSIGN.EQ.0.0D0) THEN
*.      Normal determinant ordering
        IABNUS = IOOS(IATP,IBTP,IASYM)
     &         + (IBREL-1)*NSSOA(IATP,IASYM) + IAREL - 1
        IPSFAC = 1
      ELSE IF (PSSIGN .NE. 0.0D0 ) THEN
*.      Ensure mapping to proper determinant in combination
        IF(IANUM.GE.IBNUM) THEN
*.        No need for switching around so
          IF(IASYM.EQ.IBSYM .AND. IATP. EQ. IBTP ) THEN
*.          Lower triangular packed, column wise !
            IABNUS = IOOS(IATP,IBTP,IASYM)  -1
     &             + (IBREL-1)*NSSOA(IATP,IASYM) + IAREL
     &             -  IBREL*(IBREL-1)/2
          ELSE
            IABNUS = IOOS(IATP,IBTP,IASYM)
     &             + (IBREL-1)*NSSOA(IATP,IASYM) + IAREL - 1
          END IF
          IPSFAC = 1
        ELSE IF (IBNUM .GT. IANUM ) THEN
*. Switch alpha and beta string around
          IF(IASYM.EQ.IBSYM .AND. IATP. EQ. IBTP ) THEN
*. Lower triangular packed, column wise !
            IABNUS = IOOS(IBTP,IATP,IBSYM)  -1
     &             + (IAREL-1)*NSSOB(IBTP,IBSYM) + IBREL
     &             -  IAREL*(IAREL-1)/2
          ELSE
            IABNUS = IOOS(IBTP,IATP,IBSYM)
     &             + (IAREL-1)*NSSOB(IBTP,IBSYM) + IBREL
     &             -  1
          END IF
          IPSFAC = nInt(PSSIGN)
        END IF

      END IF
*
COLD
COLD    IABNUS = IOOS(IATP,IBTP,IASYM) + (IBREL-1)*NSSOA(IATP,IASYM)
COLD &           + IAREL - 1
C?    IF(NTEST .GT. 10 ) then
C?      WRITE(6,*) ' IOOS NSSOA ',IOOS(IATP,IBTP,IASYM),
C?   &              NSSOA(IATP,IASYM)
C?    END IF
*
      IF ( NTEST .GE.200) THEN
         WRITE(6,*) ' ALPHA AND BETA STRING '
         CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
         CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
         WRITE(6,*) ' Corresponding determinant number ', IABNUS
      END IF

!      write(117,"(A,1X,I20)",advance='no')" SD_number",IABNUS ! yma
!      write(117,*)

*
      RETURN
      END

*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      Subroutine CsfInf(lSym,iSpin,MS,iSPC,iPrnt,nsym)
*
      Implicit Real*8 (A-H,O-Z)
*
#include "cstate_mclr.fh"
#include "detdim.fh"
#include "strinp_mclr.fh"
#include "strbas_mclr.fh"
#include "stinf_mclr.fh"
#include "orbinp_mclr.fh"
#include "csm.fh"
#include "cicisp_mclr.fh"
#include "csfbas_mclr.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
#include "WrkSpc.fh"
#include "Files_mclr.fh"
      integer idum(1)
*     COMMONBLOCK THE SUPPORT STORAGE OF REORDERING VECTOR
*     ON DISK
#include "csfsd.fh"

*
*
C.... Sorry about this  but this is just to tell the program
C     that no CSF<->SD coefficents is in core
      i1=-9
      iAnders=-9
      ICSM   = lSym
      ISSM   = lSym
      ICSPC  = 1
      ISSPC  = 1
      XNDET   = XISPSM(IsSM,Ispc)
      NEL    = NELCI(ISPC)
      IATP   = IASTFI(ISPC)
      IBTP   = IBSTFI(ISPC)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      MNELR1 = MNR1IC(ISPC)
      MXELR3 = MNR1IC(ISPC)
      iRefSm=lsym
*.... Obtain OOS pointer array
*     strbas_mclr.fh
      CALL GETMEM('SIOIO ','ALLO','INTE',KSIOIO,NOCTPA*NOCTPB)
      CALL IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,
     &            IWORK(KEL1(IATP)),IWORK(KEL3(IATP)),
     &            IWORK(KEL1(IBTP)),IWORK(KEL3(IBTP)),
     &            iWORK(KSIOIO),IPRNT)
      CALL GETMEM('SBLTP','ALLO','INTE',KSBLTP,NSMST)
      NOOS = NOCTPA*NOCTPB*NSMST
      CALL GETMEM('IOOS1','ALLO','INTEGER',KIOOS1,NOOS)
      CALL GETMEM('NOOS1','ALLO','INTEGER',KNOOS1,NOOS)
      CALL INTCSF(NACOB,NEL,iSpin,MS2,
     &            NORB1,NORB2,NORB3,MNELR1,MXELR3,
     &            LLCSF,1,0,PSSIGN,IPRNT,lconf,lldet)
*
*     Calculate CG COEFFICENTS ETC
*
      CALL CSDTMT(IWORK(KDFTP),IWORK(KCFTP),
     &            WORK(KDTOC),PSSIGN,IPRNT)

*
*     Calculate the reordering vector and write it to disk
*

      iA=0
      Do iSym=1,nSym
*.OOS arrayy
          CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,iWORK(KSBLTP),idum)
          CALL ZOOS(ISMOST(1,ISYM),iWORK(KSBLTP),
     &          NSMST,iWORK(KSIOIO),
     &          IWORK(KNSTSO(IATP)),IWORK(KNSTSO(IBTP)),
     &          NOCTPA,NOCTPB,idc,IWORK(KIOOS1),
     &          IWORK(KNOOS1),NCOMB,0)
*EAW           CALL CNFORD(IWORK(KICTS(1)),iWORK(KICONF(1)),
*    &          iSym,NACOB,IWORK(KDFTP),
*    &          NCNATS(1,ISYM),NEL,0,0,IDUM,IDUM,
*    &          IASTFI(ISPC),IBSTFI(ISPC),IWORK(KIOOS1),
*    &          NORB1,NORB2,NORB3,MNELR1,MXELR3,
*    &          NELEC(IASTFI(ISPC)),NELEC(IBSTFI(ISPC)),
*    &          MINOP,MAXOP,IPRNT)
           CALL CNFORD(IWORK(KICTS(1)),iWORK(KICONF(1)),
     &          iSym,NACOB,IWORK(KDFTP),
     &          NCNATS(1,ISYM),NEL,0,0,IDUM,IDUM,
     &          IASTFI(ISPC),IBSTFI(ISPC),IWORK(KIOOS1),
     &          NORB1,NORB2,NORB3,MNELR1,MXELR3,
     &          NELEC(IASTFI(ISPC)),NELEC(IBSTFI(ISPC)),
     &          MINOP,MAXOP,PSSIGN, IPRNT)
*
           CALL iDAFILE(LUCSF2SD,1,IWORK(KICTS(1)),lldet,iA)
           CALL iDAFILE(LUCSF2SD,1,IWORK(KICONF(1)),lconf,iA)
      End Do
*
*
      CALL GETMEM('KICTS','FREE','INTEGER',KICTS(1),NOOS)
      CALL GETMEM('KICONF','FREE','INTEGER',KICONF(1),NOOS)
      CALL GETMEM('IOOS1','FREE','INTEGER',KIOOS1,NOOS)
      CALL GETMEM('NOOS1','FREE','INTEGER',KNOOS1,NOOS)
      CALL GETMEM('SBLTP','FREE','INTE',KSBLTP,NSMST)
*     CALL GETMEM('KDFTP','FREE','INTE',KDFTP,NSMST)
*     strbas_mclr.fh
      CALL GETMEM('SIOIO ','FREE','INTE',KSIOIO,NOCTPA*NOCTPB)
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(MS)
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE CSFDET(NOPEN,IDET,NDET,ICSF,NCSF,CDC,PSSIGN,
     &                  IPRCSF)

* Expand csf's in terms of combinations with
* the use of the Graebenstetter method ( I.J.Q.C.10,P142(1976) )
*
* Input :
*         NOPEN : NUMBER OF OPEN ORBITALS
*         IDET  : OCCUPATION OF combinations
*         NDET  : NUMBER OF combinations
*         ICSF  : INTERMEDIATE SPIN COUPLINGS OF
*                 CSF'S IN BRANCHING DIAGRAM
* Output :
*         CDC :  NDET X NCSF MATRIX
*                GIVING EXPANSION FROM COMB'S TO CSF,S
*                CSF BASIS = Comb basis *CDC
*
* If combinations are use ( signaled by PSSIGN .ne. 0 )
* the factors are multiplies with sqrt(2), corresponding to
* a combination being 1/sqrt(2) times the sum or difference of two
* determinants
*
* The terms are not mutiplied with any sqrt(2), so the transformation is to
* the determinant normalization
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
      DIMENSION IDET(NOPEN,NDET),ICSF(NOPEN,NCSF)
      DIMENSION CDC(NDET,NCSF)

      NTEST = 0000
      NTEST = MAX(IPRCSF,NTEST)
      IF(PSSIGN.EQ.0.0D0) THEN
       CMBFAC = 1.0D0
      ELSE
       CMBFAC = SQRT(2.0D0)
      END IF
      CALL GETMEM('KSCR8','ALLO','REAL',KLMDET,NDET*NOPEN)
      CALL GETMEM('KSCR9','ALLO','REAL',KLSCSF,NDET*NOPEN)
C
C.. OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
      DO 10 JDET = 1, NDET
        CALL MSSTRN(IDET(1,JDET),WORK(KLMDET+(JDET-1)*NOPEN),NOPEN)
   10 CONTINUE

C
      DO 1000 JCSF = 1, NCSF
       IF( NTEST .GE. 105 ) WRITE(6,*) ' ....Output for CSF ',JCSF
C
C OBTAIN INTERMEDIATE COUPLINGS FOR CSF
      CALL MSSTRN(ICSF(1,JCSF),WORK(KLSCSF),NOPEN)
C
      DO 900 JDET = 1, NDET
C EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
      COEF = 1.0D0
      SIGN = 1.0D0
      JDADD = (JDET-1)*NOPEN
      DO 700 IOPEN = 1, NOPEN
C
C + + CASE
        IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.1) THEN
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.0) THEN
C + - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.1) THEN
C - + CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
          SIGN  = - SIGN
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.0) THEN
C - - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
        END IF
  700 CONTINUE
       CDC(JDET,JCSF) = SIGN * CMBFAC * SQRT(COEF)
  900 CONTINUE
 1000 CONTINUE
C
      CALL GETMEM('KSCR9','FREE','REAL',KLSCSF,NDET*NOPEN)
      CALL GETMEM('KSCR8','FREE','REAL',KLMDET,NDET*NOPEN)
*
      IF( NTEST .GE. 5 ) THEN
        WRITE(6,*)
        WRITE(6,'(A,2I2)')
     &  '  The CDC array for  NOPEN ',NOPEN
        WRITE(6,*)
        CALL WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
      END IF

C
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE SPNCOM(iwork,NOPEN,MS2,NDET,IABDET,
     &                  IABUPP,IFLAG,PSSIGN,IPRCSF)
*
* Combinations of nopen unpaired electrons.Required
* spin projection MS2/2.
* JO 21-7-84
*    IFLAG = 1 : Only combinations ( in IABDET )
*    IFLAG = 2 : combinations and upper dets
*    IFLAG = 3 : Only upper dets
* A few revisions october 1988
* Upper dets added feb 1989
* Changed to combinations June 1992
*
* If PSSIGN differs from 0, spin combinations are assumed.
* we select as the unique determinants those with first electron
* having alpha spin
*
      INTEGER ADD
      DIMENSION IABDET(NOPEN,*),IABUPP(NOPEN,*)
      DIMENSION IWORK(*)
      REAL*8 PSSIGN
*
* LENGTH OF IWORK MUST BE AT LEAST NOPEN
*
      NDET=0
      NUPPER = 0
*
* Determinants are considered as binary numbers,1=alpha,0=beta
*
      MX=2 ** NOPEN
      CALL ISETVC(IWORK,0,NOPEN+1)
* Loop over all possible binary numbers
      DO 200 I=1,MX
C.. 1 : NEXT BINARY NUMBER
        ADD=1
        J=0
  190   CONTINUE
        J=J+1
        IF(IWORK(J).EQ.1) THEN
          IWORK(J)=0
        ELSE
          IWORK(J)=1
          ADD=0
        END IF
        IF( ADD .EQ. 1 ) GOTO 190
C
C.. 2 :  CORRECT SPIN PROJECTION ?
        NALPHA=0
        DO 180 J=1,NOPEN
          NALPHA=NALPHA+IWORK(J)
  180   CONTINUE
C
        IF(2*NALPHA-NOPEN.EQ.MS2.AND.
     &    .NOT.(PSSIGN.NE.0.0D0 .AND. IWORK(1).EQ.0)) THEN
          IF (IFLAG .LT. 3 ) THEN
            NDET=NDET+1
            CALL ICOPVE(IWORK,IABDET(1,NDET),NOPEN)
          END IF
          IF (IFLAG .GT. 1 ) THEN
C UPPER DET ?
            MS2L = 0
            LUPPER = 1
C
            DO 10 IEL = 1,NOPEN
              IF (IWORK(IEL).EQ.1) THEN
                 MS2L = MS2L + 1
              ELSE
                 MS2L = MS2L - 1
              END IF
              IF( MS2L .LT. 0 ) LUPPER = 0
   10       CONTINUE
            IF( LUPPER .EQ. 1 ) THEN
              NUPPER = NUPPER + 1
              CALL ICOPVE(IWORK,IABUPP(1,NUPPER),NOPEN)
            END IF
          END IF
        END  IF
C
  200 CONTINUE
C
*     XMSD2=DBLE(MS2)/2
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRCSF)
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE CSDTMT(IDFTP,ICFTP,DTOC,PSSIGN,IPRNT)
*
* Construct list of prototype combinations in IDFTP
* Construct list of prototype CSF'S, in ICFTP
* Construct matrix expanding prototype CSF's in terms of
* prototype combinations in DTOC
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "detdim.fh"
#include "spinfo_mclr.fh"
#include "WrkSpc.fh"
      DIMENSION IDFTP(*),ICFTP(*),DTOC(*)
*./SPINFO/
*

      MULTS = MULTSP
      MS2 = MS2P
*
**.. Set up determinants and upper determinants
*
      IDTBS = 0 ! dummy initialize
      ICSBS = 0 ! dummy initialize
      DO 20 ITP = 1, NTYP
        IOPEN = MINOP+ITP-1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NDPCNT(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NCPCNT(ITP-1)
        END IF
*
        IF( IOPEN .NE. 0 ) THEN
          CALL GETMEM('ISCR7','ALLO','INTEGER',KSCR7,IOPEN+1)
*. Proto type determinants and upper determinants
          IF( MS2+1 .EQ. MULTS ) THEN
            IFLAG = 2
            CALL SPNCOM(iwork(kscr7),IOPEN,MS2,NNDET,IDFTP(IDTBS),
     &                  ICFTP(ICSBS),IFLAG,PSSIGN,IPRNT)
          ELSE
            IFLAG = 1
            CALL SPNCOM(iwork(kscr7),IOPEN,MS2,NNDET,IDFTP(IDTBS),
     &                  ICFTP(ICSBS),IFLAG,PSSIGN,IPRNT)
            IFLAG = 3
            CALL SPNCOM(iwork(kscr7),IOPEN,MULTS-1,NNDET,IDFTP(IDTBS),
     &                  ICFTP(ICSBS),IFLAG,PSSIGN,IPRNT)
          END IF
          CALL GETMEM('ISCR7','FREE','INTEGER',KSCR7,IOPEN+1)
        END IF
   20 CONTINUE
*. Matrix expressing csf's in terms of combinations
      ICDCBS =0 ! dummy initialize
      DO 30 ITP = 1, NTYP
        IOPEN = MINOP+ITP-1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
          ICDCBS =1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NDPCNT(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NCPCNT(ITP-1)
          ICDCBS = ICDCBS + NDPCNT(ITP-1)*NCPCNT(ITP-1)
        END IF
        IF(NDPCNT(ITP)*NCPCNT(ITP).EQ.0) GOTO 30
        IF(IOPEN .EQ. 0 ) THEN
          DTOC(ICDCBS) = 1.0D0
        ELSE
          CALL CSFDET(IOPEN,IDFTP(IDTBS),NDPCNT(ITP),
     &               ICFTP(ICSBS),NCPCNT(ITP),DTOC(ICDCBS),
     &               PSSIGN,IPRNT)
        END IF
   30 CONTINUE
*

      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE INTCSF(NACTOB,NACTEL,MULTP,MS2,NORB1,NORB2,NORB3,
     &                  NEL1MN,NEL3MX,
     &                  LLCSF,NCNSM,ICNSTR,PSSIGN,
     &IPRNT,lconf,lldet)

*
* Initializing routine for CSF-DET expansions of internal space
*
* Set up common block /CSFDIM/
* This gives information about the number of dets,csf's for
* each symmetry
*
* find local memory requirements for CSF routines
* Largest local memory requirements in CNFORD,CSFDET is returned in
* LLCSF
*
* DETERMINE BASE ADRESSES /CSFBAS/
*             KDFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
*             KCFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
*             KDTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
*             KICONF(I) : SPACE FOR STORING  NCNSM
*                        CONFIGURATION EXPANSIONS
*
* If PSSIGN .ne. 0, spin combinations are used !!

      IMPLICIT REAL*8 (A-H,O-Z)
#include "detdim.fh"
#include "spinfo_mclr.fh"
#include "WrkSpc.fh"
#include "csfbas_mclr.fh"
#include "dmrginfo_mclr.fh"
*
* Last modification : Sept 20 : sign and address of dets goes together
*                      in KICTS
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
      IF(PSSIGN.NE.0.0D0) THEN
       IMSCMB = 1
      ELSE
       IMSCMB = 0
      END IF
*
**. Define parameters in SPINFO
*
      MULTSP = MULTP
      MS2P = MS2
      MULTS = MULTSP

      NEL = NACTEL
* ================================
*. Allowed number of open orbitals
* ================================
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
C?    WRITE(6,*) ' MAXOP with RAS constraints :' ,MAXOP
      NTYP = MAXOP-MINOP + 1
*
      IF( NTYP .GT. MXPCTP ) THEN
        WRITE(6,*) '  NUMBER OF CONFIGURATION TYPES TO LARGE '
        WRITE(6,*) '  CHANGE PARAMETER MXPCTP TO AT LEAST ',NTYP
        WRITE(6,*) '  CURRENT VALUE OF MXPCTP ',MXPCTP
        write(6,*) ' MTYP IN LUSPIN TO SMALL '
        call Abend()
      END IF
*
      IF( NTEST .GE. 5 )
     &WRITE(6,*) ' MINOP MAXOP NTYP ',MINOP,MAXOP,NTYP
* ================================================
*. Number of sd's and csf's per configuration type
* ================================================
      DO 10 ITP = 1,NTYP
        IOPEN = MINOP+ITP - 1
        IAEL = (IOPEN + MS2 ) / 2
        IBEL = (IOPEN - MS2 ) / 2
        IF(IAEL+IBEL .EQ. IOPEN ) THEN
          NDPCNT(ITP) = IBION(IOPEN,IAEL)
          IF(IMSCMB.NE.0.AND.IOPEN.NE.0)
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
      WRITE(6,'(/A)')
     &'  Open orbitals   Combinations    CSFs '
      DO 580 IOPEN = MINOP,MAXOP,2
        ITYPE = IOPEN - MINOP + 1
        WRITE(6,'(5X,I3,10X,I6,7X,I6)')
     &  IOPEN,NDPCNT(ITYPE),NCPCNT(ITYPE)
  580 CONTINUE
      END IF
* =================================================
**. Number of Combinations and CSF's per  symmetry
* =================================================
*EAW
      CALL GETMEM('KIICL','ALLO','INTEGER',KL1,NACTOB)
      CALL GETMEM('KIIOP','ALLO','INTEGER',KL2,NACTOB)
      CALL GETMEM('KIIOC','ALLO','INTEGER',KL3,NORB1+NORB2+NORB3)
*     KL1 = 1
*     KL2 = KL1 + NACTOB
*     KL3 = KL2 + NACTOB
*EAW

      CALL CISIZE(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,NACTEL,
     &            MINOP,MAXOP,MXPCTP,MXPCSM,NCNATS,NCNASM,NDTASM,
     &            NCSASM,
     &            NDPCNT,NCPCNT,
     &            IWORK(KL1),IWORK(KL2),IWORK(KL3),IPRNT)

*EAW
      CALL GETMEM('KIIOC','FREE','INTEGER',KL3,NORB1+NORB2+NORB3)
      CALL GETMEM('KIIOP','FREE','INTEGER',KL2,NACTOB)
      CALL GETMEM('KIICL','FREE','INTEGER',KL1,NACTOB)
*EAW
* ==============================================
*   Permanent and local memory for csf routines
* ==============================================
*
*    memory for CSDTMT arrays.
*    Largest block of proto type determinants .
*    Largest number of prototype determinants
*    All configurations( of any specific symmetry )
*

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
*. local memory for CSFDET
      LCSFDT = MXPTBL + MAXOP
*. local memory for CNFORD
      LCNFOR = MAX(2*NTYP+NACTOB,(MXDT+2)*NACTEL)
*. local memory for any routine used in construction of csf basis
      LLCSF = MAX(LCSFDT,LCNFOR)
*. Memory needed to store ICONF array
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
C?      WRITE(6,*) ' MEMORY FOR HOLDING CONFS OF SYM... ',ISYM,LLCONF
        LCONF = MAX(LCONF,LLCONF)
        ILCNF = MAX(ILCNF,ILLCNF)
   30 CONTINUE

      ! notice the ILCNF number ! yma

       IF(NTEST.GE.5) THEN
       WRITE(6,'(/A,I8)')
     & '  Memory for holding largest list of configurations ',LCONF
       WRITE(6,'(/A,I8)')
     & '  Size of largest CI expansion (combinations)',LDET
       WRITE(6,'(/A,I8)')
     & '  Size of largest CI expansion (confs)',ILCNF
       END IF
       call xflush(6) !yma

!       if(doDMRG.and.doMCLR)then  ! yma
!         ndets_RGLR=LDET
!         ncsfs_RGLR=ILCNF
!       end if

*. permanent memory for csf proto type arrays

      CALL GETMEM('KDFTP','ALLO','INTEGER',KDFTP,LIDT)
      CALL GETMEM('KCFTP','ALLO','INTEGER',KCFTP,LICS)
      CALL GETMEM('KDTOC','ALLO','REAL   ',KDTOC,LDTOC)

*. Permanent arrays for reordering and phases
      IF(NCNSM .GT. MXCNSM ) THEN
        WRITE(6,'(A,2I2)')
     &  '  TROUBLE IN CSFDIM NCNSM > MXCNSM : NCNSM,MXCNSM',
     &  NCNSM,MXCNSM
        write(6,*) ' CSFDIM : NCNSM  IS GREATER THAN MXCNSM '
        Call Abend()
      END IF
      DO 60 ICNSM = 1, NCNSM
        CALL GETMEM('KICONF','ALLO','INTEGER',KICONF(ICNSM),LCONF)
        CALL GETMEM('KICTS ','ALLO','INTEGER',KICTS(ICNSM),LDET)
   60 CONTINUE
*
*
*
      lldet=ldet

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(ICNSTR)
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE CISIZE(NORB1,NORB2,NORB3,
     &                  NEL1MN,NEL3MX,NACTEL,
     &                  MINOP,MAXOP,
     &                  MXPCNT,MXPCSM,
     &                  NCNATS,NCNASM,NDTASM,NCSASM,
     &                  NDPCNT,NCPCNT,
     &                  IICL,IIOP,IIOC,IPRNT)
*
*   Number of configurations per per configuration type and symmetry
*
*  Jeppe Olsen
*         August 1990 : Improved handling of large RAS 3 space
*         Winter 1991 : Modified for LUCIA
      IMPLICIT REAL*8 (A-H,O-Z)
*.Input
      DIMENSION NDPCNT(*),NCPCNT(*)
*.Output
      DIMENSION NCNATS(MXPCNT,*),NCNASM(*),NDTASM(*),NCSASM(*)
*. Scratch
      DIMENSION IICL(*),IIOP(*),IIOC(NORB1+NORB2+NORB3)
*
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      ILOOP = 0
      ILOOP2 = 0
      NCNF = 0

*
      CALL iCOPY(MXPCSM*MXPCNT,[0],0,NCNATS,1)
      CALL iCOPY(MXPCSM,[0],0,NCSASM,1)
      CALL iCOPY(MXPCSM,[0],0,NDTASM,1)
      CALL iCOPY(MXPCSM,[0],0,NCNASM,1)

*
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
*
      IORB2F = IORB1L + 1
      IORB2L = IORB2F + NORB2 - 1
*
      IORB3F = IORB2L + 1
      IORB3L = IORB3F + NORB3 - 1
*
      NORB = NORB1 + NORB2 + NORB3
* Min number of doubly occupied orbitals in RAS 1
      MINCL1 = MAX(0,NEL1MN-NORB1)
      IF(NTEST.GE.1)  WRITE(6,*)
     &  ' Min number of doubly occupied orbitals in RAS 1',MINCL1
      DO 5000 NOP = MINOP,MAXOP,2
        ITYPE = NOP-MINOP+1
        NCL = (NACTEL-NOP)/2
        IF( NTEST .GE. 10 )
     &  WRITE(6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
C. first combination of double occupied orbitals
        CALL iCOPY(NORB,[0],0,IIOC,1)
        DO 10 ICL = 1, NCL
          IICL(ICL) = ICL
          IIOC(ICL) = 2
   10   CONTINUE
        IFRSTC = 1
C.. Loop over double occupied orbital configurations
 2000   CONTINUE
C
C. next double occupied configuration
          IF ( IFRSTC .EQ. 1. OR. NCL .EQ. 0 ) GOTO 801
C         IF ( IFRSTC .EQ. 0 .AND. NCL .NE. 0 ) THEN
C
          DO 50 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 1 ) IIOC(IORB) = 0
   50     CONTINUE
C
          IPLACE = 0
  800     IPLACE = IPLACE + 1

          IPRORB = IICL(IPLACE)
          IIOC(IPRORB) = 0
          NEWORB = IPRORB+1
          IF((IPLACE .LT. NCL .AND. NEWORB .LT. IICL(IPLACE+1))
     &      .OR.
     &      IPLACE .EQ. NCL .AND.  NEWORB.LE. NORB ) THEN
            IICL(IPLACE) = NEWORB
            IIOC(NEWORB) = 2
          ELSE IF
     &    (.NOT.(IPLACE.EQ.NCL.AND.NEWORB.GE.NORB)) THEN
            IF(IPLACE .EQ. 1 ) THEN
              IICL(1) = 1
              IIOC(1) = 2
            ELSE
              IICL(IPLACE) = IICL(IPLACE-1) + 1
              IIOC(IICL(IPLACE)) = 2
            END IF
            GOTO 800
          ELSE
C. No more inactive configurations
             GOTO 2001
          END IF
COLD      END IF
  801     CONTINUE
          IFRSTC = 0
          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next inactive configuration '
            CALL IWRTMA(IICL,1,NCL,1,NCL)
          END IF
C..         CHECK RAS1 and RAS 3
            IEL1C = 0
            IEL3C = 0
            ICL1  = 0
            DO  20 ICL = 1,NCL
              IORB = IICL(ICL)
              IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
                 IEL1C = IEL1C + 2
                 ICL1 = ICL1 + 1
              ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
                 IEL3C = IEL3C + 2
              END IF
   20       CONTINUE
            IIICHK = 1
            IF(ICL1 .LT.MINCL1.AND. IIICHK.EQ.1) THEN
* Next higher combination with a higher number of inactive orbitals
              DO 41 ICL = 1,ICL1+1
                IIOC(IICL(ICL)) = 0
                IICL(ICL) = ICL
                IIOC(ICL) = 2
   41         CONTINUE
              IPLACE=ICL1+1
              IF( IPLACE.GE.NCL) GOTO 2001
              GOTO 800
            END IF
            IF(IEL3C.GT.NEL3MX) GOTO 2000
C. Highest orbital not occupied
         MXMPTY = NORB
         IORB = NORB+1
C. begin while
   12      CONTINUE
           IORB = IORB - 1
           IF(IIOC(IORB) .EQ. 2 ) THEN
             MXMPTY = IORB-1
             IF( IORB .NE. 1 ) GOTO 12
           END IF
C. End while
C
C. first active configuration
          IORB = 0
          IOP = 0
          DO 30 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 0 ) THEN
              IOP = IOP + 1
              IF( IOP .GT. NOP ) GOTO 31
              IIOC(IORB) = 1
              IIOP(IOP) = IORB
            END IF
   30     CONTINUE
   31     CONTINUE
          IFRSTO = 1
C
C. Next open shell configuration
 1000     CONTINUE
            IF(IFRSTO.EQ.1. OR. NOP .EQ. 0 ) goto 701
COLD        IF(IFRSTO .EQ. 0 .AND. NOP .NE. 0 ) THEN
            IPLACE = 0
  700       CONTINUE
              IPLACE = IPLACE + 1
              IPRORB = IIOP(IPLACE)
              NEWORB = IPRORB + 1
              IIOC(IPRORB) = 0
  690         CONTINUE
                IF(NEWORB .LE. MXMPTY .AND. IIOC(NEWORB) .NE. 0) THEN
                   NEWORB = NEWORB + 1
                   IF(NEWORB .LE. MXMPTY ) GOTO 690
                END IF
C 691         End of loop
            IF(IPLACE .LT. NOP .AND. NEWORB .LT. IIOP(IPLACE+1)
     &      .OR.
     &      IPLACE .EQ. NOP .AND. NEWORB .LE. MXMPTY ) THEN
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
           ELSE IF (  IPLACE .NE. NOP ) THEN
              IF(IPLACE.EQ.1) THEN
                NEWORB = 1 - 1
              ELSE
                NEWORB = IIOP(IPLACE-1)
              END IF
 671          CONTINUE
                NEWORB = NEWORB + 1
              IF(IIOC(NEWORB) .NE. 0 .AND. NEWORB.LE.MXMPTY ) GOTO 671
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
              GOTO 700
           ELSE
C. No more active configurations , so
              IF( NCL .NE. 0 ) GOTO 2000
              IF( NCL .EQ. 0 ) GOTO 5001
           END IF
C          END IF
  701      CONTINUE
           IFRSTO = 0

          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next active configuration '
            CALL IWRTMA(IIOP,1,NOP,1,NOP)
          END IF
C        RAS  CONSTRAINTS
         IOKAY=1
           IEL1 = IEL1C
           IEL3 = IEL3C
C..        CHECK RAS1 and RAS3
           DO  40 IOP = 1,NOP
             IORB = IIOP(IOP)
             IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
                IEL1 = IEL1 + 1
             ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
                IEL3 = IEL3 + 1
             END IF
   40      CONTINUE
           IR3CHK  = 1
           IF(IEL3.GT.NEL3MX.AND.IR3CHK.EQ.1) THEN
*. Number of electrons in substring
             IFSTR3 = 0
             DO 5607 IOP = 1, NOP
               IF(IIOP(IOP).GE.IORB3F) THEN
                IFSTR3 = IOP
                GOTO 5608
               END IF
 5607        CONTINUE
 5608        CONTINUE
             IF(IFSTR3.NE.NOP) THEN

*. Lowest possible string with NOP electrons
               DO 5610 K = 1, IFSTR3
                 IIOC(IIOP(K)) = 0
 5610          CONTINUE
*
               KEL = 0
               KORB = 0
 5630          CONTINUE
                 KORB = KORB + 1
                 IF(IIOC(KORB).NE.2) THEN
                   KEL = KEL + 1
                   IIOC(KORB) = 1
                   IIOP(KEL) = KORB
                 END IF
               IF(KEL.NE.IFSTR3) GOTO  5630
               IPLACE = IFSTR3
               GOTO 700
             END IF
           END IF
           IF( IEL1 .LT. NEL1MN .OR. IEL3 .GT. NEL3MX ) GOTO  999
C. Spatial symmetry
         ISYM = ISYMCN_MCLR(IICL,IIOP,NCL,NOP)
C
           IF( NTEST .GE.2000)
     &     WRITE(6,*) ' ISYM : ', ISYM
           IF(NTEST.GE.1500)
     &     WRITE(6,1120) ( IIOC(I),I = 1,NORB )
 1120      FORMAT(1H0,'  configuration included ',15I3,
     &               ('                         ',15I3))
           NCNF=NCNF+1

           NCNASM(ISYM) = NCNASM(ISYM)+1
           IF(NTEST.GE.1500 )
     &     WRITE(6,1311) NCNF,(IIOC(I),I=1,NORB)
 1311      FORMAT(1H ,' configuration ',I3,20I2,/,(1H ,18X,20I2))

           NCNATS(ITYPE,ISYM)=NCNATS(ITYPE,ISYM)+1
           IF(NTEST.GE.2000) WRITE(6,3111) NCNF,ITYPE
 3111      FORMAT(1H0,'  CONFIGURATION..',I3,' IS TYPE..',I3)
C
C** LOOP OVER CONFIGURATIONS, end
C
  999  CONTINUE
       ILOOP = ILOOP + 1
       ILOOP2 = ILOOP2 + 1
       IF(ILOOP2 .EQ. 10 000 000 ) THEN
         WRITE(6,*) ' 10 million configurations generated '
         ILOOP2 = 0
       END IF
*
       IF( NOP .EQ. 0 .AND. NCL .EQ. 0 ) GOTO 5001
       IF( NOP .EQ. 0 ) GOTO 2000
       GOTO 1000
 2001 CONTINUE
 5000 CONTINUE
 5001 CONTINUE

      IF( NTEST .GE. 2 ) THEN
        WRITE(6,'(A,I8)')
     &  '  Total number of configurations generated ', NCNF
      END IF
* ================================
*. Total number of CSF's and SD's
* ================================
      NTYP = MAXOP - MINOP + 1
      DO 610 ISYM = 1, MXPCSM
        DO 600 ITYPE = 1, NTYP
          NDTASM(ISYM) = NDTASM(ISYM)
     &    + NDPCNT(ITYPE)*NCNATS(ITYPE,ISYM)
          NCSASM(ISYM) = NCSASM(ISYM)
     &    + NCPCNT(ITYPE)*NCNATS(ITYPE,ISYM)
  600   CONTINUE
  610 CONTINUE

C
      IF( NTEST .GE. 2 ) THEN
       WRITE(6,*)
       ICSM = 0
       WRITE(6,'(/A)') ' Information about actual configurations '
       WRITE(6,'( A)') ' ========================================'
       WRITE(6,'(/A)')
     & '    Symmetry     Configurations     CSFs     Combinations  '
       WRITE(6,'(A)')
     & '  =============  ============== ============ ============  '
       DO 570 ICSM = 1, MXPCSM
         IF(NCNASM(ICSM) .NE. 0 ) THEN
            WRITE(6,'(4X,I3,4X,6X,I8,6X,I8,6X,I9)')
     &      ICSM ,NCNASM(ICSM),NCSASM(ICSM),NDTASM(ICSM)
         END IF
  570  CONTINUE
      END IF

*
*
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE CNDET_MCLR(ICONF,IPDET,NDET,NEL,NORB,NOP,NCL,
     *                 IDET,IPRNT)
C
C A configuration ICONF in compressed form and a set of
C prototype determinants ,IPDET, are given .
C
C Construct the corresponding determinants in contracted  form .
C
C JEPPE OLSEN , NOVEMBER 1988
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ICONF(*   )
C     DIMENSION IPDET(NOP,NDET)
      DIMENSION IPDET(*       )
      DIMENSION IDET(NEL,*   )
C
C
C POSITIVE NUMBER  : ALPHA ORBITAL
C NEGATIVE NUMBER  : BETA  ORBITAL
C
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
      IF( NTEST .GT.200 ) THEN
        IF(NCL .NE. 0 ) THEN
          WRITE(6,*) ' DOUBLE OCCUPIED ORBITALS '
          CALL IWRTMA(ICONF,1,NCL,1,NCL)
        END IF
        IF(NOP .NE. 0 ) THEN
          WRITE(6,*) ' OPEN ORBITALS '
          CALL IWRTMA(ICONF(1+NCL),1,NOP,1,NOP)
        END IF
      END IF
C
C.. 1 DOUBLY OCCUPIED ORBITALS ARE PLACED FIRST
C
      DO 100 ICL = 1, NCL
        IBASE = 2 * (ICL-1)
        DO 90 JDET = 1, NDET
          IDET(IBASE+1,JDET) =  ICONF(ICL)
          IDET(IBASE+2,JDET) = -ICONF(ICL)
90      CONTINUE
100   CONTINUE
C
C..2  SINGLY OCCUPIED ORBITALS
C
      IADD = 2*NCL
      DO 200 JDET = 1, NDET
        DO 190 IOP = 1, NOP
          IADR = (JDET-1)*NOP + IOP
          IF( IPDET(IADR    ) .EQ. 1 ) IDET(IADD+IOP,JDET) =
     *    ICONF(NCL +IOP)
          IF( IPDET(IADR    ) .EQ. 0 ) IDET(IADD+IOP,JDET) =
     *    - ICONF(NCL +IOP)
190     CONTINUE
200   CONTINUE
C
C..3 OUTPUT
C
      IF( NTEST.GE.200) THEN
       WRITE(6,*) ' CONFIGURATION FROM DETCON '
       CALL IWRTMA(ICONF,1,NORB,1,NORB)
       WRITE(6,* ) ' PROTO TYPE DETERMINANTS '
       IF(NOP*NDET .GT. 0)
     * CALL IWRTMA(IPDET,NOP,NDET,NOP,NDET)
       IF(NEL*NDET .GT. 0 )
     * WRITE(6,*) ' CORRESPONDING DETERMINANTS '
       CALL IWRTMA(IDET,NEL,NDET,NEL,NDET)
      END IF
C
C..4  EXIT
C
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE DETSTR_MCLR(IDET,IASTR,IBSTR,NEL,NAEL,NBEL,NORB,
     *ISIGN,IWORK,IPRNT)

C
C A DETERMINANT,IDET,IS GIVEN AS A SET OF OCCUPIED SPIN ORBITALS,
C POSITIVE NUMBER INDICATES ALPHA ORBITAL AND NEGATIVE NUMBER
C INDICATES BETA ORBITAL .
C
C FIND CORRESPONDING ALPHA STRING AND BETA STRING ,
C AND DETERMINE SIGN NEEDED TO CHANGE DETERMINANT
C INTO PRODUCT OF ORDERED ALPHA STRING AND
C BETA STRING
C
C JEPPE OLSEN NOVEMBER 1988
C
      IMPLICIT REAL*8 (A-H,O-Z)

#include "dmrginfo_mclr.fh"
      DIMENSION IDET(NEL)
      DIMENSION IASTR(NAEL),IBSTR(NBEL)
      DIMENSION IWORK(*)
C
C
C
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
C
C FIRST REORDER SPIN ORBITALS IN ASCENDING SEQUENCE
C THIS WILL AUTOMATICALLY SPLIT ALPHA AND BETASTRING
C
      CALL ORDSTR_MCLR(IDET,IWORK,NEL,ISIGN,IPRNT)
C
C ALPHA STRING IS LAST NAEL ORBITALS
      CALL iCOPY(NAEL,IWORK(NBEL+1),1,IASTR,1)
C
C BETA  STRING MUST BE COMPLETELY TURNED AROUND
      DO 10 IBEL = 1, NBEL
        IBSTR(IBEL) = -IWORK(NBEL+1-IBEL)
10    CONTINUE
C SIGN CHANGE FOR SWITCH OF BETA ORBITALS
      iTmp= NBEL*(NBEL+1)/2
      ISIGN = ISIGN * (-1) ** iTmp
C
      IF( NTEST.GE.200) THEN
        WRITE(6,*) ' INPUT DETERMINANT '
        CALL IWRTMA(IDET,1,NEL,1,NEL)
        WRITE(6,*) ' CORRESPONDING ALPHA STRING '
        CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
        WRITE(6,*) ' CORRESPONDING BETA STRING '
        CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
        WRITE(6,*) ' ISIGN FOR SWITCH ', ISIGN
      END IF

!      if(doDMRG.and.doMCLR)then ! yma
!        DO I=1,NEL
!          Write(117,1110,advance='no') IDET(I)
!        end do
!        Write(117,"(A,1X,I2)",advance='no')" SIGN",ISIGN
!1110  FORMAT(1X,I5)
!      end if

C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NORB)
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE ORDSTR_MCLR(IINST,IOUTST,NELMNT,ISIGN,IPRNT)
C
C ORDER A STRING OF INTEGERS TO ASCENDING ORDER
C
C IINST : INPUT STRING IS IINST
C IOUTST : OUTPUT STRING IS IOUTST
C NELMNT : NUMBER OF INTEGERS IN STRING
C ISIGN :  SIGN OF PERMUTATION : + 1 : EVEN PERMUTATIONN
C                                - 1 : ODD  PERMUTATION
C
C THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
C ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
C
C IMPLEMENTED MORE TRANSPARENT BUBBLE SORTING INSTEAD
C               JR NOV 2006
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IINST(NELMNT),IOUTST(NELMNT)
      IF(NELMNT.EQ.0) RETURN
      ISIGN=1
      iTEMP=0

10       iPass=0
       DO I=1,NELMNT-1
       IF(IINST(I).GT.IINST(I+1)) THEN
         iTEMP=IINST(I)
         IINST(I)=IINST(I+1)
         IINST(I+1)=iTEMP
         ISIGN=-1*ISIGN
         iPass=1
         ENDIF
       ENDDO
       IF(IPASS.NE.0) GOTO 10


              DO I=1,NELMNT
              IOUTST(I)=IINST(I)
        ENDDO
C
C      NTEST = 0
C      NTEST = MAX(NTEST,IPRNT)
C
C      CALL iCOPY(NELMNT,IINST,1,IOUTST,1)
C      ISIGN = 1
C
C
c        JOE = 1
c10      I = JOE
c20      CONTINUE
c        IF(I.EQ.NELMNT) GO TO 50
c        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
c        JOE = I + 1
C30      iSWAP = IOUTST(I)
c        ISIGN = - ISIGN
c        IOUTST(I) = IOUTST(I+1)
c        IOUTST(I+1) = iSWAP
c       IF(I.EQ.1) GO TO 10
c        I = I - 1
c        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
c        GO TO 10
c40      I = I + 1
c      GO TO 20
C
C     END ORDER
C
c50    CONTINUE
c      IF( NTEST .GE.200) THEN
c        WRITE(6,*)  ' INPUT STRING ORDERED STRING ISIGN ',NELMNT
c        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
c        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
c        WRITE(6,*) ' ISIGN : ', ISIGN
c      END IF
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRNT)
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE MSSTRN(INSTRN,UTSTRN,NOPEN)
*
* A STRING IS GIVEN IN FORM A SEQUENCE OF ZEROES
* AND ONE ' S
*
* REINTERPRET THIS AS :
*
* 1 : THE INPUT STRING IS A DETERMINANT AND THE
*     1'S INDICATE ALPHA ELECTRONS AND THE
*     0'S INDICATE BETA ELECTRONS .
*     UTSTRN IS THE MS-VALUES ATE EACH VERTEX
*
* 2 : THE INPUT STRING IS A CSF GIVEN IN A
*     BRANCHING DIAGRAM, WHERE
*     1'S INDICATE UPWARDS SPIN COUPLEING
*     WHILE THE 0'S INDICATES DOWNWARDS SPIN COUPLING ,
*     REEXPRESS THIS AS S VALUES OF ALL COUPLINGS
*
* THE TWO PROCEDURES ARE IDENTICAL .
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION INSTRN(NOPEN),UTSTRN(NOPEN)
*
      UTSTRN(1) = DBLE(INSTRN(1)) - 0.5D0
      DO 10 IOPEN = 2, NOPEN
        UTSTRN(IOPEN) = UTSTRN(IOPEN-1) +DBLE(INSTRN(IOPEN))-0.5D0
10    CONTINUE
*
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE CONFG2 (NORB1,NORB2,NORB3,
     &                   NEL1MN,NEL3MX,
     &                   MINOP,MAXOP,
     &                   IREFSM,NEL,ICONF,
     &                   NCNFTP,IIOC,IIOP,IICL,IPRNT)
*
*  Generate array,ICONF,giving occupation of each configuration
*  for CI space of reference symmetry IREFSM.
*
*
*  Jeppe Olsen April 1989
*              August 1990 : Improved handling of large RAS 3 space
*
*  Turbo configuration generator
*  Iconf is ordered so all configuratiuns of the same type are
*  consecutively stored .
*  ICONF is written so closed orbitals are given first and then single
*  occupied orbitals
*
      IMPLICIT REAL*8 (A-H,O-Z)
*.Input
      DIMENSION NCNFTP(*)
*.Output
      DIMENSION ICONF(*)
*.Scratch
      DIMENSION IIOC(*),IICL(*),IIOP(*)
*
*
      NTEST = 0000

      NTEST = MAX(NTEST,IPRNT)
*
      IORB1F = 1
      IORB1L = IORB1F+NORB1-1
*
      IORB2F = IORB1L + 1
      IORB2L = IORB2F + NORB2 - 1
*
      IORB3F = IORB2L + 1
      IORB3L = IORB3F + NORB3 - 1
*
      NORB = NORB1 + NORB2 + NORB3
C
C.. Loop over types of configurations
C
      JCONF =  0
      ICFREE = 1
* Min number of doubly occupied orbitals in RAS 1
      MINCL1 = MAX(0,NEL1MN-NORB1)
      IF(NTEST.GE.1) WRITE(6,*)
     &  ' Min number of doubly occupied orbitals in RAS 1',MINCL1
      DO 5000 NOP = MINOP,MAXOP,2
        ITYPE = NOP-MINOP+1
        NCL = (NEL-NOP)/2
        IF( NTEST .GE. 10 )
     &  WRITE(6,*) ' NOP NCL ITYPE',NOP,NCL,ITYPE
C
C. first combination of double occupied orbitals
        CALL iCOPY(NORB,[0],0,IIOC,1)
        DO 10 ICL = 1, NCL
          IICL(ICL) = ICL
          IIOC(ICL) = 2
   10   CONTINUE
        IFRSTC = 1
C.. Loop over double occupied orbital configurations
 2000   CONTINUE
C
C. next double occupied configuration
          IF( IFRSTC .EQ. 1 .OR. NCL .EQ. 0 ) GOTO 801
C
          DO 50 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 1 ) IIOC(IORB) = 0
   50     CONTINUE
C
          IPLACE = 0
  800     IPLACE = IPLACE + 1

          IPRORB = IICL(IPLACE)
          IIOC(IPRORB) = 0
          NEWORB = IPRORB+1
          IF((IPLACE .LT. NCL .AND. NEWORB .LT. IICL(IPLACE+1))
     &      .OR.
     &      IPLACE .EQ. NCL .AND.  NEWORB.LE. NORB ) THEN
C
            IICL(IPLACE) = NEWORB
            IIOC(NEWORB) = 2
          ELSE IF
     &    (.NOT.(IPLACE.EQ.NCL.AND.NEWORB.GE.NORB)) THEN
C
            IF(IPLACE .EQ. 1 ) THEN
              IICL(1) = 1
              IIOC(1) = 2
            ELSE
              IICL(IPLACE) = IICL(IPLACE-1) + 1
              IIOC(IICL(IPLACE)) = 2
            END IF
            GOTO 800
          ELSE
C. No more inactive configurations
             GOTO 2001
          END IF
  801     CONTINUE
          IFRSTC = 0
C..         CHECK RAS1 and RAS 3
            IEL1C = 0
            IEL3C = 0
            ICL1  = 0
            DO  20 ICL = 1,NCL
              IORB = IICL(ICL)
              IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
                 IEL1C = IEL1C + 2
                 ICL1 = ICL1 + 1
              ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
                 IEL3C = IEL3C + 2
              END IF
   20       CONTINUE
            IIICHK = 1
            IF(ICL1 .LT.MINCL1.AND. IIICHK.EQ.1) THEN
* Next higher combination with a higher number of inactive orbitals
              DO 41 ICL = 1,ICL1+1
                IIOC(IICL(ICL)) = 0
                IICL(ICL) = ICL
                IIOC(ICL) = 2
   41         CONTINUE
              IPLACE=ICL1+1
              IF( IPLACE.GE.NCL) GOTO 2001
              GOTO 800
            END IF
            IF(IEL3C.GT.NEL3MX) GOTO 2000
C. Highest orbital not occupied
         MXMPTY = NORB
         IORB = NORB+1
C. begin while
   12      CONTINUE
           IORB = IORB - 1
           IF(IIOC(IORB) .EQ. 2 ) THEN
             MXMPTY = IORB-1
             IF( IORB .NE. 1 ) GOTO 12
           END IF
C. End while
          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next inactive configuration '
            CALL IWRTMA(IICL,1,NCL,1,NCL)
          END IF
C
C. first active configuration
          IORB = 0
          IOP = 0
          DO 30 IORB = 1, NORB
            IF(IIOC(IORB) .EQ. 0 ) THEN
              IOP = IOP + 1
              IF( IOP .GT. NOP ) GOTO 31
              IIOC(IORB) = 1
              IIOP(IOP) = IORB
            END IF
   30     CONTINUE
   31     CONTINUE
          IFRSTO = 1
C
C. Next open shell configuration
 1000     CONTINUE
            IF(IFRSTO. EQ. 1 .OR. NOP .EQ. 0 ) GOTO 701
            IPLACE = 0
  700       CONTINUE
              IPLACE = IPLACE + 1
              IPRORB = IIOP(IPLACE)
              NEWORB = IPRORB + 1
              IIOC(IPRORB) = 0

* PAM 2013: Searching for next orbital with IIOC=0:
  690         CONTINUE
              IF(NEWORB.LE.MXMPTY .and. IIOC(NEWORB) .NE. 0) THEN
                NEWORB = NEWORB + 1
                GOTO 690
              END IF

            IF(IPLACE .LT. NOP .AND. NEWORB .LT. IIOP(IPLACE+1)
     &      .OR.
     &      IPLACE .EQ. NOP .AND. NEWORB .LE. MXMPTY ) THEN
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
           ELSE IF (  IPLACE .NE. NOP ) THEN
              IF(IPLACE.EQ.1) THEN
                NEWORB = 1 - 1
              ELSE
                NEWORB = IIOP(IPLACE-1)
              END IF
 671          CONTINUE
              NEWORB = NEWORB + 1
              IF(IIOC(NEWORB) .NE. 0 .AND. NEWORB.LT.MXMPTY ) GOTO 671
              IIOP(IPLACE) = NEWORB
              IIOC(NEWORB) = 1
              GOTO 700
           ELSE
C. No more active configurations , so
              IF( NCL .NE. 0 ) GOTO 2000
              IF( NCL .EQ. 0 ) GOTO 5001
           END IF
  701      CONTINUE
           IFRSTO = 0

          IF( NTEST .GE.1500) THEN
            WRITE(6,*) ' Next active configuration '
            CALL IWRTMA(IIOP,1,NOP,1,NOP)
          END IF
C        RAS  CONSTRAINTS
         IOKAY=1
           IEL1 = IEL1C
           IEL3 = IEL3C
C..        CHECK RAS1 and RAS3
           DO  40 IOP = 1,NOP
             IORB = IIOP(IOP)
             IF(IORB1F.LE.IORB .AND. IORB.LE.IORB1L ) THEN
                IEL1 = IEL1 + 1
             ELSE IF (IORB3F .LE. IORB .AND. IORB .LE. IORB3L)THEN
                IEL3 = IEL3 + 1
             END IF
   40      CONTINUE
*. Faster routine for RAS 3, added august 1990
           IR3CHK  = 1
           IF(IEL3.GT.NEL3MX.AND.IR3CHK.EQ.1) THEN
*. Number of electrons in substring
             IFSTR3 = 0
             DO 5607 IOP = 1, NOP
               IF(IIOP(IOP).GE.IORB3F) THEN
                IFSTR3 = IOP
                GOTO 5608
               END IF
 5607        CONTINUE
 5608        CONTINUE
             IF(IFSTR3.NE.NOP) THEN

*. Lowest possible string with NOP electrons
               DO 5610 K = 1, IFSTR3
                 IIOC(IIOP(K)) = 0
 5610          CONTINUE
*
               KEL = 0
               KORB = 0
 5630          CONTINUE
                 KORB = KORB + 1
                 IF(IIOC(KORB).NE.2) THEN
                   KEL = KEL + 1
                   IIOC(KORB) = 1
                   IIOP(KEL) = KORB
                 END IF
               IF(KEL.NE.IFSTR3) GOTO  5630
               IPLACE = IFSTR3
               GOTO 700
             END IF
           END IF
           IF( IEL1 .LT. NEL1MN .OR. IEL3 .GT. NEL3MX ) GOTO  999

C. Spatial symmetry
         ISYM = ISYMCN_MCLR(IICL,IIOP,NCL,NOP)
         IF(ISYM.EQ.IREFSM) THEN
           IF(NTEST.GE.100)
     &     WRITE(6,1120) ( IIOC(I),I = 1,NORB )
 1120      FORMAT(1H0,'  configuration included ',15I3)
           JCONF=JCONF+1
C
           DO 60 ICL = 1, NCL
             ICONF(ICFREE-1+ICL) = IICL(ICL)
   60      CONTINUE
           DO 61 IOP = 1, NOP
             ICONF(ICFREE-1+NCL+IOP) = IIOP(IOP)
   61      CONTINUE
           ICFREE = ICFREE + NOP + NCL
         END IF

C
C** LOOP OVER active configurations end
C
  999   CONTINUE
        IF( NOP .EQ. 0 .AND. NCL .EQ. 0 ) GOTO 5001
        IF( NOP .EQ. 0 ) GOTO 2000
        GOTO 1000
 2001 CONTINUE
 5000 CONTINUE
 5001 CONTINUE

*
      IF( NTEST .GE. 100)THEN
        WRITE(6,'(/A,I3)') '  Configurations of symmetry ', IREFSM
        WRITE(6,*)        ' ================================='
        IBAS = 0
        DO 1200 IOPEN = MINOP,MAXOP
          ITYPE = IOPEN - MINOP + 1
          ICL = (NEL-IOPEN)/2
          IOC = IOPEN + ICL
          LICONF = NCNFTP(ITYPE)
          WRITE(6,'(/A,2I3)')
     &    '  Type with number of closed and open orbitals ',ICL,IOPEN
          WRITE(6,'(A,I7)')
     &    '  Number of configurations of this type',LICONF
          DO 1180 JCONF = 1,LICONF
            WRITE(6,'(3X,20I3)') (ICONF(IBAS+IORB),IORB=1,IOC)
            IBAS = IBAS + IOC
 1180     CONTINUE
 1200   CONTINUE
      END IF

C
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE ICISPS(IPRNT)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1991
* Last revision April 1991
      IMPLICIT REAL*8(A-H,O-Z)
*
* ===================
*.Input common blocks
* ===================
*
#include "detdim.fh"
#include "cstate_mclr.fh"
#include "strinp_mclr.fh"
#include "strbas_mclr.fh"
#include "csm.fh"
#include "stinf_mclr.fh"
*
#include "WrkSpc.fh"
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp_mclr.fh"
*
*
*.Local memory
      Call GetMem('KLBLTP','ALLO','INTE',KLBLTP,NSMST)
*     IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
*    &Call GetMem('KLCVST','INTE','ALLO',KLCVST,NSMST)
      klcvst=1

*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
*     IF(IDC.EQ.3.OR.IDC.EQ.4)
*    &CALL SIGVST(iWORK(KLCVST),NSMST)

*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
C            SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
        CALL SMOST_MCLR(NSMST,NSMCI,MXPCSM,ISMOST)

      MXSB = 0
      MXSOOB = 0
      MXCEXP = 0
      DO 100 ICI = 1, NICISP
      DO  50 ISYM = 1, NSMCI
        IATP = IASTFI(ICI)
        IBTP = IBSTFI(ICI)
CMS        write(6,*) ' NRASDT : ICI IATP IBTP ',ICI,IATP,IBTP
        IF(NAELCI(ICI).EQ.NBELCI(ICI)) THEN
          IIDC = IDC
        ELSE
          IIDC = 1
        END IF
        IF(IACTI(ICI).EQ.1) THEN
          CALL ZBLTP(ISMOST(1,ISYM),NSMST,IIDC,iWORK(KLBLTP),
     &               iWORK(KLCVST))
          CALL NRASDT(MNR1IC(ICI),MXR1IC(ICI),MNR3IC(ICI),MXR3IC(ICI),
     &         ISYM,NSMST,NOCTYP(IATP),NOCTYP(IBTP),iWORK(KEL1(IATP)),
     &         iWORK(KEL1(IBTP)),iWORK(KNSTSO(IATP)),
     &         iWORK(KNSTSO(IBTP)),
     &         iWORK(KEL3(IATP)),iWORK(KEL3(IBTP)),
     &         NCOMB,XNCOMB,MXS,MXSOO,iWORK(KLBLTP))
          XISPSM(ISYM,ICI) = XNCOMB
          MXSOOB = MAX(MXSOOB,MXSOO)
          MXSB = MAX(MXSB,MXS)
          MXCEXP = MAX(MXCEXP,NCOMB)
        ELSE
          XISPSM(ISYM,ICI) = 0.0D0
        END IF
   50 CONTINUE
  100 CONTINUE
      Call GetMem('KLBLTP','Free','INTE',KLBLTP,NSMST)
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      If (ntest.ne.0) Then
      WRITE(6,*)
     &' Number of internal combinations per symmetry '
      WRITE(6,*)
     & ' =========================================== '
      IF(NTEST.EQ.0) THEN
        MX = 1
      ELSE
        MX = NICISP
      END IF
*
      DO 200 ICI = 1, MX
        IF(IACTI(ICI).EQ.1) THEN
          WRITE(6,*) ' Internal CI space ', ICI
          CALL WRTMAT(XISPSM(1,ICI),1,NSMCI,1,NSMCI)
        END IF
  200 CONTINUE
      WRITE(6,*) ' Largest CI space                 ',MXCEXP
      WRITE(6,*) ' Largest symmetry block           ',MXSB
      WRITE(6,*) ' Largest Symmetry-type-type block ',MXSOOB
      End If
*
      RETURN
      END
*
*----------------------------------------------------------------------
*
*
*----------------------------------------------------------------------
*
      SUBROUTINE NRASDT(MNRS1,MXRS1,MNRS3,MXRS3,ITOTSM,
     &                  NSMST,NOCTPA,NOCTPB,IEL1A,IEL1B,
     &                  NSSOA,NSSOB,
     &                  IEL3A,IEL3B,NCOMB,XNCOMB,MXSB,MXSOOB,
     &                  IBLTP)
*
* Number of combimations with symmetry ITOTSM and
*       MNRS1 - MXRS1 elecs in RAS1
*       MNRS3 - MXRS3 elecs in RAS3
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
*
* Updated with IBLTP, Summer of 93
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IEL1A(*),IEL3A(*)
      DIMENSION IEL1B(*),IEL3B(*)
      DIMENSION NSSOA(NOCTPA,*),NSSOB(NOCTPB,*)
      DIMENSION IBLTP(*)
*
      MXSB = 0
      MXSOOB = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        CALL SYMCOM_MCLR(2,1,IASM,IBSM,ITOTSM)
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
             MXBTP = IATP
           ELSE
             MXBTP = NOCTPB
           END IF
           DO 100 IBTP = 1, MXBTP
             IEL1 = IEL1A(IATP)+IEL1B(IBTP)
             IEL3 = IEL3A(IATP)+IEL3B(IBTP)
             IF(IEL1.GE.MNRS1.AND.IEL1.LE.MXRS1.AND.
     &       IEL3.GE.MNRS3.AND.IEL3.LE.MXRS3 ) THEN
*. Size of unpacked block
               LTTSUP =  NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
*. Size of packed block
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 LTTSBL = NSSOA(IATP,IASM)*NSSOB(IBTP,IBSM)
               ELSE
                 LTTSBL = NSSOA(IATP,IASM)*(NSSOA(IATP,IASM)+1)/2
               END IF
               NCOMB = NCOMB + LTTSBL
               LSB = LSB + LTTSUP
               MXSOOB = MAX(MXSOOB,LTTSUP)
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 XNCOMB = XNCOMB +
     &         DBLE(NSSOA(IATP,IASM))*DBLE(NSSOB(IBTP,IBSM))
               ELSE
                 XNCOMB = XNCOMB +
     &           DBLE(NSSOA(IATP,IASM))*
     &           (DBLE(NSSOB(IBTP,IBSM))+1.0D0)/2.0D0
               END IF
             END IF
  100      CONTINUE
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
  300 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NCOMB and XNCOMB ', NCOMB,XNCOMB
      END IF
*
      RETURN
      END
