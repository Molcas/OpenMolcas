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
* Copyright (C) 1994,1996,2014, Per Ake Malmqvist                      *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE NEWFOCK(FIFA,NFIFA,CMO,NCMO)
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: usual
      use caspt2_global, only: DREF
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
#include "caspt2.fh"
      INTEGER NFIFA, NCMO
      REAL*8 FIFA(NFIFA),CMO(NCMO)

      REAL*8 D,DDVX,E,EIGVAL
      INTEGER LINT,LSC,LSC1,LSC2,LSCR,
     &        LEV1,LEV2,LEIG,LXAI,LXPQ,LXQP
      INTEGER IFGFOCK
      INTEGER I,J
      INTEGER II,IP,IQ,IR,IS,IV,IX
      INTEGER ITOT,IA,IATOT,IT,ITTOT,ITABS,IU,IUTOT,IUABS,ITU
      INTEGER NA,NA2,NA3,MA,MI,MTRES,N3,NI,NIA,NINT,NO,
     &        NAS,NASQES,NASQT,NATR,NATRES,NOSQES,NOTRES,
     &        NS,NSCR,NSCR1,NSCR2,NSCR3,NSQES,NTRES
      INTEGER ISC,ISTLT,KFIFA
      INTEGER ID,IDDVX,IDREF,IDTT,IDTU,IDUT
      INTEGER IEPS,IEPSA,IEPSI,IEPSE
      INTEGER ISYM,ISYMPQ,ISYMRS
      REAL*8 VAL,VALTU,VALUT,X

      Real*8, allocatable:: INT(:), DSQ(:), DD(:), DDTR(:),
     &                      TWOMDSQ(:), XMAT(:), SC(:)
c Purpose: Modify the standard fock matrix for experimental
c purposes. The string variable FOCKTYPE (character*8) has a
c keyword value given as input. The experimental user modifies
c this routine to suit his purposes, and links with the rest
c of the program. The Fock matrix is given as call parameter
c and returned after modification.
c To define the modified Fock matrix, a number of arrays on
c  LUONE may be useful. In addition, the active 1- and 2-
c electron density matrices, and the inactive Fock matrix
c FIMO, are available in workspace at DREF,
c PREF,FIMO, and FIFA.
c Two-electron integrals involving non-frozen, non-deleted
c orbitals, at most two secondary, are available from
c subroutines COUL and EXCH (See).
c Meaningful modifications must define a Fock operator
c which is invariant to transformations among inactives,
c among actives, and among virtuals.
c Coded 94-01-31 by Malmqvist, for CASPT2 MOLCAS-3.
c Modif 96-10-06 by Malmqvist, restructured, options added.
c Modif 14-03-19 by Malmqvist, restructured, options removed.

      ! I never meant to cause you any sorrow
      IF(FOCKTYPE.EQ.'STANDARD') RETURN

* Options MC and MC2 removed, PAM March 2014.
*CPAM96 The option FOCKTYPE='MC      ' added 961006. This option will
*C replace the active/active block with the MCSCF Fock matrix
*C while zeroing any non-diagonal blocks. This option is usually
*C quite ridiculous, but it can be used in very particular cases
*C when all active orbitals are singly occupied.
*CPAM96 The option FOCKTYPE='MC2     ' added 961006. Similar to the
*C above, but using as  active/active block the matrix
*C    D**(-1/2) FMC D**(-1/2)


CPAM A very long IF-block is replaced by a forward GOTO for clarity:
CPAM      IF((FOCKTYPE.EQ.'G1      '.OR.FOCKTYPE.EQ.'G2      '.OR.
CPAM     &   FOCKTYPE.EQ.'G3      ').AND.NASHT.GT.0) THEN
      IFGFOCK=0
      IF(FOCKTYPE.EQ.'G1      ') IFGFOCK=1
      IF(FOCKTYPE.EQ.'G2      ') IFGFOCK=1
      IF(FOCKTYPE.EQ.'G3      ') IFGFOCK=1

      IF(IFGFOCK.EQ.0) GOTO 300
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(6,*)' THE FOCK MATRIX IS MODIFIED BY KEYWORD'//
     &               ' FOCKTYPE=',FOCKTYPE
      END IF
      IF(NASHT.LE.0) GOTO 300
c
c       Determine sizes of areas for memory allocation
        NASQT=0
        NATR=0
        DO 5 ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NS=NSSH(ISYM)
          NO=NORB(ISYM)
          NASQT=NASQT+NA**2
          NATR=NATR+(NA*(NA+1))/2
    5   CONTINUE
        NINT=NOMX**2
        NSCR1=NAMX*MAX(2*NAMX,NIMX,NSMX)
        NSCR2=3*NAMX*(NAMX+1)
        NSCR3=NAMX*(3*NAMX+1)
        NSCR=NSCR1
        IF(FOCKTYPE.EQ.'G2      ') NSCR=NSCR2
        IF(FOCKTYPE.EQ.'G3      ') NSCR=NSCR3
c
c Allocate memory: Integral buffer and scratch array:
        CALL mma_allocate(INT,2*NINT,LABEL='INT')
        LINT=1
        LSCR=LINT+NINT
        CALL mma_allocate(SC,NSCR,Label='SC')
        LSC=1
c Form symmetry-packed squares of density matrix DSQ, and
c similarly (2I-DSQ)
        CALL mma_allocate(DSQ,NASQT,Label='DSQ')
        CALL mma_allocate(TWOMDSQ,NASQT,LABEL='TWOMDSQ')
c Symmetry-packed triangles of D*(2I-D)
        CALL mma_allocate(DDTR,NATR,Label='DDTR')
c Temporary use of single square symmetry-block:
        CALL mma_allocate(DD,NAMX**2,Label='DD')
C The exchange matrix, A(pq)=sum over rs of (ps,rq)*DD(rs)
        CALL mma_allocate(XMAT,NOSQT,Label='XMAT')
        NSQES=0
        DO ISYM=1,NSYM
          NA=NASH(ISYM)
          DO IT=1,NA
            ITABS=IT+NAES(ISYM)
            DO IU=1,NA
              IUABS=IU+NAES(ISYM)
              IDREF=(ITABS*(ITABS-1))/2+IUABS
              D=DREF(IDREF)
              IDTU=NSQES+IT+NA*(IU-1)
              IDUT=NSQES+IU+NA*(IT-1)
              DSQ(IDTU)=D
              DSQ(IDUT)=D
            END DO
          END DO
          DO I=1,NA*NA
            IDTU=NSQES+I
            TWOMDSQ(IDTU)=-DSQ(IDTU)
          END DO
          DO I=1,NA*NA,(NA+1)
            IDTT=NSQES+I
            TWOMDSQ(IDTT)=2.0D0-DSQ(IDTT)
          END DO
          NSQES=NSQES+NA**2
        END DO
c
c Create the matrix DDTR =D(2I-D) (triangular symmetry blocks)
C Use also temporary DD, single symmetry blocks of D*(2I-D):
        NTRES=1
        NSQES=1
        DO ISYM=1,NSYM
          NA=NASH(ISYM)
          IF(NA.GT.0) THEN
            N3=(NA*(NA+1))/2
            CALL DGEMM_('N','N',
     &                  NA,NA,NA,
     &                  1.0d0,DSQ(NSQES),NA,
     &                  TWOMDSQ(NSQES),NA,
     &                  0.0d0,DD,NA)
            CALL TRIANG(NA,DD)
            CALL DCOPY_(N3,DD,1,DDTR(NTRES),1)
            NTRES=NTRES+N3
            NSQES=NSQES+NA**2
          END IF
        END DO

C Calculation of the exchange matrix, A(pq)=sum over rs of (ps,rq)*DD(rs)
        XMAT(:)=0.0D0
        IF (IfChol) THEN
          CALL Cho_Amatrix(XMAT,CMO,NCMO,DDTR,NATR)
        ELSE
          NOSQES=0
          DO ISYMPQ=1,NSYM
            NI=NISH(ISYMPQ)
            NA=NASH(ISYMPQ)
            NS=NSSH(ISYMPQ)
            NO=NORB(ISYMPQ)
            IF(NO.GT.0) THEN
              MTRES=0
              DO ISYMRS=1,NSYM
                MI=NISH(ISYMRS)
                MA=NASH(ISYMRS)
                DO IV=1,MA
                  IR=IV+MI
                  DO IX=1,IV
                    IS=IX+MI
                    CALL EXCH(ISYMPQ,ISYMRS,ISYMPQ,ISYMRS,IR,IS,
     &                        INT,INT(LSCR))
                    IDDVX=MTRES+(IV*(IV-1))/2+IX
                    DDVX=DDTR(IDDVX)
                    IF(IR.EQ.IS) DDVX=0.5D0*DDVX
                    CALL DAXPY_(NO**2,DDVX,INT,1,
     &                                    XMAT(1+NOSQES),1)
                  END DO
                END DO
                MTRES=MTRES+(MA*(MA+1))/2
              END DO
              DO IP=2,NO
                DO IQ=1,IP-1
                  LXPQ=NOSQES+IP+NO*(IQ-1)
                  LXQP=NOSQES+IQ+NO*(IP-1)
                  VAL=0.5D0*(XMAT(LXPQ)+XMAT(LXQP))
                  XMAT(LXPQ)=VAL
                  XMAT(LXQP)=VAL
                END DO
              END DO
              NOSQES=NOSQES+NO**2
            END IF
          END DO
        END IF
c
c Determine the correction to the Fock matrix
        IF(FOCKTYPE.EQ.'G1      ') THEN
C Focktype=g1 case. A very long IF block.
          NOSQES=0
          NOTRES=0
          NASQES=0
          DO 30 ISYMPQ=1,NSYM
            NI=NISH(ISYMPQ)
            NA=NASH(ISYMPQ)
            NS=NSSH(ISYMPQ)
            NO=NORB(ISYMPQ)
            NIA=NI+NA
            NAS=NA+NS
            IF(NIA*NAS.LE.0) GO TO 31
c
c the active-inactive block
            IF(NA*NI.GT.0) THEN
              CALL DGEMM_('N','N',
     &                    NA,NI,NA,
     &                    1.0d0,TWOMDSQ(1+NASQES),NA,
     &                    XMAT(1+NOSQES+NI),NO,
     &                    0.0d0,SC,NA)
              DO IT=1,NA
                ITTOT=NI+IT
                DO II=1,NI
                  KFIFA=NOTRES+(ITTOT*(ITTOT-1))/2+II
                  ISC=IT+NA*(II-1)
                  FIFA(KFIFA)=FIFA(KFIFA)-SC(ISC)
                END DO
              END DO
            ENDIF
c
c the secondary-inactive block
            IF(NS*NI.GT.0) THEN
              DO IA=1,NS
                IATOT=NI+NA+IA
                DO II=1,NI
                  KFIFA=NOTRES+(IATOT*(IATOT-1))/2+II
                  LXAI=NOSQES+IATOT+NO*(II-1)
                  FIFA(KFIFA)=FIFA(KFIFA)-2.0D0*XMAT(LXAI)
                END DO
              END DO
            ENDIF
c
c the active-active block
            IF(NA.GT.0) THEN
              IX=NOSQES+NI+1+NO*NI
              CALL DGEMM_('N','N',
     &                    NA,NA,NA,
     &                    1.0d0,DSQ(1+NASQES),NA,
     &                    XMAT(IX),NO,
     &                    0.0d0,SC(LSC+NA*NA),NA)
              CALL DGEMM_('N','N',
     &                    NA,NA,NA,
     &                    1.0d0,SC(LSC+NA*NA),NA,
     &                    TWOMDSQ(1+NASQES),NA,
     &                    0.0d0,SC,NA)
              DO IT=1,NA
                ITTOT=NI+IT
                DO IU=1,IT
                  IUTOT=NI+IU
                  KFIFA=NOTRES+(ITTOT*(ITTOT-1))/2+IUTOT
                  VALTU=SC(IT+NA*(IU-1))
                  VALUT=SC(IU+NA*(IT-1))
                  FIFA(KFIFA)=FIFA(KFIFA)-0.5D0*(VALTU+VALUT)
                END DO
              END DO
            ENDIF
c
c the secondary-active block
            IF(NS*NA.GT.0) THEN
              IX=NOSQES+NI+NA+1+NO*NI
              CALL DGEMM_('N','N',
     &                    NS,NA,NA,
     &                    1.0d0,XMAT(IX),NO,
     &                    DSQ(1+NASQES),NA,
     &                    0.0d0,SC,NS)
              DO IA=1,NS
                IATOT=NI+NA+IA
                DO IT=1,NA
                  ITTOT=NI+IT
                  KFIFA=NOTRES+(IATOT*(IATOT-1))/2+ITTOT
                  FIFA(KFIFA)=FIFA(KFIFA)-SC(IA+NS*(IT-1))
                END DO
              END DO
            ENDIF
c
   31       CONTINUE
            NOSQES=NOSQES+NO**2
            NOTRES=NOTRES+(NO*(NO+1))/2
            NASQES=NASQES+NA**2
   30     CONTINUE
C Focktype=g1 case ends.
        ELSE
C Focktype=g2 or g3
          NOSQES=0
          NOTRES=0
          NASQES=0
          NATRES=0
          DO 130 ISYMPQ=1,NSYM
            NI=NISH(ISYMPQ)
            NA=NASH(ISYMPQ)
            NO=NORB(ISYMPQ)
            NIA=NI+NA
            NA2=NA*NA
            NA3=(NA+NA2)/2
            IF(NA.LE.0) GO TO 131
c
c Determine the matrix blocks of the correction to
c the Fock matrix and add them to the Fock matrix
c
C Form the selection matrix as a temporary square matrix.
C Compute it by spectral resolution.
C First, form a copy of the triangular D(2I-D) matrix block,
C and diagonalize it. The DDTR copy at LSC:
            CALL DCOPY_(NA3,DDTR(1+NATRES),1,SC,1)
C A unit matrix at LEV1, to become eigenvectors:
            LEV1=LSC+NA2
            CALL DCOPY_(NA2,[0.0D0],0,SC(LEV1),1)
            CALL DCOPY_(NA, [1.0D0],0,SC(LEV1),NA+1)
C A call to NIDiag diagonalizes the triangular matrix:
            CALL NIDiag(SC,SC(LEV1),NA,NA)
            CALL JACORD(SC,SC(LEV1),NA,NA)
C Make a copy of the eigenvector matrix:
            LEV2=LEV1+NA2
            CALL DCOPY_(NA2,SC(LEV1),1,SC(LEV2),1)
C Put eigenvalues at LEIG:
            LEIG=LEV2+NA2
            CALL VEIG(NA,SC,SC(LEIG))
C Now scale the second array of eigenvectors with any required
C function of the eigenvalues:
            DO J=1,NA
              EIGVAL=SC(LEIG-1+J)
              IF(FOCKTYPE.EQ.'G2      ') THEN
                X=SQRT(MAX(0.0D0,EIGVAL))
              ELSE
                X=EIGVAL
              END IF
              DO I=1,NA
                SC(LEV2-1+I+NA*(J-1))=X*SC(LEV2-1+I+NA*(J-1))
              END DO
            END DO
C Now the selection matrix can be formed, at LSC:
            CALL DGEMM_('N','T',
     &                  NA,NA,NA,
     &                  1.0d0,SC(LEV1),NA,
     &                  SC(LEV2),NA,
     &                  0.0d0,SC(LSC),NA)
C Obviously, the FOCKTYPE=G3 case can be obtained by just
C squaring the DDTR block into SC.

C Focktype=g2 or g3
            IX=NOSQES+NI+1+NO*NI
            LSC1=LSC+NA2
            LSC2=LSC1+NA2
            CALL DGEMM_('N','N',
     &                  NA,NA,NA,
     &                  1.0d0,SC(LSC),NA,
     &                  XMAT(IX),NO,
     &                  0.0d0,SC(LSC1),NA)
            CALL DGEMM_('N','N',
     &                  NA,NA,NA,
     &                  1.0d0,SC(LSC1),NA,
     &                  SC(LSC),NA,
     &                  0.0d0,SC(LSC2),NA)
            DO IT=1,NA
              ITTOT=NI+IT
              DO IU=1,IT
                IUTOT=NI+IU
                KFIFA=NOTRES+(ITTOT*(ITTOT-1))/2+IUTOT
                ITU=IT+NA*(IU-1)
                FIFA(KFIFA)=FIFA(KFIFA)-SC(LSC2-1+ITU)
              END DO
            END DO
c
  131       CONTINUE
            NOSQES=NOSQES+NO**2
            NOTRES=NOTRES+(NO*(NO+1))/2
            NASQES=NASQES+NA**2
            NATRES=NATRES+(NA*(NA+1))/2
  130     CONTINUE
        ENDIF
c
c
        CALL mma_deallocate(SC)
        CALL mma_deallocate(INT)
        CALL mma_deallocate(DSQ)
        CALL mma_deallocate(TWOMDSQ)
        CALL mma_deallocate(DD)
        CALL mma_deallocate(DDTR)
        CALL mma_deallocate(XMAT)
 300  CONTINUE
c
c     Orbital energies, EPS, EPSI,EPSA,EPSE:
      IEPS=0
      IEPSI=0
      IEPSA=0
      IEPSE=0
      ISTLT=0
      DO 400 ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        DO 401 I=1,NI
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSI=IEPSI+1
          EPSI(IEPSI)=E
  401   CONTINUE
        DO 402 I=NI+1,NI+NA
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSA=IEPSA+1
          EPSA(IEPSA)=E
  402   CONTINUE
        DO 403 I=NI+NA+1,NO
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSE=IEPSE+1
          EPSE(IEPSE)=E
  403   CONTINUE
        ISTLT=ISTLT+(NO*(NO+1))/2
  400 CONTINUE
c
c     EASUM = contract EPSA with diagonal of active dens.
      EASUM=0.0D+00
      DO 410 ISYM=1,NSYM
        NA=NASH(ISYM)
        DO 411 I=1,NA
          ITOT=NAES(ISYM)+I
          ID=(ITOT*(ITOT+1))/2
          EASUM=EASUM+EPSA(ITOT)*DREF(ID)
  411   CONTINUE
  410 CONTINUE
c
      END SUBROUTINE NEWFOCK
