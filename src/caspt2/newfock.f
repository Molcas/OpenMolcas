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
      SUBROUTINE NEWFOCK(FIFA)
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      REAL*8 FIFA(NFIFA)

      REAL*8 D,DDVX,E,EIGVAL
      INTEGER LDD,LDDTR,LDSQ,L2MDSQ,LINT,LSC,LSC1,LSC2,LSCR,
     &        LEV1,LEV2,LEIG,LXAI,LXMAT,LXPQ,LXQP
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
c Purpose: Modify the standard fock matrix for experimental
c purposes. The string variable FOCKTYPE (character*8) has a
c keyword value given as input. The experimental user modifies
c this routine to suit his purposes, and links with the rest
c of the program. The Fock matrix is given as call parameter
c and returned after modification.
c To define the modified Fock matrix, a number of arrays on
c  LUONE may be useful. In addition, the active 1- and 2-
c electron density matrices, and the inactive Fock matrix
c FIMO, are awailable in workspace at WORK(LDREF),
c WORK(LPREF),WORK(LFIMO), and WORK(LFIFA).
c Two-electron integrals involving non-frozen, non-deleted
c orbitals, at most two secondary, are available from
c subroutines COUL and EXCH (See).
c Meaningful modifications must define a Fock operator
c which is invariant to transformations among inactives,
c among actives, and among virtuals.
c Coded 94-01-31 by Malmqvist, for CASPT2 MOLCAS-3.
c Modif 96-10-06 by Malmqvist, restructured, options added.
c Modif 14-03-19 by Malmqvist, restructured, options removed.

      IF(FOCKTYPE.EQ.'STANDARD') RETURN
      IF(IFCHOL) THEN
        WRITE(6,*)' Subroutine NEWFOCK is presently unable to use'
        WRITE(6,*)' Cholesky vectors. The FOCKTYPE variable is now'
        WRITE(6,*)' changed to STANDARD, and NEWFOCK returns without'
        WRITE(6,*)' action. This will be fixed as soon as possible.'
        FOCKTYPE='STANDARD'
      END IF

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
        CALL GETMEM('FINT','ALLO','REAL',LINT,2*NINT)
        LSCR=LINT+NINT
        CALL GETMEM('FSCR','ALLO','REAL',LSC,NSCR)
c Form symmetry-packed squares of density matrix DSQ, and
c similarly (2I-DSQ)
        CALL GETMEM('DSQ','ALLO','REAL',LDSQ,NASQT)
        CALL GETMEM('2MDSQ','ALLO','REAL',L2MDSQ,NASQT)
c Symmetry-packed triangles of D*(2I-D)
        CALL GETMEM('DDTR','ALLO','REAL',LDDTR,NATR)
c Temporary use of single square symmetry-block:
        CALL GETMEM('DD','ALLO','REAL',LDD,NAMX**2)
C The exchange matrix, A(pq)=sum over rs of (ps,rq)*DD(rs)
        CALL GETMEM('XMAT','ALLO','REAL',LXMAT,NOSQT)
        NSQES=0
        DO ISYM=1,NSYM
          NA=NASH(ISYM)
          DO IT=1,NA
            ITABS=IT+NAES(ISYM)
            DO IU=1,NA
              IUABS=IU+NAES(ISYM)
              IDREF=(ITABS*(ITABS-1))/2+IUABS
              D=WORK(LDREF-1+IDREF)
              IDTU=NSQES+IT+NA*(IU-1)
              IDUT=NSQES+IU+NA*(IT-1)
              WORK(LDSQ-1+IDTU)=D
              WORK(LDSQ-1+IDUT)=D
            END DO
          END DO
          DO I=1,NA*NA
            IDTU=NSQES+I
            WORK(L2MDSQ-1+IDTU)=-WORK(LDSQ-1+IDTU)
          END DO
          DO I=1,NA*NA,(NA+1)
            IDTT=NSQES+I
            WORK(L2MDSQ-1+IDTT)=2.0D0-WORK(LDSQ-1+IDTT)
          END DO
          NSQES=NSQES+NA**2
        END DO
c
c Create the matrix DDTR =D(2I-D) (triangular symmetry blocks)
C Use also temporary DD, single symmetry blocks of D*(2I-D):
        NTRES=0
        NSQES=0
        DO ISYM=1,NSYM
          NA=NASH(ISYM)
          IF(NA.GT.0) THEN
            N3=(NA*(NA+1))/2
            CALL DGEMM_('N','N',
     &                  NA,NA,NA,
     &                  1.0d0,WORK(LDSQ+NSQES),NA,
     &                  WORK(L2MDSQ+NSQES),NA,
     &                  0.0d0,WORK(LDD),NA)
            CALL TRIANG(NA,WORK(LDD))
            CALL DCOPY_(N3,WORK(LDD),1,WORK(LDDTR+NTRES),1)
            NTRES=NTRES+N3
            NSQES=NSQES+NA**2
          END IF
        END DO

C Calculation of the exchange matrix, A(pq)=sum over rs of (ps,rq)*DD(rs)
        CALL DCOPY_(NOSQT,[0.0D0],0,WORK(LXMAT),1)
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
     &                      WORK(LINT),WORK(LSCR))
                  IDDVX=MTRES+(IV*(IV-1))/2+IX
                  DDVX=WORK(LDDTR-1+IDDVX)
                  IF(IR.EQ.IS) DDVX=0.5D0*DDVX
                  CALL DAXPY_(NO**2,DDVX,WORK(LINT),1,
     &                                  WORK(LXMAT+NOSQES),1)
                END DO
              END DO
              MTRES=MTRES+(MA*(MA+1))/2
            END DO
            DO IP=2,NO
              DO IQ=1,IP-1
                LXPQ=LXMAT+NOSQES-1+IP+NO*(IQ-1)
                LXQP=LXMAT+NOSQES-1+IQ+NO*(IP-1)
                VAL=0.5D0*(WORK(LXPQ)+WORK(LXQP))
                WORK(LXPQ)=VAL
                WORK(LXQP)=VAL
              END DO
            END DO
            NOSQES=NOSQES+NO**2
          END IF
        END DO
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
     &                    1.0d0,WORK(L2MDSQ+NASQES),NA,
     &                    WORK(LXMAT+NOSQES+NI),NO,
     &                    0.0d0,WORK(LSC),NA)
              DO IT=1,NA
                ITTOT=NI+IT
                DO II=1,NI
                  KFIFA=NOTRES+(ITTOT*(ITTOT-1))/2+II
                  ISC=IT+NA*(II-1)
                  FIFA(KFIFA)=FIFA(KFIFA)-WORK(LSC-1+ISC)
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
                  LXAI=LXMAT+NOSQES-1+IATOT+NO*(II-1)
                  FIFA(KFIFA)=FIFA(KFIFA)-2.0D0*WORK(LXAI)
                END DO
              END DO
            ENDIF
c
c the active-active block
            IF(NA.GT.0) THEN
              IX=NOSQES+NI+1+NO*NI
              CALL DGEMM_('N','N',
     &                    NA,NA,NA,
     &                    1.0d0,WORK(LDSQ+NASQES),NA,
     &                    WORK(LXMAT-1+IX),NO,
     &                    0.0d0,WORK(LSC+NA*NA),NA)
              CALL DGEMM_('N','N',
     &                    NA,NA,NA,
     &                    1.0d0,WORK(LSC+NA*NA),NA,
     &                    WORK(L2MDSQ+NASQES),NA,
     &                    0.0d0,WORK(LSC),NA)
              DO IT=1,NA
                ITTOT=NI+IT
                DO IU=1,IT
                  IUTOT=NI+IU
                  KFIFA=NOTRES+(ITTOT*(ITTOT-1))/2+IUTOT
                  VALTU=WORK(LSC-1+IT+NA*(IU-1))
                  VALUT=WORK(LSC-1+IU+NA*(IT-1))
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
     &                    1.0d0,WORK(LXMAT-1+IX),NO,
     &                    WORK(LDSQ+NASQES),NA,
     &                    0.0d0,WORK(LSC),NS)
              DO IA=1,NS
                IATOT=NI+NA+IA
                DO IT=1,NA
                  ITTOT=NI+IT
                  KFIFA=NOTRES+(IATOT*(IATOT-1))/2+ITTOT
                  FIFA(KFIFA)=FIFA(KFIFA)-WORK(LSC-1+IA+NS*(IT-1))
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
            CALL DCOPY_(NA3,WORK(LDDTR+NATRES),1,WORK(LSC),1)
C A unit matrix at LEV1, to become eigenvectors:
            LEV1=LSC+NA2
            CALL DCOPY_(NA2,[0.0D0],0,WORK(LEV1),1)
            CALL DCOPY_(NA, [1.0D0],0,WORK(LEV1),NA+1)
C A call to NIDiag diagonalizes the triangular matrix:
            CALL NIDiag(WORK(LSC),WORK(LEV1),NA,NA,0)
            CALL JACORD(WORK(LSC),WORK(LEV1),NA,NA)
C Make a copy of the eigenvector matrix:
            LEV2=LEV1+NA2
            CALL DCOPY_(NA2,WORK(LEV1),1,WORK(LEV2),1)
C Put eigenvalues at LEIG:
            LEIG=LEV2+NA2
            CALL VEIG(NA,WORK(LSC),WORK(LEIG))
C Now scale the second array of eigenvectors with any required
C function of the eigenvalues:
            DO J=1,NA
              EIGVAL=WORK(LEIG-1+J)
              IF(FOCKTYPE.EQ.'G2      ') THEN
                X=SQRT(MAX(0.0D0,EIGVAL))
              ELSE
                X=EIGVAL
              END IF
              DO I=1,NA
                WORK(LEV2-1+I+NA*(J-1))=X*WORK(LEV2-1+I+NA*(J-1))
              END DO
            END DO
C Now the selection matrix can be formed, at LSC:
            CALL DGEMM_('N','T',
     &                  NA,NA,NA,
     &                  1.0d0,WORK(LEV1),NA,
     &                  WORK(LEV2),NA,
     &                  0.0d0,WORK(LSC),NA)
C Obviously, the FOCKTYPE=G3 case can be obtained by just
C squaring the DDTR block into WORK(LSC).

C Focktype=g2 or g3
            IX=NOSQES+NI+1+NO*NI
            LSC1=LSC+NA2
            LSC2=LSC1+NA2
            CALL DGEMM_('N','N',
     &                  NA,NA,NA,
     &                  1.0d0,WORK(LSC),NA,
     &                  WORK(LXMAT-1+IX),NO,
     &                  0.0d0,WORK(LSC1),NA)
            CALL DGEMM_('N','N',
     &                  NA,NA,NA,
     &                  1.0d0,WORK(LSC1),NA,
     &                  WORK(LSC),NA,
     &                  0.0d0,WORK(LSC2),NA)
            DO IT=1,NA
              ITTOT=NI+IT
              DO IU=1,IT
                IUTOT=NI+IU
                KFIFA=NOTRES+(ITTOT*(ITTOT-1))/2+IUTOT
                ITU=IT+NA*(IU-1)
                FIFA(KFIFA)=FIFA(KFIFA)-WORK(LSC2-1+ITU)
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
        CALL GETMEM('FSCR','FREE','REAL',LSC,NSCR)
        CALL GETMEM('FINT','FREE','REAL',LINT,NINT)
        CALL GETMEM('DSQ','FREE','REAL',LDSQ,NASQT)
        CALL GETMEM('2MDSQ','FREE','REAL',L2MDSQ,NASQT)
        CALL GETMEM('DD','FREE','REAL',LDD,NAMX**2)
        CALL GETMEM('DDTR','FREE','REAL',LDDTR,NATR)
        CALL GETMEM('XMAT','FREE','REAL',LXMAT,NOSQT)
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
          EASUM=EASUM+EPSA(ITOT)*WORK(LDREF-1+ID)
  411   CONTINUE
  410 CONTINUE
c
      RETURN
      END
