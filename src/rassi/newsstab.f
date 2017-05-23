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
      INTEGER FUNCTION NEWSSTAB(LORBTAB)
      IMPLICIT NONE
      INTEGER LORBTAB
      INTEGER NTAB,ITYPE,NSSTP,NSBSTOT,KSBSMRS,KMRSSBS
      INTEGER KSBSANN,KSBSCRE
#include "symmul.fh"
#include "Morsel.fh"
#include "WrkSpc.fh"
      INTEGER I,IKMS2,IKSCR,IKSYM
      INTEGER IMS2,IPOP,ISBS,ISCR
      INTEGER ISGN,ISORB,ISPART,ISSTP
      INTEGER ISYM,J,KMS2,KOINFO,KSSTANN
      INTEGER KSSTCRE,KSSTP,KSYM,LNOSUB
      INTEGER LOINFO,LOSPN,LOSYM,LPOS,LSCR,LSCR2
      INTEGER LTAB,MORSANN,MORSCRE,MORSPOP
      INTEGER MORSSPIN,MORSSYMM,MRS,MS2
      INTEGER MXMS2,MXPOP,N
      INTEGER NEWMRS,NEWSBS,NEWSST,NO
      INTEGER NSBS,NSCR,NASPO,NASPRT
C     INTEGER NTEST,iscrmx
      INTEGER IERR,ISBS1,ISBS2,ISBS3,ISBS4,ISBS5,ISBS6,ISBS7,ISBS8
C     CHARACTER*8 STRING8
      EXTERNAL MORSANN,MORSCRE,MORSPOP,MORSSPIN,MORSSYMM

C Table type ID:
      ITYPE=19
C Pick up some data from the orbital table:
      NASPO =IWORK(LORBTAB+3)
      NSYM  =IWORK(LORBTAB+4)
      NASPRT=IWORK(LORBTAB+8)
      KOINFO=19
      LOINFO=LORBTAB-1+KOINFO
C Make temporary arrays for spin label and symmetry of each orbital
      CALL GETMEM('OrbSpn','Allo','Inte',LOSPN,NASPO)
      CALL GETMEM('OrbSym','Allo','Inte',LOSYM,NASPO)
C Make a temporary array, which gives the number of spin orbitals in
C each subpartition:
      CALL GETMEM('NOSub','Allo','Inte',LNOSUB,NASPRT)
      DO ISPART=1,NASPRT
        IWORK(LNOSUB-1+ISPART)=0
      END DO
      DO ISORB=1,NASPO
        ISYM  = IWORK(LOINFO+ 1+(ISORB-1)*8)
        IWORK(LOSYM-1+ISORB)=ISYM
        MS2   = IWORK(LOINFO+ 3+(ISORB-1)*8)
        IWORK(LOSPN-1+ISORB)=MS2
        ISPART= IWORK(LOINFO+ 6+(ISORB-1)*8)
        N=1+IWORK(LNOSUB-1+ISPART)
        IWORK(LNOSUB-1+ISPART)=N
      END DO
C We need a temporary array, NSBSSCR(NSYM,0:NPOP,-MXMS2:MXMS2,NASPRT),
C to keep the number of substrings of different kind:
      MXPOP=MORSBITS
      MXMS2=MORSBITS
      NSCR=NSYM*(1+MXPOP)*(2*MXMS2+1)*NASPRT
      CALL GETMEM('SSTScr','Allo','Inte',LSCR,NSCR)
C Addressing will be through the cumbersome formula
C ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+MS2+(2*MXMS2+1)*(ISPART-1)))
C Initialize counter of substrings:
      DO ISPART=1,NASPRT
        DO IMS2=-MXMS2,MXMS2
         DO IPOP=0,MXPOP
          DO ISYM=1,NSYM
C NSBSSCR(ISYM,IPOP,IMS2,ISPART)=0:
           ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                     (ISPART-1)))
           IWORK(LSCR-1+ISCR)=0
          END DO
         END DO
        END DO
C NSBSSCR(ISYM=1,IPOP=0,IMS2=0,ISPART)=1:
        ISCR=1+NSYM*(0+(1+MXPOP)*(MXMS2+0+(2*MXMS2+1)*(ISPART-1)))
        IWORK(LSCR-1+ISCR)=1
      END DO
C Compute number of substrings:
      ISORB=0
      DO ISPART=1,NASPRT
        NO=IWORK(LNOSUB-1+ISPART)
        DO I=1,NO
          ISORB=ISORB+1
          KSYM=IWORK(LOSYM-1+ISORB)
          KMS2=IWORK(LOSPN-1+ISORB)
          DO IPOP=I,1,-1
           DO IMS2=-IPOP,IPOP
            IKMS2=IMS2-KMS2
            IF(ABS(IKMS2).LE.IPOP-1) THEN
            DO ISYM=1,NSYM
             IKSYM=MUL(ISYM,KSYM)
             ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                      (ISPART-1)))
             N=IWORK(LSCR-1+ISCR)
             IKSCR=IKSYM+NSYM*(IPOP-1+(1+MXPOP)*(MXMS2+IKMS2+
     &                          (2*MXMS2+1)*(ISPART-1)))
             N=N+IWORK(LSCR-1+IKSCR)
             IWORK(LSCR-1+ISCR)=N
            END DO
            END IF
           END DO
          END DO
        END DO
      END DO
      NSSTP=0
      NSBSTOT=0
      DO ISPART=1,NASPRT
        NO=IWORK(LNOSUB-1+ISPART)
        DO IPOP=0,NO
         DO ISYM=1,NSYM
          DO IMS2=-NO,NO
           ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                     (ISPART-1)))
           N=IWORK(LSCR-1+ISCR)
           IF(N.GT.0) THEN
             NSSTP=NSSTP+1
             NSBSTOT=NSBSTOT+N
           END IF
          END DO
         END DO
        END DO
      END DO
C We finally know the number of substring types and the number of
C substrings. Transfer non-zero entries to the final table.

C Size of table, and offsets:
      KSSTP  =15
      KSSTANN=KSSTP   + 5*NSSTP
      KSSTCRE=KSSTANN + MORSBITS*NSSTP
      KSBSMRS=KSSTCRE + MORSBITS*NSSTP
      KMRSSBS=KSBSMRS + 2*NSBSTOT
      KSBSANN=KMRSSBS + 2*(2**MORSBITS)*NASPRT
      KSBSCRE=KSBSANN + NSBSTOT*MORSBITS
      NTAB   =KSBSCRE + NSBSTOT*MORSBITS -1
      CALL GETMEM('SbStrTab','Allo','Inte',LTAB,NTAB)
      CALL ICOPY(NTAB,0,0,IWORK(LTAB),1)
C The header data
      IWORK(LTAB+ 0)=NTAB
      IWORK(LTAB+ 1)=ITYPE
      IWORK(LTAB+ 2)=LORBTAB
      IWORK(LTAB+ 3)=NSYM
      IWORK(LTAB+ 4)=NASPRT
      IWORK(LTAB+ 5)=MORSBITS
      IWORK(LTAB+ 6)=NSSTP
      IWORK(LTAB+ 7)=NSBSTOT
      IWORK(LTAB+ 8)=KSSTANN
      IWORK(LTAB+ 9)=KSSTCRE
      IWORK(LTAB+10)=KSBSMRS
      IWORK(LTAB+11)=KMRSSBS
      IWORK(LTAB+12)=KSBSANN
      IWORK(LTAB+13)=KSBSCRE
C Fill in the Substring Type table
C Change the counter array into an array of offsets. Also allocate
C a translation table (POP,MS2,ISYM) to Substring Type.
      CALL GETMEM('SSTScr2','Allo','Inte',LSCR2,NSCR)
      ISSTP=0
      ISBS=0
      DO ISPART=1,NASPRT
        NO=IWORK(LNOSUB-1+ISPART)
        DO IPOP=0,NO
         DO ISYM=1,NSYM
          DO IMS2=-NO,NO
           ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                      (ISPART-1)))
           NSBS=IWORK(LSCR-1+ISCR)
           IWORK(LSCR-1+ISCR)=-1
           IWORK(LSCR2-1+ISCR)=-1
           IF(NSBS.GT.0) THEN
             ISSTP=ISSTP+1
             IWORK(LSCR-1+ISCR)=ISBS
             IWORK(LSCR2-1+ISCR)=ISSTP
             ISBS=ISBS+NSBS
             IWORK(LTAB-1+KSSTP+ 0+5*(ISSTP-1))=NSBS
             IWORK(LTAB-1+KSSTP+ 1+5*(ISSTP-1))=IPOP
             IWORK(LTAB-1+KSSTP+ 2+5*(ISSTP-1))=ISYM
             IWORK(LTAB-1+KSSTP+ 3+5*(ISSTP-1))=IMS2
             IWORK(LTAB-1+KSSTP+ 4+5*(ISSTP-1))=ISPART
           END IF
          END DO
         END DO
        END DO
      END DO

C Now produce all possible substrings:
CTEST      write(*,*)' Test in NEWSSTAB, producing substrings.'
CTEST      write(*,'(1x,a,8i5)')'IWORK(LOSYM-1+ISORB):',
CTEST     &                           (IWORK(LOSYM-1+ISORB),ISORB=1,NASPO)
CTEST      write(*,'(1x,a,8i5)')'IWORK(LOSPN-1+ISORB):',
CTEST     &                           (IWORK(LOSPN-1+ISORB),ISORB=1,NASPO)
      ISORB=1
      DO ISPART=1,NASPRT
        NO=IWORK(LNOSUB-1+ISPART)
        DO MRS=0,2**NO-1
          IPOP=MorsPop(MRS)
          ISYM=MorsSymm(MRS,IWORK(LOSYM-1+ISORB))
          IMS2=MorsSpin(MRS,IWORK(LOSPN-1+ISORB))
C Which substring is this?
          ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                     (ISPART-1)))
          ISBS=1+IWORK(LSCR-1+ISCR)
CTEST      write(*,'(1x,a,8i5)')'MRS,IPOP,IMS2,ISYM,ISBS:',
CTEST     &                      MRS,IPOP,IMS2,ISYM,ISBS
          IWORK(LSCR-1+ISCR)=ISBS
          ISSTP=IWORK(LSCR2-1+ISCR)
C Fill in the Substring/Morsel translation arrays:
          LPOS=LTAB-1+KSBSMRS+2*(ISBS-1)
          IWORK(LPOS)=MRS
          IWORK(LPOS+1)=ISSTP
          LPOS=LTAB-1+KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
          IWORK(LPOS)=ISBS
          IWORK(LPOS+1)=ISSTP
        END DO
        ISORB=ISORB+NO
      END DO
      CALL GETMEM('SSTScr','Free','Inte',LSCR,NSCR)
      CALL GETMEM('SSTScr2','Free','Inte',LSCR2,NSCR)
      CALL GETMEM('OrbSpn','Free','Inte',LOSPN,NASPO)
      CALL GETMEM('OrbSym','Free','Inte',LOSYM,NASPO)
      DO ISPART=1,NASPRT
       NO=IWORK(LNOSUB-1+ISPART)
       DO MRS=0,2**NO-1
        LPOS=LTAB-1+KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
        ISBS=IWORK(LPOS)
        ISSTP=IWORK(LPOS+1)
       END DO
      END DO

C Create the Substring Annihilator and Creator arrays, and
C also Substring Type Ann/Cre arrays.
CTEST      write(*,*)' Making annih and creat arrays:'
      DO ISPART=1,NASPRT
       NO=IWORK(LNOSUB-1+ISPART)
       DO MRS=0,2**NO-1
        LPOS=LTAB-1+KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
        ISBS=IWORK(LPOS)
        ISSTP=IWORK(LPOS+1)
        DO I=1,NO
         NEWMRS=MORSANN(MRS,I)
         IF(NEWMRS.NE.999999) THEN
          ISGN=1
          IF(NEWMRS.LT.0) ISGN=-1
          NEWMRS=ISGN*NEWMRS
C Position in Morsel-to-Substring Table
          LPOS=LTAB-1+KMRSSBS+2*(NEWMRS+(2**MORSBITS)*(ISPART-1))
C New substring times sign factor
          NEWSBS=ISGN*IWORK(LPOS)
          NEWSST=IWORK(LPOS+1)
C Put new substring in the table
          LPOS=LTAB-1+KSBSANN-1+I+MORSBITS*(ISBS-1)
          IWORK(LPOS)=NEWSBS
C Put new substring type in the Substring Type Annihil Table
          LPOS=LTAB-1+KSSTANN-1+I+MORSBITS*(ISSTP-1)
          IWORK(LPOS)=NEWSST
         END IF
        END DO
C Create: Very similar to the Annihilate code above.
        DO I=1,NO
         NEWMRS=MORSCRE(MRS,I)
         IF(NEWMRS.NE.999999) THEN
          ISGN=1
          IF(NEWMRS.LT.0) ISGN=-1
          NEWMRS=ISGN*NEWMRS
C Position in Morsel-to-Substring Table
          LPOS=LTAB-1+KMRSSBS+2*(NEWMRS+(2**MORSBITS)*(ISPART-1))
C New substring times sign factor
          NEWSBS=ISGN*IWORK(LPOS)
          NEWSST=IWORK(LPOS+1)
C Put new substring in the table
          LPOS=LTAB-1+KSBSCRE-1+I+MORSBITS*(ISBS-1)
          IWORK(LPOS)=NEWSBS
C Put new substring type in the Substring Type Creat Table
          LPOS=LTAB-1+KSSTCRE-1+I+MORSBITS*(ISSTP-1)
          IWORK(LPOS)=NEWSST
         END IF
        END DO
       END DO
      END DO

C Test section, may be removed later..
      IF(.FALSE.) THEN
C For all strings, check the anticommutation rules.
      IERR=0
      DO ISBS=1,NSBSTOT
C Its substring type
       LPOS=LTAB-1+KSBSMRS+2*(ISBS-1)
       ISSTP=IWORK(LPOS+1)
C Its subpartition
       ISPART=IWORK(LTAB-1+KSSTP+ 4+5*(ISSTP-1))
       NO=IWORK(LNOSUB-1+ISPART)
       DO I=1,NO
C Annihilate orbital nr I in the subpartition.
        LPOS=LTAB-1+KSBSANN-1+I+MORSBITS*(ISBS-1)
        ISBS1=IWORK(LPOS)
C Create orbital nr I in the subpartition.
        LPOS=LTAB-1+KSBSCRE-1+I+MORSBITS*(ISBS-1)
        ISBS2=IWORK(LPOS)
        IF(ISBS1.NE.0 .AND. ISBS2.NE.0) IERR=IERR+1
        IF(ISBS1.EQ.0 .AND. ISBS2.EQ.0) IERR=IERR+1
        DO J=1,I
C Same, for orbital J instead.
         LPOS=LTAB-1+KSBSANN-1+J+MORSBITS*(ISBS-1)
         ISBS3=IWORK(LPOS)
         LPOS=LTAB-1+KSBSCRE-1+J+MORSBITS*(ISBS-1)
         ISBS4=IWORK(LPOS)
C SBS5=A(I-)*SBS4:
         ISBS5=0
         IF(ISBS4.NE.0) THEN
           ISGN=1
           IF(ISBS4.LT.0) ISGN=-1
           LPOS=LTAB-1+KSBSANN-1+I+MORSBITS*(ABS(ISBS4)-1)
           ISBS5=ISGN*IWORK(LPOS)
         END IF
C SBS6=A(J+)*SBS1:
         ISBS6=0
         IF(ISBS1.NE.0) THEN
           ISGN=1
           IF(ISBS1.LT.0) ISGN=-1
           LPOS=LTAB-1+KSBSCRE-1+J+MORSBITS*(ABS(ISBS1)-1)
           ISBS6=ISGN*IWORK(LPOS)
         END IF
C SBS7=A(I+)*SBS3:
         ISBS7=0
         IF(ISBS3.NE.0) THEN
           ISGN=1
           IF(ISBS3.LT.0) ISGN=-1
           LPOS=LTAB-1+KSBSCRE-1+I+MORSBITS*(ABS(ISBS3)-1)
           ISBS7=ISGN*IWORK(LPOS)
         END IF
C SBS8=A(J-)*SBS2:
         ISBS8=0
         IF(ISBS2.NE.0) THEN
           ISGN=1
           IF(ISBS2.LT.0) ISGN=-1
           LPOS=LTAB-1+KSBSANN-1+J+MORSBITS*(ABS(ISBS2)-1)
           ISBS8=ISGN*IWORK(LPOS)
         END IF
C Now check:
         IF(I.NE.J) THEN
           IF(ISBS5+ISBS6.NE.0) IERR=IERR+1
           IF(ISBS7+ISBS8.NE.0) IERR=IERR+1
         ELSE
           IF(ISBS5+ISBS6.NE.ISBS) IERR=IERR+1
           IF(ISBS7+ISBS8.NE.ISBS) IERR=IERR+1
         END IF
        END DO
       END DO
      END DO
      WRITE(6,*)' NEWSSTAB: Nr of test errors IERR=',IERR
      END IF


      NEWSSTAB=LTAB
      CALL GETMEM('NOSub','Free','Inte',LNOSUB,NASPRT)
      RETURN
      END
