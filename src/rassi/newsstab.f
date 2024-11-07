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
      Subroutine NEWSSTAB(ORBTAB)
      use stdalloc, only: mma_allocate, mma_deallocate
      use rassi_global_arrays, only: SSTAB
      use cntrl, only: MORSBITS
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      IMPLICIT NONE
      INTEGER ORBTAB(*)

      INTEGER NTAB,ITYPE,NSSTP,NSBSTOT,KSBSMRS,KMRSSBS
      INTEGER KSBSANN,KSBSCRE
      INTEGER I,IKMS2,IKSCR,IKSYM
      INTEGER IMS2,IPOP,ISBS,ISCR
      INTEGER ISGN,ISORB,ISPART,ISSTP
      INTEGER ISYM,J,KMS2,KOINFO,KSSTANN
      INTEGER KSSTCRE,KSSTP,KSYM
      INTEGER LPOS
      INTEGER MRS,MS2
      INTEGER MXMS2,MXPOP,N
      INTEGER NEWMRS,NEWSBS,NEWSST,NO
      INTEGER NSBS,NSCR,NASPO,NASPRT
      INTEGER IERR,ISBS1,ISBS2,ISBS3,ISBS4,ISBS5,ISBS6,ISBS7,ISBS8
      INTEGER, EXTERNAL :: MORSANN,MORSCRE,MORSPOP,MORSSPIN,MORSSYMM
      Integer, allocatable:: OSPN(:), OSYM(:), NOSUB(:), SCR(:), SCR2(:)

C Table type ID:
      ITYPE=19
C Pick up some data from the orbital table:
      NASPO =ORBTAB(4)
      NSYM  =ORBTAB(5)
      NASPRT=ORBTAB(9)
      KOINFO=19
C Make temporary arrays for spin label and symmetry of each orbital
      CALL mma_allocate(OSPN,NASPO,Label='OSPN')
      CALL mma_allocate(OSYM,NASPO,Label='OSYM')
C Make a temporary array, which gives the number of spin orbitals in
C each subpartition:
      CALL mma_allocate(NOSUB,NASPRT,Label='NOSUB')
      NOSUB(:)=0
      DO ISORB=1,NASPO
        ISYM  = OrbTab(KOINFO+ 1+(ISORB-1)*8)
        OSYM(ISORB)=ISYM
        MS2   = OrbTab(KOINFO+ 3+(ISORB-1)*8)
        OSPN(ISORB)=MS2
        ISPART= OrbTab(KOINFO+ 6+(ISORB-1)*8)
        N=1+NOSUB(ISPART)
        NOSUB(ISPART)=N
      END DO
C We need a temporary array, NSBSSCR(NSYM,0:NPOP,-MXMS2:MXMS2,NASPRT),
C to keep the number of substrings of different kind:
      MXPOP=MORSBITS
      MXMS2=MORSBITS
      NSCR=NSYM*(1+MXPOP)*(2*MXMS2+1)*NASPRT
      CALL mma_allocate(SCR,NSCR,Label='SCR')
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
           SCR(ISCR)=0
          END DO
         END DO
        END DO
C NSBSSCR(ISYM=1,IPOP=0,IMS2=0,ISPART)=1:
        ISCR=1+NSYM*(0+(1+MXPOP)*(MXMS2+0+(2*MXMS2+1)*(ISPART-1)))
        SCR(ISCR)=1
      END DO
C Compute number of substrings:
      ISORB=0
      DO ISPART=1,NASPRT
        NO=NOSUB(ISPART)
        DO I=1,NO
          ISORB=ISORB+1
          KSYM=OSYM(ISORB)
          KMS2=OSPN(ISORB)
          DO IPOP=I,1,-1
           DO IMS2=-IPOP,IPOP
            IKMS2=IMS2-KMS2
            IF(ABS(IKMS2).LE.IPOP-1) THEN
            DO ISYM=1,NSYM
             IKSYM=MUL(ISYM,KSYM)
             ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                      (ISPART-1)))
             N=SCR(ISCR)
             IKSCR=IKSYM+NSYM*(IPOP-1+(1+MXPOP)*(MXMS2+IKMS2+
     &                          (2*MXMS2+1)*(ISPART-1)))
             N=N+SCR(IKSCR)
             SCR(ISCR)=N
            END DO
            END IF
           END DO
          END DO
        END DO
      END DO
      NSSTP=0
      NSBSTOT=0
      DO ISPART=1,NASPRT
        NO=NOSUB(ISPART)
        DO IPOP=0,NO
         DO ISYM=1,NSYM
          DO IMS2=-NO,NO
           ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                     (ISPART-1)))
           N=SCR(ISCR)
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
      CALL mma_allocate(SSTAB,NTAB,Label='SSTAB')
      SSTAB(:)=0
C The header data
      SSTAB( 1)=NTAB
      SSTAB( 2)=ITYPE
      SSTAB( 3)= -1 ! Not used
      SSTAB( 4)=NSYM
      SSTAB( 5)=NASPRT
      SSTAB( 6)=MORSBITS
      SSTAB( 7)=NSSTP
      SSTAB( 8)=NSBSTOT
      SSTAB( 9)=KSSTANN
      SSTAB(10)=KSSTCRE
      SSTAB(11)=KSBSMRS
      SSTAB(12)=KMRSSBS
      SSTAB(13)=KSBSANN
      SSTAB(14)=KSBSCRE
C Fill in the Substring Type table
C Change the counter array into an array of offsets. Also allocate
C a translation table (POP,MS2,ISYM) to Substring Type.
      CALL mma_allocate(SCR2,NSCR,Label='SCR2')
      ISSTP=0
      ISBS=0
      DO ISPART=1,NASPRT
        NO=NOSUB(ISPART)
        DO IPOP=0,NO
         DO ISYM=1,NSYM
          DO IMS2=-NO,NO
           ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                      (ISPART-1)))
           NSBS=SCR(ISCR)
           SCR(ISCR)=-1
           SCR2(ISCR)=-1
           IF(NSBS.GT.0) THEN
             ISSTP=ISSTP+1
             SCR(ISCR)=ISBS
             SCR2(ISCR)=ISSTP
             ISBS=ISBS+NSBS
             SSTAB(KSSTP+ 0+5*(ISSTP-1))=NSBS
             SSTAB(KSSTP+ 1+5*(ISSTP-1))=IPOP
             SSTAB(KSSTP+ 2+5*(ISSTP-1))=ISYM
             SSTAB(KSSTP+ 3+5*(ISSTP-1))=IMS2
             SSTAB(KSSTP+ 4+5*(ISSTP-1))=ISPART
           END IF
          END DO
         END DO
        END DO
      END DO

C Now produce all possible substrings:
CTEST      write(*,*)' Test in NEWSSTAB, producing substrings.'
CTEST      write(*,'(1x,a,8i5)')'OSYM(ISORB):',
CTEST     &                           (OSYM(ISORB),ISORB=1,NASPO)
CTEST      write(*,'(1x,a,8i5)')'OSPN(ISORB):',
CTEST     &                           (OSPN(ISORB),ISORB=1,NASPO)
      ISORB=1
      DO ISPART=1,NASPRT
        NO=NOSUB(ISPART)
        DO MRS=0,2**NO-1
          IPOP=MorsPop(MRS)
          ISYM=MorsSymm(MRS,OSYM(ISORB))
          IMS2=MorsSpin(MRS,OSPN(ISORB))
C Which substring is this?
          ISCR=ISYM+NSYM*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*
     &                     (ISPART-1)))
          ISBS=1+SCR(ISCR)
CTEST      write(*,'(1x,a,8i5)')'MRS,IPOP,IMS2,ISYM,ISBS:',
CTEST     &                      MRS,IPOP,IMS2,ISYM,ISBS
          SCR(ISCR)=ISBS
          ISSTP=SCR2(ISCR)
C Fill in the Substring/Morsel translation arrays:
          LPOS=KSBSMRS+2*(ISBS-1)
          SSTAB(LPOS)=MRS
          SSTAB(LPOS+1)=ISSTP
          LPOS=KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
          SSTAB(LPOS)=ISBS
          SSTAB(LPOS+1)=ISSTP
        END DO
        ISORB=ISORB+NO
      END DO
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(SCR2)
      CALL mma_deallocate(OSPN)
      CALL mma_deallocate(OSYM)
      DO ISPART=1,NASPRT
       NO=NOSUB(ISPART)
       DO MRS=0,2**NO-1
        LPOS=KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
        ISBS=SSTAB(LPOS)
        ISSTP=SSTAB(LPOS+1)
       END DO
      END DO

C Create the Substring Annihilator and Creator arrays, and
C also Substring Type Ann/Cre arrays.
CTEST      write(*,*)' Making annih and creat arrays:'
      DO ISPART=1,NASPRT
       NO=NOSUB(ISPART)
       DO MRS=0,2**NO-1
        LPOS=KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
        ISBS=SSTAB(LPOS)
        ISSTP=SSTAB(LPOS+1)
        DO I=1,NO
         NEWMRS=MORSANN(MRS,I)
         IF(NEWMRS.NE.999999) THEN
          ISGN=1
          IF(NEWMRS.LT.0) ISGN=-1
          NEWMRS=ISGN*NEWMRS
C Position in Morsel-to-Substring Table
          LPOS=KMRSSBS+2*(NEWMRS+(2**MORSBITS)*(ISPART-1))
C New substring times sign factor
          NEWSBS=ISGN*SSTAB(LPOS)
          NEWSST=SSTAB(LPOS+1)
C Put new substring in the table
          LPOS=KSBSANN-1+I+MORSBITS*(ISBS-1)
          SSTAB(LPOS)=NEWSBS
C Put new substring type in the Substring Type Annihil Table
          LPOS=KSSTANN-1+I+MORSBITS*(ISSTP-1)
          SSTAB(LPOS)=NEWSST
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
          LPOS=KMRSSBS+2*(NEWMRS+(2**MORSBITS)*(ISPART-1))
C New substring times sign factor
          NEWSBS=ISGN*SSTAB(LPOS)
          NEWSST=SSTAB(LPOS+1)
C Put new substring in the table
          LPOS=KSBSCRE-1+I+MORSBITS*(ISBS-1)
          SSTAB(LPOS)=NEWSBS
C Put new substring type in the Substring Type Creat Table
          LPOS=KSSTCRE-1+I+MORSBITS*(ISSTP-1)
          SSTAB(LPOS)=NEWSST
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
       LPOS=KSBSMRS+2*(ISBS-1)
       ISSTP=SSTAB(LPOS+1)
C Its subpartition
       ISPART=SSTAB(KSSTP+ 4+5*(ISSTP-1))
       NO=NOSUB(ISPART)
       DO I=1,NO
C Annihilate orbital nr I in the subpartition.
        LPOS=KSBSANN-1+I+MORSBITS*(ISBS-1)
        ISBS1=SSTAB(LPOS)
C Create orbital nr I in the subpartition.
        LPOS=KSBSCRE-1+I+MORSBITS*(ISBS-1)
        ISBS2=SSTAB(LPOS)
        IF(ISBS1.NE.0 .AND. ISBS2.NE.0) IERR=IERR+1
        IF(ISBS1.EQ.0 .AND. ISBS2.EQ.0) IERR=IERR+1
        DO J=1,I
C Same, for orbital J instead.
         LPOS=KSBSANN-1+J+MORSBITS*(ISBS-1)
         ISBS3=SSTAB(LPOS)
         LPOS=KSBSCRE-1+J+MORSBITS*(ISBS-1)
         ISBS4=SSTAB(LPOS)
C SBS5=A(I-)*SBS4:
         ISBS5=0
         IF(ISBS4.NE.0) THEN
           ISGN=1
           IF(ISBS4.LT.0) ISGN=-1
           LPOS=KSBSANN-1+I+MORSBITS*(ABS(ISBS4)-1)
           ISBS5=ISGN*SSTAB(LPOS)
         END IF
C SBS6=A(J+)*SBS1:
         ISBS6=0
         IF(ISBS1.NE.0) THEN
           ISGN=1
           IF(ISBS1.LT.0) ISGN=-1
           LPOS=KSBSCRE-1+J+MORSBITS*(ABS(ISBS1)-1)
           ISBS6=ISGN*SSTAB(LPOS)
         END IF
C SBS7=A(I+)*SBS3:
         ISBS7=0
         IF(ISBS3.NE.0) THEN
           ISGN=1
           IF(ISBS3.LT.0) ISGN=-1
           LPOS=KSBSCRE-1+I+MORSBITS*(ABS(ISBS3)-1)
           ISBS7=ISGN*SSTAB(LPOS)
         END IF
C SBS8=A(J-)*SBS2:
         ISBS8=0
         IF(ISBS2.NE.0) THEN
           ISGN=1
           IF(ISBS2.LT.0) ISGN=-1
           LPOS=KSBSANN-1+J+MORSBITS*(ABS(ISBS2)-1)
           ISBS8=ISGN*SSTAB(LPOS)
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


      CALL mma_deallocate(NOSUB)

      END Subroutine NEWSSTAB
