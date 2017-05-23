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
      SUBROUTINE SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,CISYG,CISD)
      IMPLICIT NONE
      REAL*8 CISYG(*),CISD(*)
      INTEGER ICNFTAB(*),ISPNTAB(*),ISSTAB(*),IFSBTAB(*)
      INTEGER NASPRT
      INTEGER MINOP,MAXOP,NACTEL,NOPEN,NCLSD
      INTEGER NCNF,NSPD,NCPL,ISPART,ISST
C     INTEGER KFSB,IBLK,ISPD,I,IPOS,IORB,ISYM
      INTEGER KFSB,IBLK,ISPD,I,IPOS,IORB
      INTEGER ISPN,NO,IOSTA,IOEND,IMORS,ISBSTR
#include "WrkSpc.fh"
      INTEGER ICNF,IEL,IEL1,IEL2,IFORM,IFSB,IOCC
      INTEGER IPART,IREST
C     INTEGER ISTARR,LDIM,LOCARR,LORBARR,LSBSET,LSSARR,LSTARR
      INTEGER        LDIM,LOCARR,LORBARR,LSBSET,LSSARR,LSTARR
      INTEGER ISORB,ISPEND,ISPSTA,ISUM,ISYGEND,ISYGSTA,NWRD
      INTEGER IWORD,IWRD,JSST,KCNF,KCNFINF,KGSLIM,KGSORB,KHSHMAP
      INTEGER KMRSSBS,KSPNINF,KSSTARR,KSSTTB,LBLK,KSPN,LSPTRA
      INTEGER LSYM,MORSBITS,MXBLK,NAPART,NBLK,NHEAD
      INTEGER NHSHMAP,NOCC,NOP,NORB,NSP,NSSTP,NSYM
C     INTEGER IERR,ICPL,KSBSMRS,JMORS,NFSB
      INTEGER IERR,     KSBSMRS,      NFSB
      INTEGER OCC2MRS
      EXTERNAL OCC2MRS

C Unbutton the configuration table:
      NACTEL=ICNFTAB(3)
      NORB  =ICNFTAB(4)
      MINOP =ICNFTAB(5)
      MAXOP =ICNFTAB(6)
      NSYM  =ICNFTAB(7)
      LSYM  =ICNFTAB(8)
      NAPART=ICNFTAB(9)
      IFORM =ICNFTAB(10)
      NHEAD=10
      KGSORB =NHEAD+1
      KGSLIM =KGSORB+(NSYM+1)*(NAPART+1)
      KCNFINF =KGSLIM+2*NAPART
CTEST      write(*,*)' Table at KGSORB:'
CTEST      do ipart=0,napart
CTEST       write(*,'(1x,10i5)')(icnftab(kgsorb+isym+(nsym+1)*ipart)
CTEST     &                                             ,isym=0,nsym)
CTEST      end do
CTEST      write(*,*)' Table at KGSLIM:'
CTEST      write(*,'(1x,10i5)')(icnftab(kgslim+2*(ipart-1)),ipart=1,napart)
CTEST      write(*,'(1x,10i5)')(icnftab(kgslim+1+2*(ipart-1)),ipart=1,napart)
C Unbutton the spin coupling table:
      KSPNINF=9
C Unbutton the Substring Table:
      NASPRT=ISSTAB(5)
      MORSBITS=ISSTAB(6)
      NSSTP =ISSTAB(7)
      KSSTTB=15
      KSBSMRS=ISSTAB(11)
      KMRSSBS=ISSTAB(12)
C Unbutton the Fock Sector Block table:
      NHEAD=7
      KSSTARR=NHEAD+1
      NFSB   =IFSBTAB(2)
      NHSHMAP=IFSBTAB(6)
      KHSHMAP=IFSBTAB(7)

C MXBLK=Largest individual SYG block of determinants:
      MXBLK=0
      DO NOPEN=MINOP,MAXOP
        NCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
        NCPL=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2)
        MXBLK=MAX(NCNF*NCPL,MXBLK)
      END DO
C A variety of small temporary arrays used:
      CALL GETMEM('TmpArr','Allo','Real',LBLK,MXBLK)
      CALL GETMEM('OrbArr','Allo','Inte',LORBARR,NACTEL)
      CALL GETMEM('OccArr','Allo','Inte',LOCARR,2*NORB)
      CALL GETMEM('STArr','Allo','Inte',LSTARR,NASPRT)
      CALL GETMEM('Dims','Allo','Inte',LDIM,NASPRT)
      CALL GETMEM('SSArr','Allo','Inte',LSSARR,NASPRT)
      CALL GETMEM('NSBSET','Allo','Inte',LSBSET,NSSTP)
C We will need later the accumulated number of substrings of
C earlier substring types:
      ISUM=0
      DO ISST=1,NSSTP
        IWORK(LSBSET-1+ISST)=ISUM
        ISUM=ISUM+ISSTAB(KSSTTB+5*(ISST-1))
      END DO
C Loop over nr of open shells.
      ISYGEND=0
      DO NOPEN=MINOP,MAXOP
        NCLSD=(NACTEL-NOPEN)/2
        IF(NCLSD.LT.0) GOTO 100
        IF(2*NCLSD+NOPEN.NE.NACTEL) GOTO 100
        NOCC=NCLSD+NOPEN
        IF(NOCC.GT.NORB) GOTO 100
        NCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
        IF(NCNF.EQ.0) GOTO 100
        KCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+1)
        NWRD=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+2)
        NCPL=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+1)
        NSPD=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2)
        KSPN=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+4)
        IF(NSPD.EQ.0) GOTO 100
        IF(NCPL.EQ.0) GOTO 100
C ISYGSTA=1st element of each block
        NBLK=NCPL*NCNF
        ISYGSTA=ISYGEND+1
        ISYGEND=ISYGEND+NBLK
C Location of spin coupling coefficients:
        LSPTRA=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+5)
C Matrix multiplication into temporary array:
        CALL  DGEMM_('N','N',NSPD,NCNF,NCPL,1.0D0,
     &               WORK(LSPTRA),NSPD,
     &               CISYG(ISYGSTA),NCPL,0.0D0,
     &               WORK(LBLK),NSPD)

C There is no phase factor in the reorder of orbitals from SYG
C to SD. But there is a fairly lengthy procedure for finding the
C correct position of each CI coefficient.
        IBLK=0
C Loop over configurations
        IWORD = 0 ! dummy initialize
        DO ICNF=1,NCNF
          IF(IFORM.EQ.1) THEN
            DO IEL=1,NOCC
              IWORK(LORBARR-1+IEL)=ICNFTAB(KCNF-1+IEL+NWRD*(ICNF-1))
            END DO
          ELSE IF(IFORM.EQ.2) THEN
            IEL2=0
            IEL1=NCLSD
            DO IORB=1,NORB
              IOCC=ICNFTAB(KCNF-1+IORB+NWRD*(ICNF-1))
              IF(IOCC.EQ.1) THEN
                IEL1=IEL1+1
                IWORK(LORBARR-1+IEL1)=IORB
              ELSE
                IEL2=IEL2+1
                IWORK(LORBARR-1+IEL2)=IORB
              END IF
            END DO
          ELSE IF(IFORM.EQ.3) THEN
            DO IEL=1,NOCC
              IWRD=(3+IEL)/4
              IREST=(3+IEL)-4*IWRD
              IF(IREST.EQ.0) THEN
                IWORD=ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
              END IF
              IORB=MOD(IWORD,256)
              IWORD=IWORD/256
              IWORK(LORBARR-1+IEL)=IORB
            END DO
          ELSE IF(IFORM.EQ.4) THEN
            IEL2=0
            IEL1=NCLSD
            DO IORB=1,NORB
              IWRD=(IORB+14)/15
              IREST=IORB+14-15*IWRD
              IF(IREST.EQ.0) THEN
                IWORD=ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
              END IF
              IOCC=MOD(IWORD,4)
              IWORD=IWORD/4
              IF(IOCC.EQ.1) THEN
                IEL1=IEL1+1
                IWORK(LORBARR-1+IEL1)=IORB
              ELSE
                IEL2=IEL2+1
                IWORK(LORBARR-1+IEL2)=IORB
              END IF
            END DO
          END IF

CTEST      write(*,'(1x,a,10i5)')'Configuration:',
CTEST     &                          (iwork(lorbarr-1+iel),iel=1,nocc)
CTEST      write(*,*)' Loop over spin determinants.'
C Loop over spin determinants
          DO ISPD=1,NSPD
CTEST      write(*,'(1x,a,20i3)')'Spin determinant:',
CTEST     &                    (ISPNTAB(KSPN-1+I+NOPEN*(ISPD-1)),I=1,nopen)
            IBLK=IBLK+1
C Construct occupation number array:
            CALL ICOPY(2*NORB,0,0,IWORK(LOCARR),1)
            DO IEL=1,NCLSD
              IORB=IWORK(LORBARR-1+IEL)
              IWORK(LOCARR-1+2*IORB-1)=1
              IWORK(LOCARR-1+2*IORB  )=1
            END DO
            DO I=1,NOPEN
C Spin of each electron is coded as 1 for alpha, 0 for beta.
              ISPN=ISPNTAB(KSPN-1+I+NOPEN*(ISPD-1))
              IEL=NCLSD+I
              ISORB=2*IWORK(LORBARR-1+IEL)-ISPN
              IWORK(LOCARR-1+ISORB)=1
            END DO
C Identify substrings:
C Loop over active partitions. Subdivide as needed into subpartitions.
CTEST      write(*,*)' Identify substrings.'
CTEST      write(*,'(1x,a,10i5)')'Occupation array:',
CTEST     &                          (iwork(locarr-1+isorb),isorb=1,2*norb)
            IOEND=0
            ISPEND=0
CTEST      write(*,'(1x,a,10i5)')'NAPART:',NAPART
            DO IPART=1,NAPART
CTEST      write(*,'(1x,a,10i5)')' In loop, IPART:',IPART
             NOP=2*ICNFTAB(KGSORB+(NSYM+1)*IPART)
CTEST      write(*,'(1x,a,10i5)')'            NOP:',NOP
CTEST             IF(NOP.EQ.0) write(*,*)' (Skip it.)'
             IF(NOP.EQ.0) GOTO 200
             NSP=(NOP+MORSBITS-1)/MORSBITS
CTEST      write(*,'(1x,a,10i5)')'            NSP:',NSP
             ISPSTA=ISPEND+1
             ISPEND=ISPEND+NSP
             DO ISPART=ISPSTA,ISPEND
CTEST      write(*,'(1x,a,10i5)')'Loop lims ISPSTA,ISPEND:',ISPSTA,ISPEND
              NO=MIN(NOP,MORSBITS)
              NOP=NOP-NO
              IOSTA=IOEND+1
              IOEND=IOEND+NO
CTEST      write(*,'(1x,a,10i5)')'IOSTA,IOEND:',IOSTA,IOEND
CTEST      write(*,'(1x,a,10i5)')'Occ array:',
CTEST     &                     (IWORK(LOCARR-1+ISORB),ISORB=IOSTA,IOEND)
              IMORS=OCC2MRS(NO,IWORK(LOCARR-1+IOSTA))
CTEST      write(*,'(1x,a,10i5)')'IMORS=',IMORS
C Position in Morsel-to-Substring table:
              IPOS=KMRSSBS+2*(IMORS+(2**MORSBITS)*(ISPART-1))
C Substring ID number
              ISBSTR=ISSTAB(IPOS)
C Test:
CTEST              JMORS=ISSTAB(KSBSMRS+2*(ISBSTR-1))
CTEST              IF(IMORS.NE.JMORS)THEN
CTEST                WRITE(*,*)' Mistranslated morsel!!'
CTEST      write(*,'(1x,a,4i12)')'IMORS->ISBSTR:',IMORS,ISBSTR
CTEST      write(*,'(1x,a,4i12)')'but ISBSTR->IMORS:',ISBSTR,JMORS
CTEST      write(*,'(1x,a,4i12)')'KMRSSBS:',KMRSSBS
CTEST      write(*,'(1x,a,4i12)')'KSBSMRS:',KSBSMRS
CTEST                CALL ABEND()
CTEST              END IF
CTEST      write(*,'(1x,a,10i5)')'ISBSTR:',ISBSTR
C Substring type ISST, nr of such substrings is NDIM
              ISST=ISSTAB(IPOS+1)
CTEST      write(*,'(1x,a,10i5)')'ISST  :',ISST
              IWORK(LSTARR-1+ISPART)=ISST
              IWORK(LDIM-1+ISPART)=ISSTAB(KSSTTB+5*(ISST-1))
              IWORK(LSSARR-1+ISPART)=ISBSTR-IWORK(LSBSET-1+ISST)
             END DO
 200         CONTINUE
            END DO
CTEST      write(*,*)' Finally, substring types and substrings:'
CTEST      write(*,'(1x,a,10i5)')'Substr types:',
CTEST     &                   (IWORK(LSTARR-1+ISPART),ISPART=1,NASPRT)
CTEST      write(*,'(1x,a,10i5)')'Substrings  :',
CTEST     &                    (IWORK(LSSARR-1+ISPART),ISPART=1,NASPRT)
CTEST      write(*,'(1x,a,10i5)')'Dimensions  :',
CTEST     &                    (IWORK(LDIM  -1+ISPART),ISPART=1,NASPRT)
C Position within FS block:
            IPOS=(IWORK(LSSARR-1+NASPRT)-1)
            DO ISPART=NASPRT-1,1,-1
              IPOS=IWORK(LDIM-1+ISPART)*IPOS+(IWORK(LSSARR-1+ISPART)-1)
            END DO
            IPOS=IPOS+1
C Identify Fock Sector Block:
CTEST      write(*,*)' Arguments in HSHGET call:'
CTEST      write(*,'(1x,a,10i5)')'Key:',
CTEST     &                   (IWORK(LSTARR-1+ISPART),ISPART=1,NASPRT)
CTEST      write(*,'(1x,a,10i5)')'Size of key:',NASPRT
CTEST      write(*,'(1x,a,10i5)')'Size of items stored:',NASPRT+2
CTEST      write(*,'(1x,a,10i5)')'Items stored at KSSTARR=',KSSTARR
CTEST      write(*,'(1x,a,10i5)')'      Map size  NHSHMAP=',NHSHMAP
CTEST      write(*,'(1x,a,10i5)')'  Map stored at KHSHMAP=',KHSHMAP
            CALL HSHGET(IWORK(LSTARR),NASPRT,NASPRT+2,IFSBTAB(KSSTARR),
     &                            NHSHMAP,IFSBTAB(KHSHMAP),IFSB)
CTEST      write(*,'(1x,a,10i5)')' Map returns index IFSB=',IFSB
CTEST      write(*,'(1x,a,10i5)')' Item stored there is  =',
CTEST     &                 (IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1)),
CTEST     &                                              ISPART=1,NASPRT+2)
C Position of this FS block in SD wave function:
            KFSB=IFSBTAB(KSSTARR+(NASPRT+2)*IFSB-1)
C Temporary check, may be removed later. See that we have picked up
C the correct FS block.
            IERR=0
            DO ISPART=1,NASPRT
              JSST=IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1))
              ISST=IWORK(LSTARR-1+ISPART)
              IF(ISST.NE.JSST) IERR=1
            END DO
            IF(IERR.NE.0) THEN
              WRITE(6,*) ' SYGTOSD Error:'//
     &                     ' Hash map returned the wrong FS block!'
              WRITE(6,'(1x,a,8I8)')'NOPEN,ICNF,ISPN:',NOPEN,ICNF,ISPN
              WRITE(6,'(1x,a,20I3)')'Configuration:',
     &                               (IWORK(LORBARR-1+IEL),IEL=1,NACTEL)
              WRITE(6,'(1x,a,20I3)')'Determinant:',
     &                            (IWORK(LOCARR-1+ISORB),ISORB=1,2*NORB)
              WRITE(6,'(1x,a,10I5)')'Substring type combination:',
     &                         (IWORK(LSTARR-1+ISPART),ISPART=1,NASPRT)
              WRITE(6,'(1x,a,10I5)')'Substring combination:',
     &                         (IWORK(LSSARR-1+ISPART),ISPART=1,NASPRT)
              WRITE(6,'(1x,a,8I8)')'Hash table says IFSB=',IFSB
              IF(IFSB.GT.0 .AND. IFSB.LE.NFSB) THEN
              WRITE(6,'(1x,a,8I8)')'but that FS block would contain',
     & (IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1)),ISPART=1,NASPRT)
              ELSE
              WRITE(6,*)'but there is no such FS block!'
              WRITE(6,*)' The FS block table follows:'
              CALL PRFSBTAB(IFSBTAB)
              END IF
              CALL ABEND()
            END IF
C Finally:
            CISD(KFSB-1+IPOS)=WORK(LBLK-1+IBLK)
C End of spin-determinant loop
          END DO
C End of loop over configurations
        END DO
C End of loop over nr of open shells
 100    CONTINUE
      END DO
      CALL GETMEM('TmpArr','Free','Real',LBLK,MXBLK)
      CALL GETMEM('OrbArr','Free','Inte',LORBARR,NACTEL)
      CALL GETMEM('OccArr','Free','Inte',LOCARR,NORB)
      CALL GETMEM('STArr','Free','Inte',LSTARR,NASPRT)
      CALL GETMEM('Dims','Free','Inte',LDIM,NASPRT)
      CALL GETMEM('SSArr','Free','Inte',LSSARR,NASPRT)
      CALL GETMEM('NSBSET','Free','Inte',LSBSET,NSSTP)
      RETURN
      END
