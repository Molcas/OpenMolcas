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
      SUBROUTINE VERTAB(NACTEL,M2SPIN,LSYM,NPART,NGASORB,NGASLIM,
     &                  ISSTAB,NFSB0,NRDETS0,NFSB,NRDETS,LFSBARR)
      IMPLICIT NONE
      INTEGER ISSTAB(*)
      INTEGER NASPRT,NACTEL,M2SPIN,LSYM
      INTEGER NPART
#include "symmul.fh"
      INTEGER NGASORB(0:NSYM,0:NPART)
      INTEGER NGASLIM(2,NPART)
      INTEGER NFSB,NFSB0,NRDETS,NRDETS0
      INTEGER NVERT,IVER1,IVER2
      INTEGER IPOS,ISPART,NPC1,M2C1,ISC1
      INTEGER ISST,KPOP,KMS2,KSYM
      INTEGER MX1,MN1,MX2,MN2,IPART,NOREM,MINEL,MAXEL
      INTEGER ISPEND,NGO,NSP,ISPSTA,MNL,MXL,MNR,MXR,NO
      INTEGER NPC2,M2C2,ISC2
      INTEGER NSSTP,KSSTP,KPOPMAX,KMS2MAX,NPOP,MS2,NSSTPTR
      INTEGER LSSTPTR,LSSTARR,N,IEND,ISTA,IADDR,ISAVE,NEXT
      INTEGER MS2MIN,MS2MAX,MSFACT,NVIDX,LVIDX,IVIDX,IVERT
      INTEGER LEVUP,LEVDWN,IADDR1,IADDR2,NARC
      INTEGER NARCVRT,NVERTAB,LVERTAB,NDWNTAB,LDWNTAB,NARCLEV,LWEIGHT
      INTEGER NFSBARR,LFSBARR,IVUP,LOWEST,ISW,NSW
      INTEGER KPOPLIM,LLTDA,LSWITCH
      INTEGER IARC,IA,ISUM,KEEP
C     INTEGER NSBS,LDUM,NDUM
      INTEGER NSBS
      INTEGER NSDBLK
#include "Morsel.fh"
#include "WrkSpc.fh"
      INTEGER LIMARR(2,0:50)
      INTEGER ITRY(50)

*----------------------------------------------------------------
* Unbutton the substring table:
      NSYM  =ISSTAB(4)
      NASPRT=ISSTAB(5)
      NSSTP =ISSTAB(7)
      KSSTP =15

CTEST      write(*,*)' Test print in VERTAB.'
CTEST      write(*,*)' NSYM  :',NSYM
CTEST      write(*,*)' NASPRT:',NASPRT
CTEST      write(*,*)' NSSTP :',NSSTP
CTEST      write(*,*)' Substrings:'
CTEST      do isst=1,nsstp
CTEST       NSBS   = ISSTAB(KSSTP+0+5*(ISST-1))
CTEST       NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
CTEST       KSYM   = ISSTAB(KSSTP+2+5*(ISST-1))
CTEST       MS2    = ISSTAB(KSSTP+3+5*(ISST-1))
CTEST       ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
CTEST       write(*,'(1x,6i8)')isst,nsbs,npop,ksym,ms2,ispart
CTEST      end do
* Now the properties of the substring types can be picked up
* as
*       NSBS   = ISSTAB(KSSTP+0+5*(ISST-1))
*       NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
*       ISYM   = ISSTAB(KSSTP+2+5*(ISST-1))
*       MS2    = ISSTAB(KSSTP+3+5*(ISST-1))
*       ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
*---------------------------------------------------------------
* Determine KPOPMAX = Max population of any subpartition
* Determine KMS2MAX = Max spin proj in any subpartition
      KPOPMAX=0
      KMS2MAX=0
      DO ISST=1,NSSTP
       NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
       MS2    = ISSTAB(KSSTP+3+5*(ISST-1))
       KPOPMAX=MAX(NPOP,KPOPMAX)
       KMS2MAX=MAX(MS2,KMS2MAX)
      END DO
CTEST      write(*,'(1x,a,8i8)')'KPOPMAX:',KPOPMAX
CTEST      write(*,'(1x,a,8i8)')'KMS2MAX:',KMS2MAX
*----------------------------------------------------------------
* Using the GAS restrictions, set up an array with min and max nr of
* electrons (cumulative) up to and including subpartition ISPART:
      MX1=0
      MN1=0
      MX2=NACTEL
      MN2=NACTEL
      DO IPART=1,NPART
       NOREM=NGASORB(0,IPART)
       MINEL=MAX(0,NGASLIM(1,IPART))
       MAXEL=MIN(NOREM,NGASLIM(2,IPART))
       MX2=MX2-MINEL
       MN2=MN2-MAXEL
      END DO
      ISPEND=0

      LIMARR(1,0)=0
      LIMARR(2,0)=0
      DO IPART=1,NPART
       NGO  =NGASORB(0,IPART)
       MINEL=MAX(0,NGASLIM(1,IPART))
       MAXEL=MIN(NGO,NGASLIM(2,IPART))
       NSP=(NGO+MORSBITS-1)/MORSBITS
       ISPSTA=ISPEND+1
       ISPEND=ISPEND+NSP
       MNL=MAX(MN1,MN2)
       MXL=MIN(MX1,MX2)
       MN1=MN1+MINEL
       MX1=MX1+MAXEL
       MN2=MN2+MAXEL
       MX2=MX2+MINEL
       MNR=MAX(MN1,MN2)
       MXR=MIN(MX1,MX2)
       NOREM=NGO
       DO ISPART=ISPSTA,ISPEND
        NO=MIN(NOREM,MORSBITS)
        NOREM=NOREM-NO
        LIMARR(1,ISPART)=MAX(MNL,MNR-NOREM)
        LIMARR(2,ISPART)=MIN(MXL+NGO-NOREM,MXR)
       END DO
      END DO
CTEST      write(*,*)' The LIMARR array:'
CTEST      write(*,'(1x,8i5)')(limarr(1,ispart),ispart=0,nasprt)
CTEST      write(*,'(1x,8i5)')(limarr(2,ispart),ispart=0,nasprt)

*----------------------------------------------------------------
* First, set up an array containing substring types.
* It would be nice to be able to quickly loop through just those
* substrings that have a limited population, and belong to a
* specific subpartition. To achieve this, we need a temporary
* array with index limits (0:KPOPMAX+1,1:NASPRT), and one with
* index limits (1:NSSTP).
      NSSTPTR=(KPOPMAX+2)*NASPRT
      CALL GETMEM('SSTPTR','ALLO','INTE',LSSTPTR,NSSTPTR)
      DO ISPART=1,NASPRT
        DO KPOP=0,KPOPMAX+1
          IWORK(LSSTPTR+KPOP+(KPOPMAX+2)*(ISPART-1))=0
        END DO
      END DO
      CALL GETMEM('SSTARR','ALLO','INTE',LSSTARR,NSSTP)
      DO ISST=1,NSSTP
       NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
       ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
       N=1+IWORK(LSSTPTR+NPOP+(KPOPMAX+2)*(ISPART-1))
       IWORK(LSSTPTR+NPOP+(KPOPMAX+2)*(ISPART-1))=N
      END DO
      IEND=0
      DO ISPART=1,NASPRT
        DO KPOP=0,KPOPMAX
          N=IWORK(LSSTPTR+KPOP+(KPOPMAX+2)*(ISPART-1))
          ISTA=IEND+1
          IEND=IEND+N
          IWORK(LSSTPTR+KPOP+(KPOPMAX+2)*(ISPART-1))=ISTA
        END DO
        IWORK(LSSTPTR+KPOPMAX+1+(KPOPMAX+2)*(ISPART-1))=IEND+1
      END DO
      DO ISST=1,NSSTP
       NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
       ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
       IADDR=IWORK(LSSTPTR+NPOP+(KPOPMAX+2)*(ISPART-1))
       IWORK(LSSTPTR+NPOP+(KPOPMAX+2)*(ISPART-1))=IADDR+1
       IWORK(LSSTARR-1+IADDR)=ISST
      END DO
      ISAVE=1
      DO ISPART=1,NASPRT
        DO KPOP=0,KPOPMAX
          NEXT=IWORK(LSSTPTR+KPOP+(KPOPMAX+2)*(ISPART-1))
          IWORK(LSSTPTR+KPOP+(KPOPMAX+2)*(ISPART-1))=ISAVE
          ISAVE=NEXT
        END DO
      END DO
CTESTC test print
CTEST      IADDR1=iwork(lsstptr+0+(kpopmax+2)*(1-1))
CTEST      IADDR2=iwork(lsstptr+kpopmax+1+(kpopmax+2)*(nasprt-1))
CTEST      write(*,'(1x,a,8i8)')'iaddr1,iaddr2:',iaddr1,iaddr2
CTEST      do IADDR=IADDR1,IADDR2-1
CTEST        isst=iwork(lsstarr-1+IADDR)
CTEST        kpop=ISSTAB(KSSTP+1+5*(ISST-1))
CTEST        ksym=ISSTAB(KSSTP+2+5*(ISST-1))
CTEST        kms2=ISSTAB(KSSTP+3+5*(ISST-1))
CTEST      write(*,'(1x,a,8i8)')'iaddr,isst,kpop,ksym,kms2:',
CTEST     &                      iaddr,isst,kpop,ksym,kms2
CTEST      end do
CTESTC End of test print

*----------------------------------------------------------------
* Here follows the construction of a graph, to be used for
* producing the FSB table.
* Use a single unique compound index for the vertices
* This way, we can use a temporary counter array to
* identify vertices.
* Subpartition is in 1..NASPRT, so we need levels 0..NASPRT
* Cumulative population is in 0..NACTEL
* Intermediate spin projection is in MS2MIN..MS2MAX:
      MS2MIN=(M2SPIN-KMS2MAX*NASPRT)/2
      MS2MAX=(M2SPIN+KMS2MAX*NASPRT)/2
* and it runs in steps of two units.
      MSFACT=1+(MS2MAX-MS2MIN)/2
* Intermediate symmetry label is in 1..NSYM
* All in all, the index
*  IND=ISCUM+NSYM*(
*    (M2CUM-MS2MIN)/2 + MSFACT*(NPCUM + (NACTEL+1)*ISPART))
* should do the trick.
* Run through the generation procedure twice. First time, we will
* determine the sizes of the vertex and downchain tables, and
* allocate them. Second time, we will fill in the values in these
* tables, and then the temporary counter table can be deallocated.
*----------------------------------
* Allocate temporary counter array:
      NVIDX=NSYM*MSFACT*(NACTEL+1)*(NASPRT+1)
      CALL GETMEM('VERIDX','ALLO','INTE',LVIDX,NVIDX)
      DO IVIDX=1,NVIDX
       IWORK(LVIDX-1+IVIDX)=0
      END DO
*----------------------------------
* Initialize top vertex:
      IPOS=LSYM+NSYM*((M2SPIN-MS2MIN)/2+MSFACT*(
     &              NACTEL+(NACTEL+1)*NASPRT))
      IWORK(LVIDX-1+IPOS)=1
      IVERT=1
*----------------------------------
* Recursively, find all other reachable vertices:
* Upper and lower level:
CTEST      CALL GETMEM('A','Check','Inte',LDUM,NDUM)
      DO LEVUP=NASPRT,1,-1
        LEVDWN=LEVUP-1
* A wasteful loop to pick out vertices on the upper level:
        DO NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
         DO M2C1=-NPC1,NPC1,2
          IF(M2C1.LT.MS2MIN) GOTO 130
          IF(M2C1.GT.MS2MAX) GOTO 130
          DO ISC1=1,NSYM
           IPOS=ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(
     &              NPC1+(NACTEL+1)*LEVUP))
           IVER1=IWORK(LVIDX-1+IPOS)
           IF(IVER1.EQ.0) GOTO 120
* This is a valid upper vertex on the upper level.
* Loop over valid arcs:
           KPOPLIM=MIN(NPC1,KPOPMAX)
           IADDR1=IWORK(LSSTPTR+0+(KPOPMAX+2)*(LEVUP-1))
           IADDR2=IWORK(LSSTPTR+KPOPLIM+1+(KPOPMAX+2)*(LEVUP-1))
           DO IADDR=IADDR1,IADDR2-1
             ISST=IWORK(LSSTARR-1+IADDR)
             KPOP=ISSTAB(KSSTP+1+5*(ISST-1))
             KSYM=ISSTAB(KSSTP+2+5*(ISST-1))
             KMS2=ISSTAB(KSSTP+3+5*(ISST-1))
* Just temporary sanity test:
             ISPART=ISSTAB(KSSTP+4+5*(ISST-1))
             IF(ISPART.NE.LEVUP) THEN
               WRITE(6,*)' THIS IS INSANE!'
               CALL ABEND()
             END IF
             NPC2=NPC1-KPOP
             IF(NPC2.LT.LIMARR(1,LEVDWN)) GOTO 110
             IF(NPC2.GT.LIMARR(2,LEVDWN)) GOTO 110
             M2C2=M2C1-KMS2
             IF(M2C2.LT.MS2MIN) GOTO 110
             IF(M2C2.GT.MS2MAX) GOTO 110
             ISC2=MUL(KSYM,ISC1)
             IF(ABS(M2C2).GT.NPC2) GOTO 110
             IF(NPC2.EQ.0 .AND. ISC2.NE.1) GOTO 110
* Inspect the lower vertex: Is it a new one?
             IPOS=ISC2+NSYM*((M2C2-MS2MIN)/2+MSFACT*(
     &              NPC2+(NACTEL+1)*LEVDWN))
             IVER2=IWORK(LVIDX-1+IPOS)
             IF(IVER2.GT.0) GOTO 110
* This is a new vertex. Register its number
             IVERT=IVERT+1
             IVER2=IVERT
             IWORK(LVIDX-1+IPOS)=IVER2
 110         CONTINUE
           END DO
* Finished looping over possible arcs.
 120      CONTINUE
         END DO
 130     CONTINUE
        END DO
       END DO
* Finished looping over possible upper vertices.
      END DO
* Finished loop over upper level, LEVUP.
*----------------------------------
* Now the index array holds a positive integer IV if the combination
* of properties is such that this vertex could be reached from above.
* We now need to keep just those that can also be reached from below.
* This is accomplished by a similar loop, now in the other direction:
*----------------------------------
CTEST      CALL GETMEM('B','Check','Inte',LDUM,NDUM)
      NARC=0
* Upper and lower level:
      DO LEVUP=1,NASPRT
        LEVDWN=LEVUP-1
* A wasteful loop to pick out vertices on the upper level:
        DO NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
         DO M2C1=-NPC1,NPC1,2
          IF(M2C1.LT.MS2MIN) GOTO 230
          IF(M2C1.GT.MS2MAX) GOTO 230
          DO ISC1=1,NSYM
           IPOS=ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(
     &              NPC1+(NACTEL+1)*LEVUP))
           IVER1=IWORK(LVIDX-1+IPOS)
           IF(IVER1.EQ.0) GOTO 220
* This is a valid upper vertex on the upper level.
* Is it reachable from below?
          NARCVRT=0
* Loop over valid arcs:
           KPOPLIM=MIN(NPC1,KPOPMAX)
           IADDR1=IWORK(LSSTPTR+0+(KPOPMAX+2)*(LEVUP-1))
           IADDR2=IWORK(LSSTPTR+KPOPLIM+1+(KPOPMAX+2)*(LEVUP-1))
           DO IADDR=IADDR1,IADDR2-1
             ISST=IWORK(LSSTARR-1+IADDR)
             KPOP=ISSTAB(KSSTP+1+5*(ISST-1))
             KSYM=ISSTAB(KSSTP+2+5*(ISST-1))
             KMS2=ISSTAB(KSSTP+3+5*(ISST-1))
* Just temporary sanity test:
             ISPART=ISSTAB(KSSTP+4+5*(ISST-1))
             IF(ISPART.NE.LEVUP) THEN
               WRITE(6,*)' THIS IS INSANE!'
               CALL ABEND()
             END IF
             NPC2=NPC1-KPOP
             IF(LEVDWN.EQ.0.AND.NPC2.GT.0) GOTO 210
             M2C2=M2C1-KMS2
             IF(M2C2.LT.MS2MIN) GOTO 210
             IF(M2C2.GT.MS2MAX) GOTO 210
             ISC2=MUL(KSYM,ISC1)
             IF(ABS(M2C2).GT.NPC2) GOTO 210
             IF(NPC2.EQ.0 .AND. ISC2.NE.1) GOTO 210
* Inspect the lower vertex: Is it reachable?
             IPOS=ISC2+NSYM*((M2C2-MS2MIN)/2+MSFACT*(
     &              NPC2+(NACTEL+1)*LEVDWN))
             IVER2=IWORK(LVIDX-1+IPOS)
             IF(IVER2.GT.0) NARCVRT=NARCVRT+1
 210         CONTINUE
           END DO
* Finished looping over possible arcs.
* Remove this vertex if it was unreachable.
          IF(NARCVRT.EQ.0) THEN
           IPOS=ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(
     &              NPC1+(NACTEL+1)*LEVUP))
           IWORK(LVIDX-1+IPOS)=0
          END IF
          NARC=NARC+NARCVRT
 220      CONTINUE
         END DO
 230     CONTINUE
        END DO
       END DO
* Finished looping over possible upper vertices.
      END DO
* Finished loop over upper level, LEVUP.
*----------------------------------
CTEST      CALL GETMEM('C','Check','Inte',LDUM,NDUM)
* Renumber the vertices:
      NVERT=0
* Upper and lower level:
      DO LEVUP=NASPRT,0,-1
* A wasteful loop to pick out vertices on the upper level:
        DO NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
         DO M2C1=-NPC1,NPC1,2
          IF(M2C1.LT.MS2MIN) GOTO 330
          IF(M2C1.GT.MS2MAX) GOTO 330
          DO ISC1=1,NSYM
           IPOS=ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(
     &              NPC1+(NACTEL+1)*LEVUP))
           IVER1=IWORK(LVIDX-1+IPOS)
           IF(IVER1.EQ.0) GOTO 320
           NVERT=NVERT+1
           IWORK(LVIDX-1+IPOS)=NVERT
 320       CONTINUE
          END DO
 330      CONTINUE
         END DO
        END DO
       END DO
*----------------------------------
* Now allocate the vertex and downchain tables, and loop again:
      NVERTAB=6*NVERT
      CALL GETMEM('VERTAB','ALLO','INTE',LVERTAB,NVERTAB)
      NDWNTAB=3*NARC
      CALL GETMEM('DWNTAB','ALLO','INTE',LDWNTAB,NDWNTAB)
      IARC =0
CTEST      write(*,'(1x,a,8i8)')'Nr of vertices NVERT=',NVERT
CTEST      write(*,'(1x,a,8i8)')'Nr of arcs     NARC =',NARC
* --------------------------------------------------------------
C Allocate LTDA(1:2,1:NASPRT), Level-to-Downarc array:
      CALL GETMEM('LTDA','Allo','Inte',LLTDA,2*(NASPRT+1))
*----------------------------------
* Upper and lower level:
CTEST      CALL GETMEM('D','Check','Inte',LDUM,NDUM)
      DO LEVUP=NASPRT,1,-1
        LEVDWN=LEVUP-1
* Level-to-Downarc array:
        IWORK(LLTDA+2*LEVUP)=IARC+1
        IWORK(LLTDA+2*LEVUP+1)=0
* Initialize counter of arcs from this upper level:
        NARCLEV=0
* A wasteful loop to pick out vertices on the upper level:
        DO NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
         DO M2C1=-NPC1,NPC1,2
          IF(M2C1.LT.MS2MIN) GOTO 430
          IF(M2C1.GT.MS2MAX) GOTO 430
          DO ISC1=1,NSYM
           IPOS=ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(
     &              NPC1+(NACTEL+1)*LEVUP))
           IVER1=IWORK(LVIDX-1+IPOS)
           IF(IVER1.EQ.0) GOTO 420
* This is a valid upper vertex on the upper level.
           IWORK(LVERTAB-1+1+6*(IVER1-1))=LEVUP
           IWORK(LVERTAB-1+2+6*(IVER1-1))=NPC1
           IWORK(LVERTAB-1+3+6*(IVER1-1))=M2C1
           IWORK(LVERTAB-1+4+6*(IVER1-1))=ISC1
           IWORK(LVERTAB-1+5+6*(IVER1-1))=IARC+1
           IWORK(LVERTAB-1+6+6*(IVER1-1))=0
* Initialize counter of arcs from this upper vertex:
           NARCVRT=0
* Loop over valid arcs:
           KPOPLIM=MIN(NPC1,KPOPMAX)
           IADDR1=IWORK(LSSTPTR+0+(KPOPMAX+2)*(LEVUP-1))
           IADDR2=IWORK(LSSTPTR+KPOPLIM+1+(KPOPMAX+2)*(LEVUP-1))
           DO IADDR=IADDR1,IADDR2-1
             ISST=IWORK(LSSTARR-1+IADDR)
             KPOP=ISSTAB(KSSTP+1+5*(ISST-1))
             KSYM=ISSTAB(KSSTP+2+5*(ISST-1))
             KMS2=ISSTAB(KSSTP+3+5*(ISST-1))
* Just temporary sanity test:
             ISPART=ISSTAB(KSSTP+4+5*(ISST-1))
             IF(ISPART.NE.LEVUP) THEN
               WRITE(6,*)' THIS IS INSANE!'
               CALL ABEND()
             END IF
             NPC2=NPC1-KPOP
             IF(LEVDWN.EQ.0.AND.NPC2.GT.0) GOTO 410
             M2C2=M2C1-KMS2
             IF(M2C2.LT.MS2MIN) GOTO 410
             IF(M2C2.GT.MS2MAX) GOTO 410
             ISC2=MUL(KSYM,ISC1)
             IF(ABS(M2C2).GT.NPC2) GOTO 410
             IF(NPC2.EQ.0 .AND. ISC2.NE.1) GOTO 410
* Is it a valid arc? See if a lower vertex exists.
             IPOS=ISC2+NSYM*((M2C2-MS2MIN)/2+MSFACT*(
     &              NPC2+(NACTEL+1)*LEVDWN))
             IVER2=IWORK(LVIDX-1+IPOS)
             IF(IVER2.EQ.0) GOTO 410
* A valid arc has been found. The upper vertex IVER1 is
* joined to the lower vertex IVER2 by an arc associated
* with the substring type ISST.
             NARCVRT=NARCVRT+1
             NARCLEV=NARCLEV+1
             IARC=IARC+1
             IWORK(LDWNTAB-1+1+3*(IARC-1))=ISST
             IWORK(LDWNTAB-1+2+3*(IARC-1))=IVER1
             IWORK(LDWNTAB-1+3+3*(IARC-1))=IVER2
 410         CONTINUE
           END DO
* Finished looping over possible arcs.
          IWORK(LVERTAB-1+6+6*(IVER1-1))=NARCVRT
 420      CONTINUE
         END DO
 430     CONTINUE
        END DO
       END DO
* Finished looping over vertices on this level.
       IWORK(LLTDA+2*LEVUP+1)=NARCLEV
      END DO
* Finished loop over upper level, LEVUP.
      IWORK(LVERTAB-1+1+6*(NVERT-1))=0
      IWORK(LVERTAB-1+2+6*(NVERT-1))=0
      IWORK(LVERTAB-1+3+6*(NVERT-1))=0
      IWORK(LVERTAB-1+4+6*(NVERT-1))=1
      IWORK(LVERTAB-1+5+6*(NVERT-1))=IARC
      IWORK(LVERTAB-1+6+6*(NVERT-1))=0
      IWORK(LLTDA)=IARC
      IWORK(LLTDA+1)=0
CTEST      CALL GETMEM('E','Check','Inte',LDUM,NDUM)
*---------------------------------------------------------
* The graph is finished. we do no longer need these arrays:
      CALL GETMEM('SSTPTR','FREE','INTE',LSSTPTR,NSSTPTR)
      CALL GETMEM('SSTARR','FREE','INTE',LSSTARR,NSSTP)
      CALL GETMEM('VERIDX','FREE','INTE',LVIDX,NVIDX)
*---------------------------------------------------------
* A graph has been constructed. We may construct a weight array:
      CALL GETMEM('WEIGHT','ALLO','INTE',LWEIGHT,NVERT)
      DO IVERT=1,NVERT-1
       IWORK(LWEIGHT-1+IVERT)=0
      END DO
      IWORK(LWEIGHT-1+NVERT)=1
CTEST      CALL GETMEM('G','Check','Inte',LDUM,NDUM)
      DO LEVUP=1,NASPRT
        LEVDWN=LEVUP-1
* Loop over arcs leading down from upper level:
        NARCLEV=IWORK(LLTDA+2*LEVUP+1)
        DO IA=1,NARCLEV
          IARC=IWORK(LLTDA+2*LEVUP)-1+IA
          ISST =IWORK(LDWNTAB-1+1+3*(IARC-1))
          IVER1=IWORK(LDWNTAB-1+2+3*(IARC-1))
          IVER2=IWORK(LDWNTAB-1+3+3*(IARC-1))
          ISUM=IWORK(LWEIGHT-1+IVER1)+IWORK(LWEIGHT-1+IVER2)
          IWORK(LWEIGHT-1+IVER1)=ISUM
        END DO
      END DO
C No more use for Level-to-Downarc array:
      CALL GETMEM('LTDA','Free','Inte',LLTDA,2*(NASPRT+1))
CTEST      CALL GETMEM('H','Check','Inte',LDUM,NDUM)
*---------------------------------------------------------
* Now, for any vertex iv, iwork(lweight-1+iv) is the total number of
* downwalks that can reach vertex nr nvert (The bottom vertex).
* Traverse the graph. Use a swich array. The graph is traversed
* from the top. At each vertex, one of the available arcs downwards
* is selected. Which arc to select is given by ISWITCH(LEVUP), where
* LEVUP is the level of the upper vertex.
      CALL GETMEM('Switch','Allo','Inte',LSWITCH,NASPRT)
* Initial values: Lowest walk
      DO LEVUP=1,NASPRT
        IWORK(LSWITCH-1+LEVUP)=1
      END DO
* The total number of walks in the graph:
      NFSB0=IWORK(LWEIGHT)
      NFSBARR=(NASPRT+2)*NFSB0
      CALL GETMEM('FSBARR','ALLO','INTE',LFSBARR,NFSBARR)
      NRDETS0=0
      NFSB=0
      NRDETS=0
CTEST      write(*,*)' A list of all generated FS blocks:'
 500  CONTINUE
* Construct this walk. While constructing it, also find out if
* it has a successor.
      IVUP=1
      LOWEST=NASPRT+1
      DO LEVUP=NASPRT,1,-1
       ISW=IWORK(LSWITCH-1+LEVUP)
* Select arc nr isw from those available to vertex ivup.
       IARC=IWORK(LVERTAB-1+5+6*(IVUP-1))-1+ISW
* Could it be incremented?
       NSW=IWORK(LVERTAB-1+6+6*(IVUP-1))
CTEST      write(*,'(1x,a,8i8)')'IVUP,ISW,NSW:',IVUP,ISW,NSW
       IF(NSW.GT.ISW) LOWEST=LEVUP
* Consult the downarc table:
       ISST =IWORK(LDWNTAB-1+1+3*(IARC-1))
       IVER1=IWORK(LDWNTAB-1+2+3*(IARC-1))
       IVER2=IWORK(LDWNTAB-1+3+3*(IARC-1))
CTEST      write(*,'(1x,a,8i8)')'ISST,IVER1,IVER2:',ISST,IVER1,IVER2
       IF(IVER1.NE.IVUP) THEN
         WRITE(6,*)' OOOOPS! THIS SHOULD NOT HAPPEN.'
         CALL ABEND()
       END IF
       ITRY(LEVUP)=ISST
       IVUP=IVER2
      END DO
* Check GAS restrictions:
      NSDBLK=1
      ISPEND=0
      KEEP=1
      DO IPART=1,NPART
       NGO  =NGASORB(0,IPART)
       NSP=(NGO+MORSBITS-1)/MORSBITS
       ISPSTA=ISPEND+1
       ISPEND=ISPEND+NSP
       NOREM=NGO
       ISUM=0
       DO ISPART=ISPSTA,ISPEND
        ISST=ITRY(ISPART)
        NSBS   = ISSTAB(KSSTP+0+5*(ISST-1))
        NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
        ISUM=ISUM+NPOP
        NSDBLK=NSDBLK*NSBS
       END DO
       IF(ISUM.LT.NGASLIM(1,IPART)) KEEP=0
       IF(ISUM.GT.NGASLIM(2,IPART)) KEEP=0
      END DO
      NRDETS0=NRDETS0+NSDBLK
      IF(KEEP.EQ.1) THEN
        NFSB=NFSB+1
        DO ISPART=1,NASPRT
         IWORK(LFSBARR-1+ISPART+(NASPRT+2)*(NFSB-1))=ITRY(ISPART)
        END DO
        IWORK(LFSBARR+NASPRT+(NASPRT+2)*(NFSB-1))=NSDBLK
        IWORK(LFSBARR+NASPRT+1+(NASPRT+2)*(NFSB-1))=NRDETS+1
        NRDETS=NRDETS+NSDBLK
CTEST      write(*,'(1x,i8,5x,10i5)') NFSB,(IWORK(LFSBARR-1+
CTEST     &                  LEVUP+NASPRT*(NFSB-1)),LEVUP=1,NASPRT),NSDBLK

* Next walk: Lowest increasable switch value is at level lowest.
      END IF
      IF(LOWEST.GT.NASPRT) GOTO 600
      DO LEVUP=1,LOWEST-1
       IWORK(LSWITCH-1+LEVUP)=1
      END DO
      IWORK(LSWITCH-1+LOWEST)=1+IWORK(LSWITCH-1+LOWEST)
      GOTO 500

 600  CONTINUE
      CALL GETMEM('Switch','Free','Inte',LSWITCH,NASPRT)
      CALL GETMEM('VERTAB','FREE','INTE',LVERTAB,NVERTAB)
      CALL GETMEM('DWNTAB','FREE','INTE',LDWNTAB,NDWNTAB)
      CALL GETMEM('WEIGHT','FREE','INTE',LWEIGHT,NVERT)
* Here we are. No more FS blocks can be constructed.
CTEST      write(*,*)' Total nr of Slater determinants, NRDETS=',NRDETS
CTEST      write(*,*)' Returning from VERTAB.'

      RETURN
      END
