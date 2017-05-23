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
      INTEGER FUNCTION NEWCNFTAB(NEL,NORB,MINOP,MAXOP,LSYM,NGAS,
     &                           NGASORB,NGASLIM,IFORM)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "symmul.fh"
      DIMENSION NGASLIM(2,NGAS)
      DIMENSION NGASORB(NSYM,NGAS)

      NEWCNFTAB=-1
CTEST      write(*,*)' Just entered NEWCNFTAB.'
CTEST      write(*,'(1x,a,10i8)')'NEL,NORB,MINOP,MAXOP:',NEL,NORB,MINOP,MAXOP
CTEST      write(*,'(1x,a,10i8)')'NSYM,LSYM,NGAS:',NSYM,LSYM,NGAS
C Note how input parameter LSYM is used: If non-zero, only those configurations
C with symmetry label LSYM are selected. But if LSYM=0, they are all selected.
C We must figure out sizes before allocating the new configuration table.
C Set up a table NCNF1(NSYM,NPOS) with NPOS=((NEL+1)*(NEL+2))/2
      NNCNF1=NSYM*((NEL+1)*(NEL+2))/2
      CALL GETMEM('NCnf1','Allo','Inte',LNCNF1,NNCNF1)
C We need also a table NCNF2, temporarily. Need to know mx nr of active orbitals
C in any GAS subspace:
      MXO=0
       DO IGAS=1,NGAS
        ISUM=0
        DO ISYM=1,NSYM
         ISUM=ISUM+NGASORB(ISYM,IGAS)
        END DO
        MXO=MAX(MXO,ISUM)
      END DO
      NNCNF2=NSYM*((MXO+1)*(MXO+2))/2
      CALL GETMEM('NCnf2','Allo','Inte',LNCNF2,NNCNF2)
      CALL NRCNF1(NEL,NORB,NGAS,NGASLIM,NGASORB,IWORK(LNCNF1),
     &            MXO,IWORK(LNCNF2))
      CALL GETMEM('NCnf2','Free','Inte',LNCNF2,NNCNF2)
C NCNF1(ISYM,IPOS) contains the number of possible configurations for symmetry
C label ISYM, nr of closed shells NCLS, and nr of open shells NOPN. The latter
C are combined as pair index IPOS=1+NOPN+(NOCC*(NOCC+1))/2 with NOCC=NCLS+NOPN

C Header (See below for contents):
      NTAB=10
C NGASORB array:
      NTAB=NTAB+(NSYM+1)*(NGAS+1)
C NGASLIM array:
      NTAB=NTAB+2*NGAS
C Save offset to INFO table for later use:
      KINFO=NTAB+1
C INFO array:
      NTAB=NTAB+3*NSYM*(MAXOP-MINOP+1)
C Save offset to configuration arrays for later use:
      KCNFSTA=NTAB+1
C Configuration arrays:
      DO NOPN=MINOP,MIN(2*NORB-NEL,NEL,MAXOP)
       NCLS=(NEL-NOPN)/2
       IF(NCLS.LT.0) GOTO 19
       IF(2*NCLS+NOPN.NE.NEL) GOTO 19
       NOCC=NCLS+NOPN
       IF(NOCC.GT.NORB) GOTO 19
       DO ISYM=1,NSYM
        NCNF=0
        IF(LSYM.GE.1 .AND. LSYM.LE.NSYM) THEN
          IPOS=1+NOPN+(NOCC*(NOCC+1))/2
          NCNF=IWORK(LNCNF1-1+ISYM+NSYM*(IPOS-1))
        END IF
        LENCNF=NOCC
        IF(IFORM.EQ.2) LENCNF=NORB
        IF(IFORM.EQ.3) LENCNF=(NOCC+3)/4
        IF(IFORM.EQ.4) LENCNF=(NORB+14)/15
        NTAB=NTAB+NCNF*LENCNF
       END DO
  19   CONTINUE
      END DO
C Sizes and offsets are known. Now, we can allocate the table:
      CALL GETMEM('CnfTab','Allo','Inte',LTAB,NTAB)
      NEWCNFTAB=LTAB
C Enter header:
      IWORK(LTAB+ 0)=NTAB
      IWORK(LTAB+ 1)=37
      IWORK(LTAB+ 2)=NEL
      IWORK(LTAB+ 3)=NORB
      IWORK(LTAB+ 4)=MINOP
      IWORK(LTAB+ 5)=MAXOP
      IWORK(LTAB+ 6)=NSYM
      IWORK(LTAB+ 7)=LSYM
      IWORK(LTAB+ 8)=NGAS
      IWORK(LTAB+ 9)=IFORM
C Enter copy of NGASORB array:
      DO IGAS=1,NGAS
       ISUM=0
       DO ISYM=1,NSYM
        L=LTAB+10+ISYM+(NSYM+1)*IGAS
        IWORK(L)=NGASORB(ISYM,IGAS)
        ISUM=ISUM+NGASORB(ISYM,IGAS)
       END DO
       IWORK(LTAB+10+(NSYM+1)*IGAS)=ISUM
      END DO
      DO ISYM=0,NSYM
       ISUM=0
       DO IGAS=1,NGAS
        L=LTAB+10+ISYM+(NSYM+1)*IGAS
        ISUM=ISUM+IWORK(L)
       END DO
       IWORK(LTAB+10+ISYM)=ISUM
      END DO
      L=LTAB+9+(NSYM+1)*(NGAS+1)
C Enter copy of NGASLIM array:
      DO IGAS=1,NGAS
       L=L+1
       IWORK(L)=NGASLIM(1,IGAS)
       L=L+1
       IWORK(L)=NGASLIM(2,IGAS)
      END DO
C Construct and enter INFO table.
      LINFO=LTAB-1+KINFO
C The INFO table has a relative pointer to configuration arrays:
      KCNFEND=KCNFSTA-1
      DO NOPN=MINOP,MAXOP
       NCLS=(NEL-NOPN)/2
       IFPOSS=1
       IF(NCLS.LT.0) IFPOSS=0
       IF(2*NCLS+NOPN.NE.NEL) IFPOSS=0
       IF(NCLS+NOPN.GT.NORB) IFPOSS=0
       IF(IFPOSS.EQ.0) THEN
        DO ISYM=1,NSYM
C No such configuration is possible.
C INFO(1,ISYM,NOPN)=NCNF
         IWORK(LINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))=0
C INFO(2,ISYM,NOPN)=NTAB+1
         IWORK(LINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))=-1
C INFO(3,ISYM,NOPN)=LENCNF
         IWORK(LINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))=0
CTEST      write(*,'(1x,a,2I5)')' Test NOPN,ISYM:',NOPN,ISYM
CTEST      write(*,*)' Impossible configuration.'
CTEST      write(*,*)' Info array:'
CTEST      write(*,'(1x,3I8)')(IWORK(LINFO+I+
CTEST     &                   3*(ISYM-1+NSYM*(NOPN-MINOP))),I=0,2)
        END DO
       ELSE
        NOCC=NCLS+NOPN
        DO ISYM=1,NSYM
         NCNF=0
C If LSYM=0, all symmetry labels will be accepted. Else, only
C those with ISYM=LSYM.
         IF(LSYM.EQ.0 .OR. ISYM.EQ.LSYM) THEN
           IPOS=1+NOPN+(NOCC*(NOCC+1))/2
           NCNF=IWORK(LNCNF1-1+ISYM+NSYM*(IPOS-1))
         END IF
         IF(NCNF.EQ.0) THEN
C INFO(1,ISYM,NOPN)=NCNF
           IWORK(LINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))=0
C INFO(2,ISYM,NOPN)=NTAB+1
           IWORK(LINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))=-1
C INFO(3,ISYM,NOPN)=LENCNF
           IWORK(LINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))=0
         ELSE
           LENCNF=NOCC
           IF(IFORM.EQ.2) LENCNF=NORB
           IF(IFORM.EQ.3) LENCNF=(NOCC+3)/4
           IF(IFORM.EQ.4) LENCNF=(NORB+14)/15
C The relative pointer into this configuration array:
           KCNFSTA=KCNFEND+1
           KCNFEND=KCNFEND+NCNF*LENCNF
C INFO(1,ISYM,NOPN)=NCNF
           IWORK(LINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))=NCNF
C INFO(2,ISYM,NOPN)=NTAB+1
           IWORK(LINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))=KCNFSTA
C INFO(3,ISYM,NOPN)=LENCNF
           IWORK(LINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))=LENCNF
CTEST      write(*,'(1x,a,2I5)')' Test NOPN,ISYM:',NOPN,ISYM
CTEST      write(*,*)' Possible configuration.'
CTEST      write(*,*)' Info array:'
CTEST      write(*,'(1x,3I8)')(IWORK(LINFO+I+
CTEST     &                   3*(ISYM-1+NSYM*(NOPN-MINOP))),I=0,2)
           DO I=1,NCNF*LENCNF
             IWORK(LTAB-1+KCNFSTA-1+I)=0
           END DO
         END IF
        END DO
       END IF
      END DO
C The NCNF1 array is no longer needed.
      CALL GETMEM('NCnf1','Free','Inte',LNCNF1,NNCNF1)

C Finally, only now when we know where to store each (ISYM,NOPN) block of
C configurations, can we compute the actual configuration arrays:
      CALL MKCONF(IWORK(LTAB))

      RETURN
      END
      SUBROUTINE NRCNF1(MAXEL,NORB,NGAS,NGASLIM,
     &                  NGASORB,NCNF1,MXTMP,NCNF2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "symmul.fh"
      DIMENSION NCNF1( NSYM, ((MAXEL+1)*(MAXEL+2))/2 )
      DIMENSION NGASLIM(2,NGAS),NGASORB(NSYM,NGAS)
      DIMENSION NCNF2(NSYM, ((MXTMP+1)*(MXTMP+2))/2 )
C Returns the array NCNF1, which contains the number of
C configurations with the following criteria:
C    Orbital indices range from 1..NORB
C    GAS restrictions described by NGASLIM, NGASORB
C    Total number of electrons is at most MAXEL
C    NCLS closed-shell and NOPN open-shell orbitals
C    Symmetry label LSYM
C for all possible values of NCLS,NOPN and LSYM, stored as
C          NCNF1(LSYM,IPOS)
C with IPOS=(NOCC*(NOCC+1))/2+NOPN+1, NOCC=NCLS+NOPN,
C provided that 0<=NCLS, 0<=NOPN, and NOCC<=MIN(NORB,MAXEL).
C MXTMP=Max nr of orbitals in one GAS partition.
C Prerequisite: The orbital symmetry labels stored in ISM.
C               The GAS restriction arrays
C Method: Induction over GAS partitions.

      MAXOCC=MIN(MAXEL,NORB)
C Initialize:
      DO IPOS=1,((MAXOCC+1)*(MAXOCC+2))/2
       DO ISYM=1,NSYM
        NCNF1(ISYM,IPOS)=0
       END DO
      END DO
      NCNF1(1,1)=1
      CALL GETMEM('OrbSym','Allo','Inte',LISM,NORB)
C Max nr of occupied orbitals so far:
      NOCCMX=0
      DO IGAS=1,NGAS
        MXOCCOLD=NOCCMX
C Nr of orbitals in this partition
        NO=0
        II=0
        DO ISYM=1,NSYM
          NG=NGASORB(ISYM,IGAS)
          NO=NO+NG
          DO I=1,NG
           II=II+1
           IWORK(LISM-1+II)=ISYM
          END DO
        END DO
        NELMN=MAX(0,NGASLIM(1,IGAS))
        NELMX=MIN(2*NO,NGASLIM(2,IGAS))
C Possibilities within this partition:
CTEST      write(*,*)' In GAS partition nr IGAS=',IGAS
CTEST      write(*,*)' Prepare to call NRCNF2.'
CTEST      write(*,'(1x,a,8I5)')' Input data is: NO=',NO
CTEST      write(*,'(1x,a,8I5)')'              NSYM=',NSYM
CTEST      write(*,'(1x,a,8I5)')'               ISM=',
CTEST     &                       (IWORK(LISM-1+I),I=1,NO)
CTEST      write(*,*)' Call NRCNF2.'
        CALL NRCNF2(NO,IWORK(LISM),NCNF2)
CTEST      write(*,*)' Back from NRCNF2.'
CTEST      write(*,*)' The NCNF2 array:'
CTEST      ij=0
CTEST      do i=0,no
CTEST       do j=0,i
CTEST        ij=ij+1
CTEST        write(*,'(1x,2i5,5x,8i5)')i-j,j,(ncnf2(isym,ij),isym=1,nsym)
CTEST       end do
CTEST      end do

CTEST      write(*,*)' GAS partition nr IGAS=',IGAS
CTEST      write(*,*)' Before this partition, NCNF1 array is:'
CTEST      ij=0
CTEST      do i=0,noccmx
CTEST       do j=0,i
CTEST        ij=ij+1
CTEST        write(*,'(1x,2i5,5x,8i5)')i-j,j,(ncnf1(isym,ij),isym=1,nsym)
CTEST       end do
CTEST      end do

C Occupied and open within earlier+new partition:
CTEST        write(*,*)' New NOCCNW up to:',MIN(MAXOCC,MXOCCOLD+NELMX)
        DO NOCCNW=MIN(MAXOCC,MXOCCOLD+NELMX),0,-1
        DO NOPNNW=0,NOCCNW
         NCLSNW=NOCCNW-NOPNNW
         IPOSNW=(NOCCNW*(NOCCNW+1))/2+NOPNNW+1
         DO ISYMNW=1,NSYM
          NEW=0
CTEST        write(*,'(1x,a,3I4)')' NCLSNW,NOPNNW,ISYMNW:',
CTEST     &                                    NCLSNW,NOPNNW,ISYMNW
CTEST        write(*,*)' Add the following contributions:'
CTEST        write(*,*)' Limits on NOCC:',NELMN/2,MIN(NELMX,NO)
C Occupied and open within the GAS partition:
          DO NOCC=NELMN/2,MIN(NELMX,NO)
CTEST        write(*,*)' Limits on NOPN:',MAX(0,2*NOCC-NELMX),
CTEST     &                    MIN(2*NOCC-NELMN,NOCC,NELMX)
          DO NOPN=MAX(0,2*NOCC-NELMX),
     &                    MIN(2*NOCC-NELMN,NOCC,NELMX)
           NCLS=NOCC-NOPN
           IPOS=(NOCC*(NOCC+1))/2+NOPN+1
           DO ISYM=1,NSYM
            NY=NCNF2(ISYM,IPOS)
            IF(NY.EQ.0) GOTO 19
            NCLSOLD=NCLSNW-NCLS
            IF(NCLSOLD.LT.0) GOTO 19
            NOPNOLD=NOPNNW-NOPN
            IF(NOPNOLD.LT.0) GOTO 19
            NOCCOLD=NCLSOLD+NOPNOLD
            IF(NOCCOLD.GT.MXOCCOLD) GOTO 19
            ISYMOLD=MUL(ISYM,ISYMNW)
            IPOSOLD=(NOCCOLD*(NOCCOLD+1))/2+NOPNOLD+1
            NX=NCNF1(ISYMOLD,IPOSOLD)
            IF(NX.EQ.0) GOTO 19
C Test print:
CTEST        write(*,'(1x,a,5I4)')' NCLS,NOPN,ISYM,NCNF1,NCNF2:',
CTEST     &                            NCLS,NOPN,ISYM,NX,NY
C End of test prints
            NEW=NEW+NCNF1(ISYMOLD,IPOSOLD)*NCNF2(ISYM,IPOS)
            NOCCMX=MAX(NOCCNW,NOCCMX)
  19        CONTINUE
           END DO
          END DO
          END DO
          NCNF1(ISYMNW,IPOSNW)=NEW
         END DO
        END DO
        END DO

CTEST      write(*,*)' GAS partition nr IGAS=',IGAS
CTEST      write(*,*)' After this partition, NCNF1 array is:'
CTEST      ij=0
CTEST      do i=0,noccmx
CTEST       do j=0,i
CTEST        ij=ij+1
CTEST        write(*,'(1x,2i5,5x,8i5)')i-j,j,(ncnf1(isym,ij),isym=1,nsym)
CTEST       end do
CTEST      end do

      END DO
      CALL GETMEM('OrbSym','Free','Inte',LISM,NORB)
      RETURN
      END
      SUBROUTINE NRCNF2(NORB,ISM,NCNF2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "symmul.fh"
      DIMENSION NCNF2(NSYM, ((NORB+1)*(NORB+2))/2 )
      DIMENSION ISM(NORB)
C Returns the array NCNF2, which contains the number of
C (sub-)configurations with NCLS closed-shell and NOPN open-shell
C orbitals and having symmetry label LSYM, stored as
C          NCNF2(LSYM,IPOS)
C with IPOS=(NOCC*(NOCC+1))/2+NOPN+1, NOCC=NCLS+NOPN,
C provided that 0<=NCLS, 0<=NOPN, and NOCC<=NORB.
C Prerequisite: The orbital symmetry labels stored in ISM.
C Method: Induction

      DO IPOS=1,((NORB+1)*(NORB+2))/2
       DO ISYM=1,NSYM
        NCNF2(ISYM,IPOS)=0
       END DO
      END DO
      NCNF2(1,1)=1
      DO L=1,NORB
CTEST      write(*,*)' NRCNF2 level L=',L
CTEST      write(*,*)'  ISM(L)=',ISM(L)
        DO NOCC=L,1,-1
          DO NOPN=0,NOCC
            NCLS=NOCC-NOPN
CTEST      write(*,'(1x,a,8i5)')' NCLS,NOPN=',NCLS,NOPN
            IPOS1=(NOCC*(NOCC+1))/2+NOPN+1
            IPOS2=IPOS1-NOCC
            IPOS3=IPOS2-1
CTEST      write(*,'(1x,a,8i5)')' IPOS1,IPOS2,IPOS3:',IPOS1,IPOS2,IPOS3
CTEST      write(*,'(1x,a,8i5)')
CTEST     &         ' From last level:',(NCNF2(I,IPOS1),I=1,NSYM)
CTEST      if(ncls.gt.0) then
CTEST      write(*,'(1x,a,8i5)')
CTEST     &         ' Adding closed  :',(NCNF2(I,IPOS2),I=1,NSYM)
CTEST      end if
CTEST      if(nopn.gt.0) then
CTEST      write(*,'(1x,a,8i5)')
CTEST     &   ' Adding open    :',(NCNF2(MUL(ISM(L),I),IPOS3),I=1,NSYM)
CTEST      write(*,'(1x,a,8i5)')
CTEST     &   ' (Descrambled)  :',(NCNF2(I,IPOS3),I=1,NSYM)
CTEST      end if
            DO ISYM=1,NSYM
              NEW=NCNF2(ISYM,IPOS1)
              IF(NCLS.GT.0) NEW=NEW+NCNF2(ISYM,IPOS2)
              JSYM=MUL(ISM(L),ISYM)
              IF(NOPN.GT.0) NEW=NEW+NCNF2(JSYM,IPOS3)
              NCNF2(ISYM,IPOS1)=NEW
            END DO
CTEST      write(*,'(1x,a,8i5)')
CTEST     &         ' New result     :',(NCNF2(I,IPOS1),I=1,NSYM)
          END DO
        END DO
      END DO
      RETURN
      END
      SUBROUTINE MKCONF(ICNFTAB)
      IMPLICIT NONE

      INTEGER ICNFTAB(*)
      INTEGER NEL,MINOP,MAXOP,LSYM
      INTEGER MXPRT
      PARAMETER (MXPRT=150)
      INTEGER LIMPOP(2,MXPRT),LIMOP(2,MXPRT),LIMCL(2,MXPRT)
      INTEGER IOPDST(MXPRT),ICLDST(MXPRT),IOC(MXPRT),ICNF(MXPRT)
      INTEGER LIM1,LIM2,LIM1SUM,LIM2SUM,IGAS,NOR,IERR
      INTEGER MNOP,MXOP,MNCL,MXCL,NOPN,NCLS
      INTEGER INIT1,INIT2,M,MORE,NOP1,NCL2
      INTEGER ISYM,NORB,ICONF
      INTEGER IFORM,IR,ITYPE,IW,KCNFSTA,KGASLIM
      INTEGER KGASORB,KINFO,KPOS,LCLS,LENCNF,LISM,LOPN,NCNF
      INTEGER NCNFSYM(8),NGAS,NOCC,NTAB
      INTEGER I,J,K,IOFF,IO,IORB,N,NCL,NOP
#include "WrkSpc.fh"
#include "symmul.fh"
C     INTEGER MIN,MAX
      INTRINSIC MIN,MAX
      INTEGER IPOW4(0:15),IPOW256(0:3)
      DATA IPOW4 / 1,4,16,64,256,1024,4096,16384,65536,262144,1048576,
     &             4194304,16777216,67108864,268435456,1073741824 /
      DATA IPOW256 / 1,256,65536,16777216 /

      ITYPE=ICNFTAB(2)
      IF(ITYPE.NE.37) THEN
        WRITE(6,*)'MKCONF error: This is not a configuration table!'
        CALL ABEND()
      END IF
C Unbutton the CNF table.
      NTAB =ICNFTAB(1)
      NEL  =ICNFTAB(3)
      NORB =ICNFTAB(4)
      MINOP=ICNFTAB(5)
      MAXOP=ICNFTAB(6)
      NSYM =ICNFTAB(7)
      LSYM =ICNFTAB(8)
      NGAS =ICNFTAB(9)
      IFORM=ICNFTAB(10)
      KGASORB=11
      KGASLIM=KGASORB+(NSYM+1)*(NGAS+1)
C Check and refine the GAS limits:
      IF(NGAS.LE.0 .OR. NGAS.GT.MXPRT) THEN
         WRITE(6,*)' MKCONF ERROR: Nr of GAS partitions is out of'
         WRITE(6,*)' bounds. NGAS must be .GT.0 and .LT. MXPRT=',MXPRT
         WRITE(6,*)' Input argument NGAS is ',NGAS
         CALL ABEND()
      END IF
      LIM1SUM=0
      LIM2SUM=0
      NORB=0
      DO IGAS=1,NGAS
       NOR=ICNFTAB(KGASORB+(NSYM+1)*IGAS)
       IF(NOR.LT.0) IERR=1
       LIM1=MAX(0,ICNFTAB(KGASLIM+2*(IGAS-1)))
       LIM2=MIN(NEL,ICNFTAB(KGASLIM+1+2*(IGAS-1)),2*NOR)
       LIM1SUM=LIM1SUM+LIM1
       LIM2SUM=LIM2SUM+LIM2
       LIMPOP(1,IGAS)=LIM1
       LIMPOP(2,IGAS)=LIM2
       NORB=NORB+NOR
      END DO
      IERR=0
      DO IGAS=1,NGAS
       LIM1=MAX(LIMPOP(1,IGAS),NEL-(LIM2SUM-LIMPOP(2,IGAS)))
       LIM2=MIN(LIMPOP(2,IGAS),NEL-(LIM1SUM-LIMPOP(1,IGAS)))
       LIMPOP(1,IGAS)=LIM1
       LIMPOP(2,IGAS)=LIM2
       IF(LIM1.GT.LIM2) IERR=1
      END DO
      IF(IERR.GT.0) THEN
         WRITE(6,*)' MKCONF ERROR: The input GAS restrictions are'
         WRITE(6,*)' impossible to meet. No configurations are'
         WRITE(6,*)' generated. The program stops here.'
         WRITE(6,'(1X,A,I2)')' Number of GAS partitions:',NGAS
         WRITE(6,'(1X,A,50I3)')' Partition:',(IGAS,IGAS=1,NGAS)
         WRITE(6,'(1X,A,50I3)')' NGASORB:  ',
     &             (ICNFTAB(KGASORB+(NSYM+1)*IGAS),IGAS=1,NGAS)
         WRITE(6,'(1X,A,50I3)')'NGASLIM(1):',
     &              (ICNFTAB(KGASLIM  +2*(IGAS-1)),IGAS=1,NGAS)
         WRITE(6,'(1X,A,50I3)')'NGASLIM(2):',
     &              (ICNFTAB(KGASLIM+1+2*(IGAS-1)),IGAS=1,NGAS)
         WRITE(6,'(1X,A,I2)')' Number of electrons:',NEL
         CALL ABEND()
      END IF
      IF(NEL.LT.0 .OR. NEL.GT.2*NORB) THEN
         WRITE(6,*)' MKCONF ERROR: Nr of electrons is out of bounds.'
         WRITE(6,*)' NEL must be .GT.0 and .LT. 2*NORB=',2*NORB
         WRITE(6,*)' Input argument NEL is ',NEL
         CALL ABEND()
      END IF
C Array for orbital symmetry:
      CALL GETMEM('OrbSym','Allo','Inte',LISM,NORB)
C Initialize table with orbital symmetry.
      IORB=0
      DO IGAS=1,NGAS
       DO ISYM=1,NSYM
        N=ICNFTAB(KGASORB+ISYM+(NSYM+1)*IGAS)
        DO K=1,N
         IORB=IORB+1
         IWORK(LISM-1+IORB)=ISYM
        END DO
       END DO
      END DO
C INFO table inside ICNFTAB:
      KINFO=KGASLIM+2*NGAS
C Note: Nr of conf, their position and length can now be accessed as:
C      NCNF   =ICNFTAB(KINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))
C      KCNFSTA=ICNFTAB(KINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
C      LENCNF =ICNFTAB(KINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))
C Make a list of possible number of open shells in each partition:
      DO IGAS=1,NGAS
       NOR=ICNFTAB(KGASORB+(NSYM+1)*IGAS)
       MNOP=1
       MXOP=0
       DO I=LIMPOP(1,IGAS),LIMPOP(2,IGAS)
        MNOP=MIN(MNOP,MOD(I,2))
        MXOP=MAX(MXOP,MIN(I,2*NOR-I))
       END DO
       LIMOP(1,IGAS)=MNOP
       LIMOP(2,IGAS)=MXOP
      END DO
C Counter of configurations:
      ICONF=0
C Loop over the requested range of open shells:
      DO NOPN=MINOP,MAXOP
       NCLS=(NEL-NOPN)/2
       IF(NCLS.LT.0) GOTO 120
       IF(2*NCLS+NOPN.NE.NEL) GOTO 120
C Size of each entry in the configuration table:
       NOCC=NCLS+NOPN
       LENCNF=NOCC
       IF(IFORM.EQ.2) LENCNF=NORB
       IF(IFORM.EQ.3) LENCNF=(NOCC+3)/4
       IF(IFORM.EQ.4) LENCNF=(NORB+14)/15
C Counter of configurations/symmetry for this nr of open shells:
       DO ISYM=1,NSYM
        NCNFSYM(ISYM)=0
       END DO
C Loop over all ways of distributing NOPN open shells among
C the partitions. First make a start distribution:
       INIT1=NGAS
       NOP1=NOPN
  10   CONTINUE
C Create the lexically lowest distribution with NOP1 open
C shells among the INIT1 lowest partitions, and
C increment the next higher partition (if any).
C Let M=Max tot nr of open shells in lower partitions.
       IF(INIT1.LT.NGAS) IOPDST(INIT1+1)=IOPDST(INIT1+1)+1
       M=NOP1
       DO IGAS=INIT1,1,-1
        M=M-LIMOP(1,IGAS)
       END DO
       IF(M.LT.0) GOTO 120
C But actually, we start with zero. So M open shells must be
C distributed in excess of the allowed minimum, among the INIT1
C partitions.
       DO IGAS=1,INIT1
        MORE=MIN(LIMOP(2,IGAS)-LIMOP(1,IGAS),M)
        IOPDST(IGAS)=LIMOP(1,IGAS)+MORE
        M=M-MORE
       END DO
       IF(M.GT.0) GOTO 120
C At this point of the code, all possible distributions of
C open shells will be generated. Use them.
CTEST      write(*,*)' Try this IOPDST distribution:'
CTEST      write(*,'(1x,a,50i3)')
CTEST     &          ' IOPDST(IGAS):',(IOPDST(IGAS),IGAS=1,NGAS)
C First, use it to generate a table of limits for the distribution
C of closed shells:
      DO IGAS=1,NGAS
       NOP=IOPDST(IGAS)
       NOR=ICNFTAB(KGASORB+(NSYM+1)*IGAS)-NOP
       LIM1=MAX(0,(LIMPOP(1,IGAS)-NOP)/2)
       LIM2=MIN(NOR,(LIMPOP(2,IGAS)-NOP)/2)
       IF(LIM1.GT.LIM2) GOTO 110
       MNCL=LIM2
       MXCL=LIM1
       DO I=LIM1,LIM2
        N=2*I+NOP
        IF(N.GE.LIMPOP(1,IGAS) .AND. N.LE.LIMPOP(2,IGAS)) THEN
         MNCL=MIN(MNCL,I)
         MXCL=MAX(MXCL,I)
        END IF
       END DO
       IF(MNCL.GT.MXCL) GOTO 110
       LIMCL(1,IGAS)=MNCL
       LIMCL(2,IGAS)=MXCL
      END DO
CTEST      write(*,'(1x,a,20i3)')
CTEST     &          'LIMCL(1,IGAS):',(LIMCL(1,IGAS),IGAS=1,NGAS)
CTEST      write(*,'(1x,a,50i3)')
CTEST     &          'LIMCL(2,IGAS):',(LIMCL(2,IGAS),IGAS=1,NGAS)

C Loop over all possible ways of distributing NCLS closed shells
C among  the partitions, subject to restrictions.
C In order to create the start distribution:
      INIT2=NGAS
      NCL2=(NEL-NOPN)/2
  20  CONTINUE
C Create the lexically lowest distribution with NCL2 closed shells
C among the INIT2 lowest partitions, and increment the next higher
C partition (if any).
      IF(INIT2.LT.NGAS) ICLDST(INIT2+1)=ICLDST(INIT2+1)+1
      M=NCL2
      DO IGAS=INIT2,1,-1
       M=M-LIMCL(1,IGAS)
      END DO
      IF(M.LT.0) GOTO 110
      DO IGAS=1,INIT2
       MORE=MIN(LIMCL(2,IGAS)-LIMCL(1,IGAS),M)
       ICLDST(IGAS)=LIMCL(1,IGAS)+MORE
       M=M-MORE
      END DO
      IF(M.GT.0) GOTO 110
C Here follows code to use this population distribution.
C Initialize the configuration subarrays of partitions nr
C 1..NGAS, within this population distribution:
CTEST      write(*,*)' New distribution of open/closed shells among'
CTEST      write(*,*)' the GAS spaces is:'
CTEST      write(*,'(1x,a,50i3)')
CTEST     &          ' IOPDST(IGAS):',(IOPDST(IGAS),IGAS=1,NGAS)
CTEST      write(*,'(1x,a,50i3)')
CTEST     &          ' ICLDST(IGAS):',(ICLDST(IGAS),IGAS=1,NGAS)
      IORB=0
      DO IGAS=1,NGAS
       NCL=ICLDST(IGAS)
       DO I=1,NCL
        IORB=IORB+1
        IOC(IORB)=2
       END DO
       NOP=IOPDST(IGAS)
       DO I=1,NOP
        IORB=IORB+1
        IOC(IORB)=1
       END DO
       NOR=ICNFTAB(KGASORB+(NSYM+1)*IGAS)
       DO I=1,NOR-NCL-NOP
        IORB=IORB+1
        IOC(IORB)=0
       END DO
      END DO
  30  CONTINUE
C Here finally we will get all possible configurations, restricted
C by the population arrays. Screening by combined symmetry:
CTEST      write(*,'(1x,a,50I3)')' New config:',(ioc(io),io=1,norb)
CTEST      write(*,*)' Where should it go??'
      ISYM=1
      DO IO=1,NORB
       IF(IOC(IO).EQ.1) ISYM=MUL(IWORK(LISM-1+IO),ISYM)
      END DO
CTEST      write(*,*)' Symmetry ISYM=',ISYM
C Skip if wrong symmetry:
      IF (LSYM.GT.0 .AND. ISYM.NE.LSYM) GOTO 99
CTEST      write(*,*)' Keep it!'
      ICONF=ICONF+1
CTEST      write(*,*)' This is configuration nr ICONF=',ICONF
C Put this configuration into the ICNFTAB table.
C First, determine where it should go:
      N=NCNFSYM(ISYM)
CTEST      write(*,*)' Earlier nr of conf in this symm & nopn:',N
      NCNFSYM(ISYM)=N+1
CTEST      write(*,*)' New value NCNFSYM(ISYM)=',NCNFSYM(ISYM)
      KCNFSTA=ICNFTAB(KINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
      KPOS=KCNFSTA+N*LENCNF
      IF(KPOS+LENCNF-1.GT.NTAB) THEN
        WRITE(6,*)' MKCONF error: Table overflow.'
        WRITE(6,*)'KCNFSTA:',KCNFSTA
        WRITE(6,*)' LENCNF:',LENCNF
        WRITE(6,*)'   KPOS:',KPOS
        WRITE(6,*)'   NTAB:',NTAB
        CALL ABEND()
      END IF
C Put together configuration array in standard format:
        LCLS=1
        LOPN=NCLS+1
        DO IO=1,NORB
         N=IOC(IO)
         IF(N.EQ.1) THEN
           ICNF(LOPN)=IO
           LOPN=LOPN+1
         ELSE IF(N.EQ.2) THEN
           ICNF(LCLS)=IO
           LCLS=LCLS+1
         END IF
        END DO
C Add this configuration to the configuration table:
CTEST      write(*,*)' Add configuration to conf table.'
        IF(IFORM.EQ.1) THEN
          DO I=1,NOCC
            ICNFTAB(KPOS-1+I)=ICNF(I)
          END DO
        ELSE IF(IFORM.EQ.3) THEN
          DO I=1,NOCC
           IW=(3+I)/4
           IR=(3+I)-4*IW
           IF(IR.EQ.0) THEN
            ICNFTAB(KPOS-1+IW)=ICNF(I)
           ELSE
            ICNFTAB(KPOS-1+IW)=ICNFTAB(KPOS-1+IW)+
     &                   IPOW256(IR)*ICNF(I)
           END IF
          END DO
        ELSE
          IF(IFORM.EQ.2) THEN
            DO I=1,NORB
              ICNFTAB(KPOS-1+I)=IOC(I)
            END DO
          ELSE IF(IFORM.EQ.4) THEN
            DO I=1,NORB
             IW=(14+I)/15
             IR=(14+I)-15*IW
             IF(IR.EQ.0) THEN
              ICNFTAB(KPOS-1+IW)=IOC(I)
             ELSE
              ICNFTAB(KPOS-1+IW)=ICNFTAB(KPOS-1+IW)+
     &                     IPOW4(IR)*IOC(I)

             END IF
            END DO
          END IF
        END IF

  99  CONTINUE
C Get next configuration.
CTEST      write(*,*)' 99 CONTINUE: Get next conf (if any).'
      IOFF=0
      DO IGAS=1,NGAS
       NOR=ICNFTAB(KGASORB+(NSYM+1)*IGAS)
CTEST      write(*,*)'   Is there a new permut in IGAS',IGAS
C Try to find next permutation within this partition:
       DO K=2,NOR
        IF(IOC(IOFF+K-1).GT.IOC(IOFF+K)) THEN
         DO I=1,(K-1)/2
          J=IOC(IOFF+I)
          IOC(IOFF+I)=IOC(IOFF+K-I)
          IOC(IOFF+K-I)=J
         END DO
         DO I=K-1,1,-1
          IF(IOC(IOFF+I).GT.IOC(IOFF+K)) THEN
           J=IOC(IOFF+I)
           IOC(IOFF+I)=IOC(IOFF+K)
           IOC(IOFF+K)=J
CTEST      write(*,*)'   Yes there is! GOTO 30.'
C OK, the next permutation has been obtained.
           GOTO 30
          END IF
         END DO
        END IF
       END DO
CTEST      write(*,*)'   No, there is not. Reset IGAS',IGAS
C Not possible. Reset permutation in this partition, and
C then try the next one:
       N=ICLDST(IGAS)
       M=IOPDST(IGAS)
       DO IO=1,N
        IOC(IOFF+IO)=2
       END DO
       DO IO=N+1,N+M
        IOC(IOFF+IO)=1
       END DO
       DO IO=N+M+1,NOR
        IOC(IOFF+IO)=0
       END DO
       IOFF=IOFF+NOR
      END DO
CTEST      write(*,*)' There is no more, with this distrib.'
C All failed. No more configurations with this distribution of
C closed and open shells.
C Next ICLDST distribution. First find the first increasable index:
      M=0
      NCL2=-1
      DO IGAS=1,NGAS
       INIT2=IGAS-1
       IF(M.GT.0 .AND. ICLDST(IGAS).LT.LIMCL(2,IGAS)) GOTO 20
       M=M+ICLDST(IGAS)-LIMCL(1,IGAS)
       NCL2=NCL2+ICLDST(IGAS)
      END DO
CTEST      write(*,*)' No more ICLDST distribution.'
C No more ICLDST distribution is possible.
 110  CONTINUE
C Next IOPDST distribution. First find the first increasable index:
C That is the first partition with less than LIMOP(2,IGAS) open
C shells, above partitions with nonzero excess number M.
      M=0
      NOP1=-1
      DO IGAS=1,NGAS
       INIT1=IGAS-1
       IF(M.GT.0 .AND. IOPDST(IGAS).LT.LIMOP(2,IGAS)) GOTO 10
       M=M+IOPDST(IGAS)-LIMOP(1,IGAS)
       NOP1=NOP1+IOPDST(IGAS)
      END DO
CTEST      write(*,*)' No more IOPDST distribution.'
CTEST      write(*,*)' There is no more distribution possible.'
CTEST      write(*,*)' Test NCNFSYM(ISYM) vs. INFO:'
CTEST      write(*,'(1x,a,8I5)')'ICNFTAB(KINFO+...):',
CTEST     &  (ICNFTAB(KINFO+3*(ISYM-1+NSYM*(NOPN-MINOP))),ISYM=1,NSYM)
CTEST      write(*,'(1x,a,8I5)')'     NCNFSYM(ISYM):',
CTEST     &                                (NCNFSYM(ISYM),ISYM=1,NSYM)
C Temporary check: Has everything worked perfectly??
      IERR=0
      DO ISYM=1,NSYM
        N=ICNFTAB(KINFO  +3*(ISYM-1+NSYM*(NOPN-MINOP)))
        IF(NCNFSYM(ISYM).NE.N) IERR=1
      END DO
      IF(IERR.NE.0) GOTO 900
C No more IOPDST distribution is possible. Next NOPN value:
 120  CONTINUE
      END DO

      CALL GETMEM('OrbSym','Free','Inte',LISM,NORB)
      RETURN
 900  CONTINUE
      WRITE(6,*)' MKCNF ERROR: Unforeseen calamity.'
      WRITE(6,*)' At end of loop over NOPN, the number of'
      WRITE(6,*)' configurations generated does not match'
      WRITE(6,*)' that which was allocated.'
      WRITE(6,*)' INFO table in ICNFTAB says:'
      WRITE(6,*)
      WRITE(6,*)'  NOPN ISYM       Nr of conf Start point'//
     &             '  Words/config'
      DO NOPN=MINOP,MAXOP
       NCLS=(NEL-NOPN)/2
       NOCC=NCLS+NOPN
       DO ISYM=1,NSYM
        NCNF=ICNFTAB(KINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))
        KCNFSTA=ICNFTAB(KINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
        LENCNF=ICNFTAB(KINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))
        WRITE(6,'(1X,2I4,5X,3I12)') NOPN,ISYM,NCNF,KCNFSTA,LENCNF
       END DO
      END DO
      CALL ABEND()
      END
