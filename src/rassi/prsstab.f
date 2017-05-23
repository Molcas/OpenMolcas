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
      SUBROUTINE PrSSTab(LSSTAB)
      IMPLICIT NONE
      INTEGER LSSTAB,LPOS
#include "WrkSpc.fh"
      INTEGER NSSTP,ISSTP,KSSTP,NSBS,NPOP,ISYM,MS2,ISPART
      INTEGER KSSTANN,KSSTCRE,KSBSMRS,KMRSSBS,KSBSANN,KSBSCRE
      INTEGER I,ISBS,NSBSTOT,NMORS,IMRS,NASPRT
      INTEGER NRSBST,ISBSSTA,ISBSEND,IMRSSTA,IMRSEND
#include "Morsel.fh"
      WRITE(6,*)
      WRITE(6,*)'============================================='
      WRITE(6,*)' Substring table printout.'
      WRITE(6,*)' (A) Header and TOC:'
      WRITE(6,'(a,i16)')'     Workspace pointer:',LSSTAB
      WRITE(6,'(a,i16)')'            Table size:',IWORK(LSSTAB)
      WRITE(6,'(a,i16)')'       Table type code:',IWORK(LSSTAB+1)
      WRITE(6,'(a,i16)')' Orbital Table        :',IWORK(LSSTAB+2)
      WRITE(6,'(a,i16)')' Nr of symmetry labels:',IWORK(LSSTAB+3)
      WRITE(6,'(a,i16)')' Nr of active subpart :',IWORK(LSSTAB+4)
      WRITE(6,'(a,i16)')' Nr of bits/morsel    :',IWORK(LSSTAB+5)
      WRITE(6,'(a,i16)')' Nr of Substring Types:',IWORK(LSSTAB+6)
      WRITE(6,'(a,i16)')' Nr of Substrings     :',IWORK(LSSTAB+7)
      WRITE(6,'(a,i16)')' Substr Type Ann Table:',IWORK(LSSTAB+8)
      WRITE(6,'(a,i16)')' Substr Type Cre Table:',IWORK(LSSTAB+9)
      WRITE(6,'(a,i16)')' Substr/Morsel Table  :',IWORK(LSSTAB+10)
      WRITE(6,'(a,i16)')' Morsel/Substr Table  :',IWORK(LSSTAB+11)
      WRITE(6,'(a,i16)')' Substring Annih Table:',IWORK(LSSTAB+12)
      WRITE(6,'(a,i16)')' Substring Creat Table:',IWORK(LSSTAB+13)
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (B) Substring Types:'
      WRITE(6,*)' NSBS  = Nr of substrings of each type.'
      WRITE(6,*)' NPOP  = Electron population of each type.'
      WRITE(6,*)' ISYM  = Combined symmetry of each type.'
      WRITE(6,*)' MS2   = Twice spin proj of each type.'
      WRITE(6,*)' ISPART= Subpartition of each type.'
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)'           NSBS   NPOP   ISYM    MS2  ISBPRT'
      NSSTP=IWORK(LSSTAB+6)
      NSBSTOT =IWORK(LSSTAB+7)
      KSSTP  =15
      DO ISSTP=1,NSSTP
        NSBS   = IWORK(LSSTAB-1+KSSTP+0+5*(ISSTP-1))
        NPOP   = IWORK(LSSTAB-1+KSSTP+1+5*(ISSTP-1))
        ISYM   = IWORK(LSSTAB-1+KSSTP+2+5*(ISSTP-1))
        MS2    = IWORK(LSSTAB-1+KSSTP+3+5*(ISSTP-1))
        ISPART = IWORK(LSSTAB-1+KSSTP+4+5*(ISSTP-1))
        WRITE(6,'(1x,6I7)') ISSTP,NSBS,NPOP,ISYM,MS2,ISPART
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (C1) Morsels to Substrings:'
      NASPRT=IWORK(LSSTAB+4)
      KSBSMRS=IWORK(LSSTAB+10)
      KMRSSBS=IWORK(LSSTAB+11)
      NMORS=2**MORSBITS
      DO IMRSSTA=0,NMORS-1,15
        IMRSEND=MIN(NMORS-1,IMRSSTA+14)
        WRITE(6,*)
        WRITE(6,'(1x,a,3x,15I5)')'MRS:',(IMRS,IMRS=IMRSSTA,IMRSEND)
        DO ISPART=1,NASPRT
        WRITE(6,'(1x,a,i3,15I5)')'SBS:',ISPART,(IWORK(LSSTAB-1+KMRSSBS
     &        +2*(IMRS+NMORS*(ISPART-1))),IMRS=IMRSSTA,IMRSEND)
        END DO
      END DO

      NRSBST=IWORK(LSSTAB+7)
      WRITE(6,*)' (C2) Substrings to Morsels:'
      DO ISBSSTA=1,NRSBST,15
        ISBSEND=MIN(NRSBST,ISBSSTA+14)
        WRITE(6,*)
        WRITE(6,'(1x,a,15I5)')'SBS:',(ISBS,ISBS=ISBSSTA,ISBSEND)
        WRITE(6,'(1x,a,15I5)')'MRS:',(IWORK(LSSTAB-1+KSBSMRS
     &                +2*(ISBS-1)),ISBS=ISBSSTA,ISBSEND)
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (D1) Substring Type Annihilation Table:'
      KSSTANN=IWORK(LSSTAB+8)
      KSSTCRE=IWORK(LSSTAB+9)
      DO ISSTP=1,NSSTP
       LPOS=LSSTAB-1+KSSTANN+MORSBITS*(ISSTP-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISSTP,(IWORK(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)' (D2) Substring Type Creation Table:'
      DO ISSTP=1,NSSTP
       LPOS=LSSTAB-1+KSSTCRE+MORSBITS*(ISSTP-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISSTP,(IWORK(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (E1) Substring Annihilation Table:'
      KSBSANN=IWORK(LSSTAB+12)
      KSBSCRE=IWORK(LSSTAB+13)
      DO ISBS=1,NSBSTOT
       LPOS=LSSTAB-1+KSBSANN+MORSBITS*(ISBS-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISBS,(IWORK(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (E2) Substring Creation Table:'
      DO ISBS=1,NSBSTOT
       LPOS=LSSTAB-1+KSBSCRE+MORSBITS*(ISBS-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISBS,(IWORK(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)'============================================='
      RETURN
      END
