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
      SUBROUTINE PrSSTab(SSTAB)
      use cntrl, only: MORSBITS
      IMPLICIT NONE
      INTEGER SSTAB(*)

      INTEGER LPOS
      INTEGER NSSTP,ISSTP,KSSTP,NSBS,NPOP,ISYM,MS2,ISPART
      INTEGER KSSTANN,KSSTCRE,KSBSMRS,KMRSSBS,KSBSANN,KSBSCRE
      INTEGER I,ISBS,NSBSTOT,NMORS,IMRS,NASPRT
      INTEGER NRSBST,ISBSSTA,ISBSEND,IMRSSTA,IMRSEND
      WRITE(6,*)
      WRITE(6,*)'============================================='
      WRITE(6,*)' Substring table printout.'
      WRITE(6,*)' (A) Header and TOC:'
      WRITE(6,'(a,i16)')'            Table size:',SSTAB(1)
      WRITE(6,'(a,i16)')'       Table type code:',SSTAB(2)
      WRITE(6,'(a,i16)')' Orbital Table        :',SSTAB(3)
      WRITE(6,'(a,i16)')' Nr of symmetry labels:',SSTAB(4)
      WRITE(6,'(a,i16)')' Nr of active subpart :',SSTAB(5)
      WRITE(6,'(a,i16)')' Nr of bits/morsel    :',SSTAB(6)
      WRITE(6,'(a,i16)')' Nr of Substring Types:',SSTAB(7)
      WRITE(6,'(a,i16)')' Nr of Substrings     :',SSTAB(8)
      WRITE(6,'(a,i16)')' Substr Type Ann Table:',SSTAB(9)
      WRITE(6,'(a,i16)')' Substr Type Cre Table:',SSTAB(10)
      WRITE(6,'(a,i16)')' Substr/Morsel Table  :',SSTAB(11)
      WRITE(6,'(a,i16)')' Morsel/Substr Table  :',SSTAB(12)
      WRITE(6,'(a,i16)')' Substring Annih Table:',SSTAB(13)
      WRITE(6,'(a,i16)')' Substring Creat Table:',SSTAB(14)
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (B) Substring Types:'
      WRITE(6,*)' NSBS  = Nr of substrings of each type.'
      WRITE(6,*)' NPOP  = Electron population of each type.'
      WRITE(6,*)' ISYM  = Combined symmetry of each type.'
      WRITE(6,*)' MS2   = Twice spin proj of each type.'
      WRITE(6,*)' ISPART= Subpartition of each type.'
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)'           NSBS   NPOP   ISYM    MS2  ISBPRT'
      NSSTP=SSTAB(7)
      NSBSTOT =SSTAB(8)
      KSSTP  =15
      DO ISSTP=1,NSSTP
        NSBS   = SSTAB(KSSTP+0+5*(ISSTP-1))
        NPOP   = SSTAB(KSSTP+1+5*(ISSTP-1))
        ISYM   = SSTAB(KSSTP+2+5*(ISSTP-1))
        MS2    = SSTAB(KSSTP+3+5*(ISSTP-1))
        ISPART = SSTAB(KSSTP+4+5*(ISSTP-1))
        WRITE(6,'(1x,6I7)') ISSTP,NSBS,NPOP,ISYM,MS2,ISPART
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (C1) Morsels to Substrings:'
      NASPRT=SSTAB(5)
      KSBSMRS=SSTAB(11)
      KMRSSBS=SSTAB(12)
      NMORS=2**MORSBITS
      DO IMRSSTA=0,NMORS-1,15
        IMRSEND=MIN(NMORS-1,IMRSSTA+14)
        WRITE(6,*)
        WRITE(6,'(1x,a,3x,15I5)')'MRS:',(IMRS,IMRS=IMRSSTA,IMRSEND)
        DO ISPART=1,NASPRT
        WRITE(6,'(1x,a,i3,15I5)')'SBS:',ISPART,(SSTAB(KMRSSBS
     &        +2*(IMRS+NMORS*(ISPART-1))),IMRS=IMRSSTA,IMRSEND)
        END DO
      END DO

      NRSBST=SSTAB(8)
      WRITE(6,*)' (C2) Substrings to Morsels:'
      DO ISBSSTA=1,NRSBST,15
        ISBSEND=MIN(NRSBST,ISBSSTA+14)
        WRITE(6,*)
        WRITE(6,'(1x,a,15I5)')'SBS:',(ISBS,ISBS=ISBSSTA,ISBSEND)
        WRITE(6,'(1x,a,15I5)')'MRS:',(SSTAB(KSBSMRS
     &                +2*(ISBS-1)),ISBS=ISBSSTA,ISBSEND)
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (D1) Substring Type Annihilation Table:'
      KSSTANN=SSTAB(9)
      KSSTCRE=SSTAB(10)
      DO ISSTP=1,NSSTP
       LPOS=KSSTANN+MORSBITS*(ISSTP-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISSTP,(SSTAB(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)' (D2) Substring Type Creation Table:'
      DO ISSTP=1,NSSTP
       LPOS=KSSTCRE+MORSBITS*(ISSTP-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISSTP,(SSTAB(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (E1) Substring Annihilation Table:'
      KSBSANN=SSTAB(13)
      KSBSCRE=SSTAB(14)
      DO ISBS=1,NSBSTOT
       LPOS=KSBSANN+MORSBITS*(ISBS-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISBS,(SSTAB(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' (E2) Substring Creation Table:'
      DO ISBS=1,NSBSTOT
       LPOS=KSBSCRE+MORSBITS*(ISBS-1)
       WRITE(6,'(1x,I8,5X,8I8)')
     &                      ISBS,(SSTAB(LPOS-1+I),I=1,MORSBITS)
      END DO
      WRITE(6,*)'============================================='

      END SUBROUTINE PrSSTab
