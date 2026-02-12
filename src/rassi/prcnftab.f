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
      SUBROUTINE PRCNFTAB(CnfTab,MXPRT)
      use definitions, only: iwp, u6
      IMPLICIT none
      Integer(kind=iwp), intent(in):: CnfTab(*), MxPrt
      CHARACTER(LEN=144) TEXT

      Integer(kind=iwp) NPRT, NTAB, ITYPE, NEL, NORB, MINOP, MAXOP,
     &                  NSYM, LSYM, NGAS, IFORM, ICNF, IERR, IGAS,
     &                  ISUM, ISYM, KCNFSTA, KSTA, LENCNF, LENGTH,
     &                  LGASLIM, LGASORB, LIM1, LIM2, LINFO, NCLS,
     &                  NCNF, NOPN
C Sanity test:
      NPRT=MIN(MXPRT,10000)
C header:
      NTAB =CnfTab( 1)
      ITYPE=CnfTab( 2)
      NEL  =CnfTab( 3)
      NORB =CnfTab( 4)
      MINOP=CnfTab( 5)
      MAXOP=CnfTab( 6)
      NSYM =CnfTab( 7)
      LSYM =CnfTab( 8)
      NGAS =CnfTab( 9)
      IFORM=CnfTab(10)
      WRITE(u6,*)'---------------------------------------------------'
      WRITE(u6,*)'       Configuration Table Printout'
      WRITE(u6,*)' Table header:'
      WRITE(u6,'(1x,a,i16)')' Table size             NTAB=',NTAB
      WRITE(u6,'(1x,a,i16)')' Nr of electrons         NEL=',NEL
      WRITE(u6,'(1x,a,i16)')' Nr of orbitals         NORB=',NORB
      WRITE(u6,'(1x,a,i16)')' Min nr of open shells MINOP=',MINOP
      WRITE(u6,'(1x,a,i16)')' Max nr of open shells MAXOP=',MAXOP
      WRITE(u6,'(1x,a,i16)')' Point group order      NSYM=',NSYM
      WRITE(u6,'(1x,a,i16)')' Selected symmetry      LSYM=',LSYM
      WRITE(u6,'(1x,a,i16)')' Nr of GAS restrictions NGAS=',NGAS
      WRITE(u6,'(1x,a,i16)')' Configuration format  IFORM=',IFORM
      IF(ITYPE.NE.37) THEN
        WRITE(u6,*)' PRCNFTAB error: This is not a configuration table!'
        CALL ABEND()
      END IF
      IERR=0
      IF(NEL  .LT. 0) IERR=IERR+1
      IF(NEL  .GT.99) IERR=IERR+1
      IF(NORB .LT. 0) IERR=IERR+1
      IF(NORB .GT.199) IERR=IERR+1
      IF(MINOP.LT. 0) IERR=IERR+1
      IF(MAXOP.LT.MINOP) IERR=IERR+1
      IF(NSYM .LT. 1) IERR=IERR+1
      IF(NSYM .GT. 8) IERR=IERR+1
      IF(LSYM .LT. 0) IERR=IERR+1
      IF(LSYM .GT.NSYM) IERR=IERR+1
      IF(NGAS .LT. 0) IERR=IERR+1
      IF(NGAS .GT.19) IERR=IERR+1
      IF(IFORM.LT.1) IERR=IERR+1
      IF(IFORM.GT.4) IERR=IERR+1
      IF(IERR.GT.0) THEN
        WRITE(6,*)' PRCNFTAB error: Those values are unacceptable!'
        CALL ABEND()
      END IF
      LGASORB=11
      LGASLIM=LGASORB+(NSYM+1)*(NGAS+1)
      IF(NGAS.GT.0) THEN
        WRITE(u6,*)
        WRITE(u6,*)' Orbitals by symmetry in GAS partitions:'
        WRITE(u6,'(A,8I5)')' Space       Tot     ',(ISYM,ISYM=1,NSYM)
        WRITE(u6,'(1X,A,5X,I5,5X,8I5)')'Active',
     &          CnfTab(LGASORB),(CnfTab(LGASORB+ISYM),ISYM=1,NSYM)
        DO IGAS=1,NGAS
         WRITE(u6,'(1X,A,I2,5X,I5,5X,8I5)')'GAS',IGAS,
     &          CnfTab(LGASORB+(NSYM+1)*IGAS),
     &         (CnfTab(LGASORB+ISYM+(NSYM+1)*IGAS),ISYM=1,NSYM)
        END DO
        WRITE(u6,*)' GAS orbital partitions, min and max population:'
        DO IGAS=1,NGAS
         ISUM=CnfTab(LGASORB+(NSYM+1)*IGAS)
         LIM1=CnfTab(LGASLIM  +2*(IGAS-1))
         LIM2=CnfTab(LGASLIM+1+2*(IGAS-1))
         WRITE(u6,'(1X,A,I2,5X,3I5)')'GAS',IGAS,ISUM,LIM1,LIM2
        END DO
      END IF
      LINFO=LGASLIM+2*NGAS
      WRITE(u6,*)
      WRITE(u6,*)' INFO table starts at LINFO=',LINFO
      DO NOPN=MINOP,MAXOP
       NCLS=(NEL-NOPN)/2
       DO ISYM=1,NSYM
        NCNF   =CnfTab(LINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))
        KCNFSTA=CnfTab(LINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
        LENCNF =CnfTab(LINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))
        IF(NCNF.NE.0) THEN
         WRITE(u6,*)
         WRITE(u6,*)'  NOPN ISYM       Nr of conf Start point'//
     &             '  Words/config'
         WRITE(u6,'(1X,2I4,5X,3I12)') NOPN,ISYM,NCNF,KCNFSTA,LENCNF
         IF(NCNF.GT.NPRT) THEN
           WRITE(u6,'(1X,A,I5)')
     &             ' The first NPRT configurations. NPRT=',NPRT
         ELSE
           WRITE(u6,*)' Configurations:'
         END IF
         KSTA=KCNFSTA
         DO ICNF=1,MIN(NCNF,NPRT)
           CALL CNF2TXT(IFORM,NORB,NCLS,NOPN,CnfTab(KSTA),LENGTH,TEXT)
           KSTA=KSTA+LENCNF
           IF(LENGTH.LE.72) THEN
            WRITE(u6,'(8X,A)')TEXT(1:LENGTH)
           ELSE
            WRITE(u6,'(8X,A)')TEXT(1:72)
            WRITE(u6,'(8X,A)')TEXT(73:LENGTH)
           END IF
         END DO
         IF(NCNF.GT.NPRT) THEN
           WRITE(u6,*)' ( ...and more. This list was truncated.)'
         END IF
        END IF
       END DO
      END DO
      END SUBROUTINE PRCNFTAB
