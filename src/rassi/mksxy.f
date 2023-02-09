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
* Copyright (C) 1987, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE MKSXY(CMO1,CMO2,SXY)
      use OneDat, only: sNoNuc, sNoOri
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SXY(NSXY),CMO1(NCMO),CMO2(NCMO)
      character(len=8) :: LABEL
C  PURPOSE: FORM THE OVERLAP MATRIX SXY FOR ORBITAL BASES CMO1, CMO2.
C  CODED 1987-02-18, P-AA M.
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"
C  CALCULATE SIZE AND ALLOCATE A FIELD SZZ FOR OVERLAP MATRIX
C  IN COMMON BASIS SET (TRIANGULAR), SSQ TEMPORARY STORAGE
C  FOR EACH OF ITS SYMMETRY BLOCKS (SQUARE), AND PROD FOR
C  INTERMEDIATE MATRIX PRODUCTS.
      NSZZ=0
      NSSQ=0
      NPROD=0
      DO ISY=1,NSYM
        NO=NOSH(ISY)
        NB=NBASF(ISY)
        NSZZ=NSZZ+(NB*(NB+1))/2
        NSSQ=MAX(NSSQ,NB**2)
        NPROD=MAX(NPROD,NO*NB)
      end do
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NSZZ)
      CALL GETMEM('SSQ   ','ALLO','REAL',LSSQ,NSSQ)
      CALL GETMEM('PROD  ','ALLO','REAL',LPROD,NPROD)
C  READ OVERLAP MATRIX SZZ:
      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE MKSXY ***'
        WRITE(6,*)'     OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF
C  LOOP OVER SYMMETRIES:
      LSZZ1=LSZZ
      ISXY=1
      ICMO=1
      DO ISY=1,NSYM
        NB=NBASF(ISY)
        IF(NB.EQ.0) cycle
        NO=NOSH(ISY)
        if (NO /= 0) then
        CALL SQUARE(WORK(LSZZ1),WORK(LSSQ),1,NB,NB)
C  PROD:=SSQ*CMO2
        CALL DGEMM_('N','N',NB,NO,NB,1.0D0,WORK(LSSQ),NB,CMO2(ICMO),NB,
     &             0.0D0,WORK(LPROD),NB)
C  SXY:=(CMO1(TRANSP))*PROD
        CALL DGEMM_('T','N',NO,NO,NB,1.0D0,CMO1(ICMO),NB,WORK(LPROD),NB,
     &             0.0D0,SXY(ISXY),NO)
        ISXY=ISXY+NO**2
        ICMO=ICMO+NO*NB
        end if
        LSZZ1=LSZZ1+(NB*(NB+1))/2
      end do
      CALL GETMEM('      ','FREE','REAL',LSZZ,NSZZ)
      CALL GETMEM('      ','FREE','REAL',LSSQ,NSSQ)
      CALL GETMEM('      ','FREE','REAL',LPROD,NPROD)
      RETURN
      END
