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
      SUBROUTINE LINPUR(CMO)
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: BName, IXSYM
      use general_data, only: NSYM,NBAS,NORB
      use Molcas, only: LenIn

      IMPLICIT None
      Real*8 CMO(*)

      CHARACTER(LEN=2) LCHAR
      Real*8 WGTLMB(0:9)
      LOGICAL IFTEST
* Define mxsym etc.
      Integer, Allocatable:: LMB(:)
      Integer i, IB, IBAS, IBASES, ICMOES, IO, IORB, IORBES, ISSLAB,
     &        ISYM, L, LCOUNT, LEXIST, LMX, MNL, MXL, NB, NBTOT, NO,
     &        NONZ
      REAL*8 WGT, WMX

* Set IFTEST=.true. to get supsym input generated in the output
* for further use, or for testing.
      IFTEST=.false.

      IF(IFTEST) WRITE(6,*)'SUPSYM'

* Set up array with lambda for each basis function:
      NBTOT=0
      DO ISYM=1,NSYM
       NBTOT=NBTOT+NBAS(ISYM)
      END DO
      CALL mma_allocate(LMB,NBTOT,Label='LMB')
      DO IBAS=1,NBTOT
       LCHAR=BName(IBAS)(LenIn+4:LenIn+5)
       IF(LCHAR.EQ.'  ') THEN
         L=0
       ELSE IF(LCHAR.EQ.'x ') THEN
         L=1
       ELSE IF(LCHAR.EQ.'y ') THEN
         L=-1
       ELSE IF(LCHAR.EQ.'z ') THEN
         L=0
       ELSE
         READ(LCHAR,'(I2)') L
         IF (BName(IBAS)(LenIn+5:LenIn+5).EQ.'-') THEN
           L=-L
         END IF
       END IF
       LMB(IBAS)=ABS(L)
      END DO
      ICMOES=0
      IBASES=0
      IORBES=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       NO=NORB(ISYM)
       IF(NO.EQ.0) GOTO 100
       DO IO=1,NO
        IORB=IORBES+IO
        DO L=0,9
         WGTLMB(L)=0.0D0
        END DO
        DO IB=1,NB
         IBAS=IBASES+IB
         L=LMB(IBAS)
         WGT=CMO(ICMOES+IB+NB*(IO-1))**2
         WGTLMB(L)=WGTLMB(L)+WGT
        END DO
        LMX=0
        WMX=WGTLMB(0)
        DO L=0,9
         IF(WGTLMB(L).GT.WMX) THEN
           LMX=L
           WMX=WGTLMB(L)
         END IF
        END DO
        IXSYM(IORB)=LMX
       END DO
* We have now a provisional IXSYM array. How many different
* L values appear in it?
       MNL=9
       MXL=0
       NONZ=0
       DO L=0,9
        LEXIST=0
        DO IO=1,NO
         IORB=IORBES+IO
         IF(L.EQ.IXSYM(IORB))THEN
          LEXIST=1
          MNL=MIN(L,MNL)
          MXL=MAX(L,MXL)
          GOTO 19
         END IF
        END DO
  19    CONTINUE
        NONZ=NONZ+LEXIST
       END DO
* There are NONZ different values, so we want NONZ-1 special supsym
* labels for this symmetry. Reuse LMB for orbital numbers:
* This will be the supsym label:
       IF (IFTEST) WRITE(6,*) NONZ-1
       ISSLAB=0
       DO L=MNL,MXL
        LCOUNT=0
        DO IO=1,NO
         IORB=IORBES+IO
         IF(L.EQ.IXSYM(IORB))THEN
          LCOUNT=LCOUNT+1
          LMB(LCOUNT)=IO
         END IF
        END DO
        IF(LCOUNT.GT.0) THEN
* Replace provisional IXSYM value with correct label:
          DO IO=1,NO
           IORB=IORBES+IO
           IF(IXSYM(IORB).eq.L) IXSYM(IORB)=ISSLAB
          END DO
* Lowest L = label zero = do not specify in input:
          IF (IFTEST.and.(ISSLAB.GT.0)) THEN
           WRITE(6,'(1x,I3,16I5,(/,5X,16I5))') LCOUNT,
     &         (LMB(i),i=1,LCOUNT)
          END IF
          ISSLAB=ISSLAB+1
        END IF
       END DO

       ICMOES=ICMOES+NO*NB
       IORBES=IORBES+NO
 100   CONTINUE
       IBASES=IBASES+NB
      END DO
      CALL mma_deallocate(LMB)

      END SUBROUTINE LINPUR
