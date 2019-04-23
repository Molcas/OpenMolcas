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
      SUBROUTINE SPHPUR_m(CMO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMO(*)
      CHARACTER*1 LCHAR
      DIMENSION WGTLQN(0:9)
      LOGICAL IFTEST
* Define mxsym etc.
#include "rasdim.fh"
* general.fh defines NSYM,NBAS,NORB:
#include "general.fh"
* rasscf.fh defines NAME:
#include "rasscf.fh"
* angtp.fh defines ITABMX,ANGTP
#include "angtp.fh"

#include "WrkSpc.fh"

* Set IFTEST=.true. to get supsym input generated in the output
* for further use, or for testing.
      IFTEST=.false.
      IF(IFTEST) WRITE(6,*)'SUPSYM'

* Set up array with angular quant num for each basis function:
      NBTOT=0
      DO ISYM=1,NSYM
       NBTOT=NBTOT+NBAS(ISYM)
      END DO
      CALL GETMEM('LQN','ALLO','INTE',LLQN,NBTOT)
      DO IBAS=1,NBTOT
       LCHAR=NAME(IBAS)(LENIN3:LENIN3)
       L=-999999
       DO ITP=0,ITABMX
         IF(LCHAR.EQ.ANGTP(ITP)) L=ITP
       END DO
       IWORK(LLQN-1+IBAS)=L
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
         WGTLQN(L)=0.0D0
        END DO
        DO IB=1,NB
         IBAS=IBASES+IB
         L=IWORK(LLQN-1+IBAS)
         WGT=CMO(ICMOES+IB+NB*(IO-1))**2
         WGTLQN(L)=WGTLQN(L)+WGT
        END DO
        LMX=0
        WMX=WGTLQN(0)
        DO L=0,9
         IF(WGTLQN(L).GT.WMX) THEN
           LMX=L
           WMX=WGTLQN(L)
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
* labels for this symmetry. Reuse IWORK(LLQN) for orbital numbers:
* This will be the supsym label:
       IF (IFTEST) WRITE(6,*) NONZ-1
       ISSLAB=0
       DO L=MNL,MXL
        LCOUNT=0
        DO IO=1,NO
         IORB=IORBES+IO
         IF(L.EQ.IXSYM(IORB))THEN
          LCOUNT=LCOUNT+1
          IWORK(LLQN-1+LCOUNT)=IO
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
     &         (IWORK(LLQN+i),i=0,LCOUNT-1)
          END IF
          ISSLAB=ISSLAB+1
        END IF
       END DO

       ICMOES=ICMOES+NO*NB
       IORBES=IORBES+NO
 100   CONTINUE
       IBASES=IBASES+NB
      END DO
      CALL GETMEM('LQN','FREE','INTE',LLQN,NBTOT)
      RETURN
      END
