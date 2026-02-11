!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE SPHPUR(CMO)
      use define_af, only: iTabMx, AngTp
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: BName, IXSYM
      use general_data, only: NSYM,NBAS,NORB
      use Molcas, only: LenIn

      IMPLICIT None
      Real*8 CMO(*)

      CHARACTER(LEN=1) LCHAR
      Real*8 WGTLQN(0:9)
      LOGICAL IFTEST
      Integer, Allocatable:: LQN(:)
      Integer I, IB, IBAS, IBASES, ICMOES, IO, IORB, IORBES, ISSLAB,    &
     &        ISYM, ITP, L, LCOUNT, LEXIST, LMX, MNL, MXL, NB, NBTOT,   &
     &        NO, NONZ
      REAL*8 WGT, WMX

! Set IFTEST=.true. to get supsym input generated in the output
! for further use, or for testing.
      IFTEST=.false.
      IF(IFTEST) WRITE(6,*)'SUPSYM'

! Set up array with angular quant num for each basis function:
      NBTOT=0
      DO ISYM=1,NSYM
       NBTOT=NBTOT+NBAS(ISYM)
      END DO
      CALL mma_allocate(LQN,NBTOT,Label='LQN')
      DO IBAS=1,NBTOT
       LCHAR=BName(IBAS)(LenIn+3:LenIn+3)
       L=-999999
       DO ITP=0,ITABMX
         IF(LCHAR.EQ.ANGTP(ITP)) L=ITP
       END DO
       LQN(IBAS)=L
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
         L=LQN(IBAS)
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
! We have now a provisional IXSYM array. How many different
! L values appear in it?
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
! There are NONZ different values, so we want NONZ-1 special supsym
! labels for this symmetry. Reuse IWORK(LLQN) for orbital numbers:
! This will be the supsym label:
       IF (IFTEST) WRITE(6,*) NONZ-1
       ISSLAB=0
       DO L=MNL,MXL
        LCOUNT=0
        DO IO=1,NO
         IORB=IORBES+IO
         IF(L.EQ.IXSYM(IORB))THEN
          LCOUNT=LCOUNT+1
          LQN(LCOUNT)=IO
         END IF
        END DO
        IF(LCOUNT.GT.0) THEN
! Replace provisional IXSYM value with correct label:
          DO IO=1,NO
           IORB=IORBES+IO
           IF(IXSYM(IORB).eq.L) IXSYM(IORB)=ISSLAB
          END DO
! Lowest L = label zero = do not specify in input:
          IF (IFTEST.and.(ISSLAB.GT.0)) THEN
           WRITE(6,'(1x,I3,16I5,(/,5X,16I5))') LCOUNT,                  &
     &         (LQN(i),i=1,LCOUNT)
          END IF
          ISSLAB=ISSLAB+1
        END IF
       END DO

       ICMOES=ICMOES+NO*NB
       IORBES=IORBES+NO
 100   CONTINUE
       IBASES=IBASES+NB
      END DO
      CALL mma_deallocate(LQN)
      END SUBROUTINE SPHPUR
