!u**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Per Ake Malmqvist                                *
!               2021, Yoshio Nishimoto                                 *
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      SUBROUTINE DERSPE(NLEV,NG3,DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
      use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
      use sguga, only: LEVEL
      use caspt2_module, only: NACTEL, ISCF
      use caspt2_module, only: ETA
      use Constants, only: Zero, One, Two
      use definitions, only: wp, iwp, byte, u6

      implicit none

      integer(kind=iwp), intent(in) :: NLEV, NG3
      real(kind=wp), intent(in) :: DF1(NLEV,NLEV),                      &
     &  DF2(NLEV,NLEV,NLEV,NLEV), DF3(NG3), G1(NLEV,NLEV),              &
     &  G2(NLEV,NLEV,NLEV,NLEV), G3(NG3)
      integer(kind=byte), intent(inout) :: idxG3(6,NG3)
      real(kind=wp), intent(inout) :: DEPSA(NLEV,NLEV)

! SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
! OR CLOSED-SHELL SCF CASE.

      integer(kind=iwp) :: I, NLEV2, NLEV4, iG3, nTask, ID, iTask,      &
     &  IND1, IND2, IT1, IU1, LU1, IT2, IU2, LU2, IT3, IU3, IND3, LU3,  &
     &  LT, IT, IU, LU, IV, LV
      real(kind=wp) :: ESUM, DESUM, OCC

      ESUM=Zero
      DESUM=Zero
      DO I=1,NLEV
        ESUM=ESUM+ETA(I)
      END DO
! ISCF=1 for closed-shell, =2 for hispin
      OCC=Two
      IF(ISCF == 2) OCC=One

!     if (NACTEL == 1 .or. NACTEL == 2) NG3 = 0
      if (NACTEL /= 1) then
        if (NACTEL /= 2) then
          write(u6,*) 'I have not implemented for non-standard Psi0, ', &
     &      'when A and C subspaces contribute to the energy, ',        &
     &      'in particular'
          write(u6,*) 'I cannot debug, ',                               &
     &      'because I do not know when it happens'
!         call abend()

          NLEV2=NLEV**2
          NLEV4=NLEV**4

          iG3=0
          nTask=NLEV4
! SVC20100908 initialize the series of tasks
          Call Init_Tsk(ID, nTask)

          do
#ifdef _MOLCAS_MPP_
            IF ((NG3-iG3) < NLEV2) exit
#endif
            IF (.NOT.Rsv_Tsk(ID,iTask)) exit

            IND1=MOD(iTask-1,NLEV2)+1
            IND2=((iTask-IND1)/(NLEV2))+1
            IF(IND2 > IND1) cycle

            IT1=MOD(IND1-1,NLEV)+1
            IU1=(IND1-IT1)/NLEV+1
            LU1=LEVEL(IU1)
            IT2=MOD(IND2-1,NLEV)+1
            IU2=(IND2-IT2)/NLEV+1
            LU2=LEVEL(IU2)

            DO IT3=1,NLEV
             DO IU3=1,NLEV
              IND3=IT3+NLEV*(IU3-1)
              IF(IND3 > IND2) cycle
              LU3=LEVEL(IU3)
!             VAL=G1(IT1,IU1)*G1(IT2,IU2)*G1(IT3,IU3)

! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
! Add here the necessary Kronecker deltas times 2-body matrix
! elements and lower, so we get a true normal-ordered density matrix
! element.

! <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
! = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
! -D(T3,U2)*(G2(T1,U1,T2,U3)+D(T2,U1)*G1(T1,U3))
! -D(T2,U1)*G2(T1,U2,T3,U3)
! -D(T3,U1)*G2(T2,U2,T1,U3)

!             IF(IT3 == IU2) THEN
!               VAL=VAL-G2(IT1,IU1,IT2,IU3)
!               IF(IT2 == IU1) THEN
!                 VAL=VAL-G1(IT1,IU3)
!               END IF
!             END IF
!             IF(IT2 == IU1) THEN
!               VAL=VAL-G2(IT1,IU2,IT3,IU3)
!             END IF
!             IF(IT3 == IU1) THEN
!               VAL=VAL-G2(IT2,IU2,IT1,IU3)
!             END IF

! VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
              iG3=iG3+1
              idxG3(1,iG3)=INT(iT1,kind=byte)
              idxG3(2,iG3)=INT(iU1,kind=byte)
              idxG3(3,iG3)=INT(iT2,kind=byte)
              idxG3(4,iG3)=INT(iU2,kind=byte)
              idxG3(5,iG3)=INT(iT3,kind=byte)
              idxG3(6,iG3)=INT(iU3,kind=byte)
!             G3(iG3)=VAL
!             F3(iG3)=(ESUM*OCC-ETA(LU1)-ETA(LU2)-ETA(LU3))*VAL
              DESUM = DESUM + OCC*G3(iG3)*DF3(iG3)
              DEPSA(LU1,LU1) = DEPSA(LU1,LU1) - OCC*G3(iG3)*DF3(iG3)
              DEPSA(LU2,LU2) = DEPSA(LU2,LU2) - OCC*G3(iG3)*DF3(iG3)
              DEPSA(LU3,LU3) = DEPSA(LU3,LU3) - OCC*G3(iG3)*DF3(iG3)
              END DO
            END DO

!SVC: The master node now continues to only handle task scheduling,
!     needed to achieve better load balancing. So it exits from the task
!     list.  It has to do it here since each process gets at least one
!     task.
          end do

! SVC2010: no more tasks, wait here for the others.
          CALL Free_Tsk(ID)
        end if
        DO IT=1,NLEV
         LT=LEVEL(IT)
         DO IU=1,NLEV
          LU=LEVEL(IU)
!         G2(IT,IT,IU,IU)=G1(IT,IT)*G1(IU,IU)
!         IF(IU == IT) THEN
!          G2(IT,IT,IU,IU)=G2(IT,IT,IU,IU)-G1(IT,IU)
!         ELSE
!          G2(IT,IU,IU,IT)=-G1(IT,IT)
!         END IF
!         F2(IT,IT,IU,IU)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IT,IU,IU)
!         F2(IT,IU,IU,IT)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IU,IU,IT)
          DESUM = DESUM + OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IU)
          DESUM = DESUM + OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IT)
          DO IV=1,NLEV
          LV=LEVEL(IV)
          DEPSA(LT,LV)=DEPSA(LT,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IV,IU,IU)
          DEPSA(LU,LV)=DEPSA(LU,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IV)
          DEPSA(LT,LV)=DEPSA(LT,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IV)
          DEPSA(LU,LV)=DEPSA(LU,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IV,IT)
          END DO
         END DO
        END DO
      end if

      DO IT=1,NLEV
!       G1(IT,IT)=OCC
        LT=LEVEL(IT)
!       F1(IT,IT)=(ESUM*OCC-ETA(LT))*G1(IT,IT)
        DESUM = DESUM + OCC*G1(IT,IT)*DF1(IT,IT)
        Do IU=1, NLEV
          LU=LEVEL(IU)
          DEPSA(LT,LU) = DEPSA(LT,LU) - OCC*G1(IT,IT)*DF1(IT,IU)
        End Do
      END DO

      END subroutine DERSPE
