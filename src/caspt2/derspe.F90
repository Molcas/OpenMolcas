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

subroutine DERSPE(NLEV,NG3,DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
! SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
! OR CLOSED-SHELL SCF CASE.

use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use sguga, only: LEVEL
use caspt2_module, only: ISCF, NACTEL
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, byte, u6

implicit none
integer(kind=iwp), intent(in) :: NLEV, NG3
real(kind=wp), intent(in) :: DF1(NLEV,NLEV), DF2(NLEV,NLEV,NLEV,NLEV), DF3(NG3), G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV), G3(NG3)
integer(kind=byte), intent(inout) :: idxG3(6,NG3)
real(kind=wp), intent(inout) :: DEPSA(NLEV,NLEV)
integer(kind=iwp) :: ID, iG3, IND1, IND2, IND3, IT, IT1, IT2, IT3, iTask, IU, IU1, IU2, IU3, IV, LT, LU, LU1, LU2, LU3, LV, NLEV2, &
                     NLEV4, nTask
real(kind=wp) :: DESUM, OCC

!ESUM = sum(ETA(1:NLEV))
DESUM = Zero
! ISCF=1 for closed-shell, =2 for hispin
OCC = Two
if (ISCF == 2) OCC = One

!if ((NACTEL == 1) .or. (NACTEL == 2)) NG3 = 0
if (NACTEL /= 1) then
  if (NACTEL /= 2) then
    write(u6,*) 'I have not implemented for non-standard Psi0, when A and C subspaces contribute to the energy, in particular'
    write(u6,*) 'I cannot debug, because I do not know when it happens'
    !call abend()

    NLEV2 = NLEV**2
    NLEV4 = NLEV**4

    iG3 = 0
    nTask = NLEV4
    ! SVC20100908 initialize the series of tasks
    call Init_Tsk(ID,nTask)

    do
#     ifdef _MOLCAS_MPP_
      if ((NG3-iG3) < NLEV2) exit
#     endif
      if (.not. Rsv_Tsk(ID,iTask)) exit

      IND1 = mod(iTask-1,NLEV2)+1
      IND2 = ((iTask-IND1)/(NLEV2))+1
      if (IND2 > IND1) cycle

      IT1 = mod(IND1-1,NLEV)+1
      IU1 = (IND1-IT1)/NLEV+1
      LU1 = LEVEL(IU1)
      IT2 = mod(IND2-1,NLEV)+1
      IU2 = (IND2-IT2)/NLEV+1
      LU2 = LEVEL(IU2)

      do IT3=1,NLEV
        do IU3=1,NLEV
          IND3 = IT3+NLEV*(IU3-1)
          if (IND3 > IND2) cycle
          LU3 = LEVEL(IU3)
          !VAL = G1(IT1,IU1)*G1(IT2,IU2)*G1(IT3,IU3)

          ! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
          ! Add here the necessary Kronecker deltas times 2-body matrix
          ! elements and lower, so we get a true normal-ordered density matrix
          ! element.

          ! <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
          ! = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
          ! -D(T3,U2)*(G2(T1,U1,T2,U3)+D(T2,U1)*G1(T1,U3))
          ! -D(T2,U1)*G2(T1,U2,T3,U3)
          ! -D(T3,U1)*G2(T2,U2,T1,U3)

          !if (IT3 == IU2) then
          !  VAL = VAL-G2(IT1,IU1,IT2,IU3)
          !  if (IT2 == IU1) VAL = VAL-G1(IT1,IU3)
          !end if
          !if (IT2 == IU1) VAL = VAL-G2(IT1,IU2,IT3,IU3)
          !if (IT3 == IU1) VAL = VAL-G2(IT2,IU2,IT1,IU3)

          ! VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
          iG3 = iG3+1
          idxG3(1,iG3) = int(iT1,kind=byte)
          idxG3(2,iG3) = int(iU1,kind=byte)
          idxG3(3,iG3) = int(iT2,kind=byte)
          idxG3(4,iG3) = int(iU2,kind=byte)
          idxG3(5,iG3) = int(iT3,kind=byte)
          idxG3(6,iG3) = int(iU3,kind=byte)
          !G3(iG3) = VAL
          !F3(iG3) = (ESUM*OCC-ETA(LU1)-ETA(LU2)-ETA(LU3))*VAL
          DESUM = DESUM+OCC*G3(iG3)*DF3(iG3)
          DEPSA(LU1,LU1) = DEPSA(LU1,LU1)-OCC*G3(iG3)*DF3(iG3)
          DEPSA(LU2,LU2) = DEPSA(LU2,LU2)-OCC*G3(iG3)*DF3(iG3)
          DEPSA(LU3,LU3) = DEPSA(LU3,LU3)-OCC*G3(iG3)*DF3(iG3)
        end do
      end do

      !SVC: The master node now continues to only handle task scheduling,
      !     needed to achieve better load balancing. So it exits from the task
      !     list.  It has to do it here since each process gets at least one
      !     task.
    end do

    ! SVC2010: no more tasks, wait here for the others.
    call Free_Tsk(ID)
  end if
  do IT=1,NLEV
    LT = LEVEL(IT)
    do IU=1,NLEV
      LU = LEVEL(IU)
      !G2(IT,IT,IU,IU) = G1(IT,IT)*G1(IU,IU)
      !if (IU == IT) then
      ! G2(IT,IT,IU,IU) = G2(IT,IT,IU,IU)-G1(IT,IU)
      !else
      ! G2(IT,IU,IU,IT) = -G1(IT,IT)
      !end if
      !F2(IT,IT,IU,IU) = (ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IT,IU,IU)
      !F2(IT,IU,IU,IT) = (ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IU,IU,IT)
      DESUM = DESUM+OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IU)
      DESUM = DESUM+OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IT)
      do IV=1,NLEV
        LV = LEVEL(IV)
        DEPSA(LT,LV) = DEPSA(LT,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IV,IU,IU)
        DEPSA(LU,LV) = DEPSA(LU,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IV)
        DEPSA(LT,LV) = DEPSA(LT,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IV)
        DEPSA(LU,LV) = DEPSA(LU,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IV,IT)
      end do
    end do
  end do
end if

do IT=1,NLEV
  !G1(IT,IT) = OCC
  LT = LEVEL(IT)
  !F1(IT,IT) = (ESUM*OCC-ETA(LT))*G1(IT,IT)
  DESUM = DESUM+OCC*G1(IT,IT)*DF1(IT,IT)
  do IU=1,NLEV
    LU = LEVEL(IU)
    DEPSA(LT,LU) = DEPSA(LT,LU)-OCC*G1(IT,IT)*DF1(IT,IU)
  end do
end do

end subroutine DERSPE
