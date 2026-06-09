!***********************************************************************
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
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine SPECIAL(G1,G2,G3,F1,F2,F3,idxG3,nLev,mG3)
! SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
! OR CLOSED-SHELL SCF CASE.

use caspt2_global, only: SGS
use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use caspt2_module, only: ETA, iSCF, NG3
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: nLev, mG3
real(kind=wp), intent(out) :: G1(nLev,nLev), G2(nLev,nLev,nLev,nLev), G3(mG3), F1(nLev,nLev), F2(nLev,nLev,nLev,nLev), F3(mG3)
integer(kind=byte), intent(out) :: idxG3(6,mG3)
integer(kind=iwp) :: ID, IG3, IND1, IND2, IND3, IT, IT1, IT2, IT3, ITASK, IU, IU1, IU2, IU3, LT, LU, LU1, LU2, LU3, NLEV2, NLEV4, &
                     NTASK
real(kind=wp) :: ESUM, Occ, Val

G1(:,:) = Zero
G2(:,:,:,:) = Zero
G3(:) = Zero
F1(:,:) = Zero
F2(:,:,:,:) = Zero
F3(:) = Zero

ESUM = sum(ETA(1:NLEV))
! ISCF=1 for closed-shell, =2 for hispin
OCC = Two
if (ISCF == 2) OCC = One
do IT=1,nLev
  G1(IT,IT) = OCC
  LT = SGS%LEVEL(IT)
  F1(IT,IT) = (ESUM*OCC-ETA(LT))*G1(IT,IT)
end do

! This code can be activated once the subsequent code can handle it.
! As it is now the code process G3 regardless of the setting here.
! If nG3=0 this results in NAN processing. Until this is explicitly
! taken case of this code should remain commented out.
!if (NACTEL == 1) then
!  NG3 = 0
!  return
!end if

do IT=1,nLev
  LT = SGS%LEVEL(IT)
  do IU=1,nLev
    LU = SGS%LEVEL(IU)
    G2(IT,IT,IU,IU) = G1(IT,IT)*G1(IU,IU)
    if (IU == IT) then
      G2(IT,IT,IU,IU) = G2(IT,IT,IU,IU)-G1(IT,IU)
    else
      G2(IT,IU,IU,IT) = -G1(IT,IT)
    end if
    F2(IT,IT,IU,IU) = (ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IT,IU,IU)
    F2(IT,IU,IU,IT) = (ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IU,IU,IT)
  end do
end do

!if (NACTEL == 2) then
!  NG3 = 0
!  return
!end if

NLEV2 = NLEV**2
NLEV4 = NLEV**4

iG3 = 0
nTask = NLEV4
! SVC20100908 initialize the series of tasks
call Init_Tsk(ID,nTask)

Outer: do
# ifdef _MOLCAS_MPP_
  if ((NG3-iG3) < NLEV2) exit Outer
# endif

  Inner: do
    if (.not. Rsv_Tsk(ID,iTask)) exit Outer
    IND1 = mod(iTask-1,NLEV2)+1
    IND2 = ((iTask-IND1)/(NLEV2))+1
    if (IND2 <= IND1) exit Inner
  end do Inner

  IT1 = mod(IND1-1,nLev)+1
  IU1 = (IND1-IT1)/nLev+1
  LU1 = SGS%LEVEL(IU1)
  IT2 = mod(IND2-1,nLev)+1
  IU2 = (IND2-IT2)/nLev+1
  LU2 = SGS%LEVEL(IU2)

  do IT3=1,NLEV
    do IU3=1,NLEV
      IND3 = IT3+nLev*(IU3-1)
      if (IND3 > IND2) cycle
      LU3 = SGS%LEVEL(IU3)
      VAL = G1(IT1,IU1)*G1(IT2,IU2)*G1(IT3,IU3)

      ! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
      ! Add here the necessary Kronecker deltas times 2-body matrix
      ! elements and lower, so we get a true normal-ordered density matrix
      ! element.

      ! <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
      ! = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
      ! -D(T3,U2)*(G2(T1,U1,T2,U3)+D(T2,U1)*G1(T1,U3))
      ! -D(T2,U1)*G2(T1,U2,T3,U3)
      ! -D(T3,U1)*G2(T2,U2,T1,U3)

      if (IT3 == IU2) then
        VAL = VAL-G2(IT1,IU1,IT2,IU3)
        if (IT2 == IU1) VAL = VAL-G1(IT1,IU3)
      end if
      if (IT2 == IU1) VAL = VAL-G2(IT1,IU2,IT3,IU3)
      if (IT3 == IU1) VAL = VAL-G2(IT2,IU2,IT1,IU3)

      ! VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
      iG3 = iG3+1
      idxG3(1,iG3) = int(iT1,kind=byte)
      idxG3(2,iG3) = int(iU1,kind=byte)
      idxG3(3,iG3) = int(iT2,kind=byte)
      idxG3(4,iG3) = int(iU2,kind=byte)
      idxG3(5,iG3) = int(iT3,kind=byte)
      idxG3(6,iG3) = int(iU3,kind=byte)
      G3(iG3) = VAL
      F3(iG3) = (ESUM*OCC-ETA(LU1)-ETA(LU2)-ETA(LU3))*VAL

    end do
  end do

  !SVC: The master node now continues to only handle task scheduling,
  !     needed to achieve better load balancing. So it exits from the task
  !     list.  It has to do it here since each process gets at least one
  !     task.

end do Outer

! SVC2010: no more tasks, wait here for the others.
call Free_Tsk(ID)

NG3 = iG3

end subroutine SPECIAL
