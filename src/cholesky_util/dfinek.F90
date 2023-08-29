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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine DfineK(K_Lap,R,InitR,Demand)
!-----------------------------------------------------------------------
! Function : Define optimal K value
!
! Demand : Select the accuracy which you want
!        = MILLI ... 10E-03 accuracy (K = 2-)
!        = MICRO ... 10E-06 accuracy (K = 4-)
!        = NANO  ... 10E-09 accuracy (K = 7-)
!        = PICO  ... 10E-12 accuracy (K =11-)
!-----------------------------------------------------------------------

use ReMez_mod, only: IW
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: K_Lap
real(kind=wp), intent(in) :: R
integer(kind=iwp), intent(in) :: InitR
character(len=8), intent(in) :: Demand
integer(kind=iwp) :: IdxR, IVal
real(kind=wp) :: ErrVal, R_Val
logical(kind=iwp) :: Trial
integer(kind=iwp), parameter :: IdxMin(20) = [1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,4,6,7,8], &
                                KList1(31) = [2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3], &
                                KList2(31) = [3,4,5,6,6,7,7,7,7,7,7,7,8,8,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10], &
                                KList3(31) = [5,7,8,9,10,10,11,11,11,11,12,12,13,14,14,14,15,15,15,15,15,16,17,17,17,18,18,18,18, &
                                              18,18], &
                                KList4(31) = [6,10,10,12,13,14,14,15,15,16,16,16,18,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20, &
                                              20,20,20]
real(kind=wp), parameter :: ErLst1(20) = [Zero,8.752e-3_wp,5.052e-3_wp,1.700e-3_wp,6.428e-4_wp,2.646e-4_wp,1.163e-4_wp, &
                                          5.392e-5_wp,2.611e-5_wp,1.312e-5_wp,6.807e-6_wp,3.630e-6_wp,1.984e-6_wp,1.108e-6_wp, &
                                          6.311e-7_wp,3.659e-7_wp,2.155e-7_wp,1.289e-7_wp,7.811e-8_wp,4.794e-8_wp], &
                            ErLst2(20) = [Zero,Zero,Zero,6.258e-6_wp,4.243e-6_wp,7.379e-6_wp,9.841e-6_wp,8.724e-6_wp,7.468e-6_wp, &
                                          9.296e-6_wp,3.844e-6_wp,1.582e-6_wp,6.481e-7_wp,2.648e-7_wp,1.079e-7_wp,4.388e-8_wp, &
                                          1.780e-8_wp,7.213e-9_wp,2.918e-9_wp,1.179e-9_wp], &
                            ErLst3(20) = [Zero,Zero,Zero,Zero,Zero,7.741e-9_wp,1.716e-9_wp,4.045e-9_wp,6.822e-9_wp,9.323e-9_wp, &
                                          3.357e-9_wp,4.301e-9_wp,7.587e-9_wp,7.555e-9_wp,7.134e-9_wp,8.484e-9_wp,7.213e-9_wp, &
                                          8.602e-9_wp,8.602e-9_wp,9.700e-9_wp], &
                            ErLst4(20) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,9.021e-12_wp,6.492e-13_wp,5.389e-12_wp, &
                                          5.735e-12_wp,9.405e-12_wp,6.723e-12_wp,5.298e-12_wp,1.050e-12_wp,3.128e-12_wp, &
                                          7.099e-12_wp,9.371e-12_wp], &
                            ErLst5(20) = [8.556e-2_wp,1.785e-2_wp,5.052e-3_wp,1.700e-3_wp,6.428e-4_wp,2.646e-4_wp,1.163e-4_wp, &
                                          5.392e-5_wp,2.611e-5_wp,1.312e-5_wp,6.807e-6_wp,3.630e-6_wp,1.984e-6_wp,1.108e-6_wp, &
                                          6.311e-7_wp,3.659e-7_wp,2.155e-7_wp,1.289e-7_wp,7.811e-8_wp,4.794e-8_wp], &
                            RList(30) = [2.0e0_wp,5.0e0_wp,1.0e1_wp,2.0e1_wp,3.0e1_wp,4.0e1_wp,5.0e1_wp,6.0e1_wp,7.0e1_wp, &
                                         8.0e1_wp,9.0e1_wp,1.0e2_wp,2.0e2_wp,3.0e2_wp,4.0e2_wp,5.0e2_wp,6.0e2_wp,7.0e2_wp, &
                                         8.0e2_wp,9.0e2_wp,1.0e3_wp,2.0e3_wp,3.0e3_wp,4.0e3_wp,5.0e3_wp,6.0e3_wp,7.0e3_wp, &
                                         8.0e3_wp,9.0e3_wp,1.0e4_wp]
character(len=*), parameter :: Micro = 'MICRO   ', Milli = 'MILLI   ', Nano = 'NANO    ', Pico = 'PICO    '

! ===== Check a larger R value =====

Trial = .false.

write(IW,'(A,A8,A)') 'Demanded accuracy is ',Demand,'.'

if (InitR == 31) then
  if (Demand == Milli) then
    K_Lap = 3
    ErrVal = ErLst5(K_Lap)
    write(IW,'(/A,E11.4E2)') ' This K guarantees the error less than ',ErrVal
  else if (Demand == Micro) then
    K_Lap = 11
    ErrVal = ErLst5(K_Lap)
    write(IW,'(/A,E11.4E2)') ' This K guarantees the error less than ',ErrVal
  else
    K_Lap = 20
    ErrVal = ErLst5(K_Lap)
    if (Demand == Nano) then
      if (R <= 3.0e4_wp) then
        K_Lap = 19
        ErrVal = ErLst3(14)
        write(IW,'(/A,E11.4E2)') ' This K guarantees the error less than ',ErrVal
        return
      else if (R <= 1.0e5_wp) then
        ErrVal = ErLst3(15)
        write(IW,'(/A,E11.4E2)') ' This K guarantees the error less than ',ErrVal
        return
      end if
    end if
    write(IW,'(/A)') '!!! Caution !!!'
    write(IW,'(A,E11.4E2,A)') 'In this R value, we can only guarantee',ErrVal,' accuracy.'
  end if
  return
end if

if (InitR <= 8) Trial = .true.

IVal = InitR+1

if (Demand == Milli) then
  K_Lap = KList1(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < Zero) IVal = IdxR
  end if
  ErrVal = ErLst1(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst1(K_Lap),' .'
else if (Demand == Micro) then
  K_Lap = KList2(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < Zero) IVal = IdxR
  end if
  ErrVal = ErLst2(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst2(K_Lap),' .'
else if (Demand == Nano) then
  K_Lap = KList3(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < Zero) IVal = IdxR
  end if
  ErrVal = ErLst3(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst3(K_Lap),' .'
else if (Demand == Pico) then
  K_Lap = KList4(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < Zero) IVal = IdxR
  end if
  ErrVal = ErLst4(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst4(K_Lap),' .'
end if

return

end subroutine DfineK
