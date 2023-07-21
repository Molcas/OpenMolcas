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

use ReMez_mod

implicit real*8(A-H,O-Z)
parameter(MxList=31,MxK=20,ZERO=0.0D+00)
character*8 Milli, Micro, Nano, Pico, Demand
real*8 ErLst1(MxK), ErLst2(MxK), ErLst3(MxK), ErLst4(MxK), ErLst5(MxK), RList(30)
integer IdxMin(MxK), KList1(MxList), KList2(MxList), KList3(MxList), KList4(MxList)
data Milli,Micro/'MILLI   ','MICRO   '/
data Nano,Pico/'NANO    ','PICO    '/
data IdxMin/1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,4,6,7,8/
data KList1/2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/
data KList2/3,4,5,6,6,7,7,7,7,7,7,7,8,8,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10/
data KList3/5,7,8,9,10,10,11,11,11,11,12,12,13,14,14,14,15,15,15,15,15,16,17,17,17,18,18,18,18,18,18/
data KList4/6,10,10,12,13,14,14,15,15,16,16,16,18,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20/
data ErLst1/0.000D+00,8.752D-03,5.052D-03,1.700D-03,6.428D-04,2.646D-04,1.163D-04,5.392D-05,2.611D-05,1.312D-05,6.807D-06, &
            3.630D-06,1.984D-06,1.108D-06,6.311D-07,3.659D-07,2.155D-07,1.289D-07,7.811D-08,4.794D-08/
data ErLst2/0.000D+00,0.000D+00,0.000D+00,6.258D-06,4.243D-06,7.379D-06,9.841D-06,8.724D-06,7.468D-06,9.296D-06,3.844D-06, &
            1.582D-06,6.481D-07,2.648D-07,1.079D-07,4.388D-08,1.780D-08,7.213D-09,2.918D-09,1.179D-09/
data ErLst3/0.000D+00,0.000D+00,0.000D+00,0.000D+00,0.000D+00,7.741D-09,1.716D-09,4.045D-09,6.822D-09,9.323D-09,3.357D-09, &
            4.301D-09,7.587D-09,7.555D-09,7.134D-09,8.484D-09,7.213D-09,8.602D-09,8.602D-09,9.700D-09/
data ErLst4/0.000D+00,0.000D+00,0.000D+00,0.000D+00,0.000D+00,0.000D+00,0.000D+00,0.000D+00,0.000D+00,9.021D-12,6.492D-13, &
            5.389D-12,5.735D-12,9.405D-12,6.723D-12,5.298D-12,1.050D-12,3.128D-12,7.099D-12,9.371D-12/
data ErLst5/8.556D-02,1.785D-02,5.052D-03,1.700D-03,6.428D-04,2.646D-04,1.163D-04,5.392D-05,2.611D-05,1.312D-05,6.807D-06, &
            3.630D-06,1.984D-06,1.108D-06,6.311D-07,3.659D-07,2.155D-07,1.289D-07,7.811D-08,4.794D-08/
data RList/2.0D+00,5.0D+00,1.0D+01,2.0D+01,3.0D+01,4.0D+01,5.0D+01,6.0D+01,7.0D+01,8.0D+01,9.0D+01,1.0D+02,2.0D+02,3.0D+02, &
           4.0D+02,5.0D+02,6.0D+02,7.0D+02,8.0D+02,9.0D+02,1.0D+03,2.0D+03,3.0D+03,4.0D+03,5.0D+03,6.0D+03,7.0D+03,8.0D+03, &
           9.0D+03,1.0D+04/
logical Trial

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
      if (R <= 3.0D+04) then
        K_Lap = 19
        ErrVal = ErLst3(14)
        write(IW,'(/A,E11.4E2)') ' This K guarantees the error less than ',ErrVal
        goto 100
      else if (R <= 1.0D+05) then
        ErrVal = ErLst3(15)
        write(IW,'(/A,E11.4E2)') ' This K guarantees the error less than ',ErrVal
        goto 100
      end if
    end if
    write(IW,'(/A)') '!!! Caution !!!'
    write(IW,'(A,E11.4E2,A)') 'In this R value, we can only guarantee',ErrVal,' accuracy.'
  end if
100 continue
  return
end if

if (InitR <= 8) Trial = .true.

IVal = InitR+1

if (Demand == Milli) then
  K_Lap = KList1(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < ZERO) IVal = IdxR
  end if
  ErrVal = ErLst1(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst1(K_Lap),' .'
else if (Demand == Micro) then
  K_Lap = KList2(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < ZERO) IVal = IdxR
  end if
  ErrVal = ErLst2(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst2(K_Lap),' .'
else if (Demand == Nano) then
  K_Lap = KList3(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < ZERO) IVal = IdxR
  end if
  ErrVal = ErLst3(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst3(K_Lap),' .'
else if (Demand == Pico) then
  K_Lap = KList4(IVal)
  if (Trial) then
    IdxR = IdxMin(K_Lap)
    R_Val = RList(IdxR)
    if ((R-R_Val) < ZERO) IVal = IdxR
  end if
  ErrVal = ErLst4(K_Lap)
  write(IW,'(/A,E11.4E2,A)') ' This K guarantees the error less than ',ErLst4(K_Lap),' .'
end if

return

end subroutine DfineK
