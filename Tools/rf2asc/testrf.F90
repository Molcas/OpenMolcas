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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This program tests the runfile utilities.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
! Written: July 2003                                                   *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

program TestRF

use RunFile_data, only: TypDbl, TypInt, TypStr
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: Mxdata=64
real(kind=wp), parameter :: Step=0.75_wp, Thr=1.0e-12_wp
integer(kind=iwp) :: i, iOpt, iRc, iSeed, Loop, nDataA, nDataAx, nDataB, nDataBx, nDataC, nDataCx, nDataD, nDataDx, nDataE, &
                     nDataEx, nDataF, nDataFx, RecTypA, RecTypB, RecTypC, RecTypD, RecTypE, RecTypF
integer(kind=iwp) :: DataC(MxData), DataD(MxData)
real(kind=wp) :: DataA(MxData), DataB(MxData)
character :: DataE(MxData), DataF(MxData)
logical(kind=iwp) :: UseOld
real(kind=wp), external :: Random_Molcas

call IniMem()
call Init_LinAlg()
call PrgmInit('TestRF')
call NameRun('RUNFILE')

UseOld = .false.
if (.not. UseOld) then

  iSeed = 12345
  nDataA = 0
  nDataB = 0
  nDataC = 0
  nDataD = 0
  nDataE = 0
  nDataF = 0

  do i=1,MxData
    DataA(i) = real(i,kind=wp)
    DataC(i) = i
    DataE(i) = char(64+mod(i,48))
  end do
  DataB(:) = -DataA(:)
  DataD(:) = -DataC(:)
  DataF(:) = DataE(:)

  iRc = 0
  iOpt = 0
  do Loop=1,16
    if (Random_Molcas(iSeed) > Step) then
      nDataA = int((MxData-4)*Random_Molcas(iSeed)+4)
      write(u6,*) 'Adding A',nDataA
      call dWrRun('Amat',DataA,nDataA)
    end if
    if (Random_Molcas(iSeed) > Step) then
      nDataB = int((MxData-4)*Random_Molcas(iSeed)+4)
      write(u6,*) 'Adding B',nDataB
      call dWrRun('Bmat',DataB,nDataB)
    end if
    if (Random_Molcas(iSeed) > Step) then
      nDataC = int((MxData-4)*Random_Molcas(iSeed)+4)
      write(u6,*) 'Adding C',nDataC
      call iWrRun('Cmat',DataC,nDataC)
    end if
    if (Random_Molcas(iSeed) > Step) then
      nDataD = int((MxData-4)*Random_Molcas(iSeed)+4)
      write(u6,*) 'Adding D',nDataD
      call iWrRun('Dmat',DataD,nDataD)
    end if
    if (Random_Molcas(iSeed) > Step) then
      call NameRun('RUNXXX')
      i = int((MxData-4)*Random_Molcas(iSeed)/4)
      nDataE = 4*i+4
      write(u6,*) 'Adding E',nDataE
      call cWrRun('Emat',DataE,nDataE)
      call NameRun('RUNFILE')
    end if
    if (Random_Molcas(iSeed) > Step) then
      i = int((MxData-4)*Random_Molcas(iSeed)/4)
      nDataF = 4*i+4
      write(u6,*) 'Adding F',nDataF
      call cWrRun('Fmat',DataF,nDataF)
    end if
  end do

end if

DataA(:) = Zero
DataB(:) = Zero
DataC(:) = 0
DataD(:) = 0
DataE(:) = char(0)
DataF(:) = char(0)

nDataAx = nDataA
nDataBx = nDataB
nDataCx = nDataC
nDataDx = nDataD
nDataEx = nDataE
nDataFx = nDataF

nDataA = 0
nDataB = 0
nDataC = 0
nDataD = 0
nDataE = 0
nDataF = 0

call ffRun('Amat',nDataA,RecTypA)
call ffRun('Bmat',nDataB,RecTypB)
call ffRun('Cmat',nDataC,RecTypC)
call ffRun('Dmat',nDataD,RecTypD)
call NameRun('RUNXXX')
call ffRun('Emat',nDataE,RecTypE)
call NameRun('RUNFILE')
call ffRun('Fmat',nDataF,RecTypF)

if (UseOld) then
  nDataAx = nDataA
  nDataBx = nDataB
  nDataCx = nDataC
  nDataDx = nDataD
  nDataEx = nDataE
  nDataFx = nDataF
end if

if (nDataA /= nDataAx) then
  write(u6,*) 'nDataA:',nDataA,nDataAx
  stop
end if

if (nDataB /= nDataBx) then
  write(u6,*) 'nDataB:',nDataB,nDataBx
  stop
end if

if (nDataC /= nDataCx) then
  write(u6,*) 'nDataC:',nDataC,nDataCx
  stop
end if

if (nDataD /= nDataDx) then
  write(u6,*) 'nDataD:',nDataD,nDataDx
  stop
end if

if (nDataE /= nDataEx) then
  write(u6,*) 'nDataE:',nDataE,nDataEx
  stop
end if

if (nDataF /= nDataFx) then
  write(u6,*) 'nDataF:',nDataF,nDataFx
  stop
end if

if (RecTypA /= TypDbl) then
  write(u6,*) 'RecTypA:',RecTypA,TypDbl
  stop
end if

if (RecTypB /= TypDbl) then
  write(u6,*) 'RecTypB:',RecTypB,TypDbl
  stop
end if

if (RecTypC /= TypInt) then
  write(u6,*) 'RecTypC:',RecTypC,TypInt
  stop
end if

if (RecTypD /= TypInt) then
  write(u6,*) 'RecTypD:',RecTypD,TypInt
  stop
end if

if (RecTypE /= TypStr) then
  write(u6,*) 'RecTypE:',RecTypE,TypStr
  stop
end if

if (RecTypF /= TypStr) then
  write(u6,*) 'RecTypF:',RecTypF,TypStr
  stop
end if

if (nDataA > 0) then
  write(u6,*) 'Testing A',nDataA
  call dRdRun('Amat',DataA,nDataA)
  do i=1,nDataA
    if (abs(DataA(i)-i) > Thr) write(u6,*) 'A:',i,DataA(i)
  end do
end if

if (nDataB > 0) then
  write(u6,*) 'Testing B',nDataB
  call dRdRun('Bmat',DataB,nDataB)
  do i=1,nDataB
    if (abs(DataB(i)+i) > Thr) write(u6,*) 'B:',i,DataB(i)
  end do
end if

if (nDataC > 0) then
  write(u6,*) 'Testing C',nDataC
  call iRdRun('Cmat',DataC,nDataC)
  do i=1,nDataC
    if (DataC(i)-i /= 0) write(u6,*) 'C:',i,DataC(i)
  end do
end if

if (nDataD > 0) then
  write(u6,*) 'Testing D',nDataD
  call iRdRun('Dmat',DataD,nDataD)
  do i=1,nDataD
    if (DataD(i)+i /= 0) write(u6,*) 'D:',i,DataD(i)
  end do
end if

if (nDataE > 0) then
  call NameRun('RUNXXX')
  write(u6,*) 'Testing E',nDataE
  call cRdRun('Emat',DataE,nDataE)
  do i=1,nDataE
    if (ichar(DataE(i))-64-mod(i,48) /= 0) write(u6,*) 'E:',i,DataE(i)
  end do
  call NameRun('RUNFILE')
end if

if (nDataF > 0) then
  write(u6,*) 'Testing F',nDataF
  call cRdRun('Fmat',DataF,nDataF)
  do i=1,nDataF
    if (ichar(DataF(i))-64-mod(i,48) /= 0) write(u6,*) 'F:',i,DataF(i)
  end do
end if

call DumpRun(iRc,iOpt)

call GetMem('Finish','LIST','REAL',I,0)
call GetMem('Finish','TERM','REAL',I,0)

end program TestRF
