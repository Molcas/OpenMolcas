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

subroutine PrInp_FFPT
!***********************************************************************
!                                                                      *
!     Objective: Write the title page on the standard output unit      *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "input.fh"
character*120 PrLine, BlLine, StLine
character*72 data, Line
character*8 Fmt1, Fmt2
character*4 Com, Sub1, Sub2, Parm
integer StrnLn
logical clear, nice

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!                                                                      *
!----------------------------------------------------------------------*

!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*
lLine = len(PrLine)
do i=1,lLine
  BlLine(i:i) = ' '
  StLine(i:i) = '*'
end do
lPaper = 132
left = (lPaper-lLine)/2
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
!     Print the project title                                          *
!----------------------------------------------------------------------*
if (mTit > 0) then
  write(6,*)
  nLine = mTit+5
  do i=1,nLine
    PrLine = BlLine
    if (i == 1 .or. i == nLine) PrLine = StLine
    if (i == 3) PrLine = 'Project:'
    if (i >= 4 .and. i <= nLine-2) then
      PrLine = Title(i-3)
    end if
    call Center(PrLine)
    write(6,Fmt1) '*'//PrLine//'*'
  end do
  write(6,*)
end if

!----------------------------------------------------------------------*
!     Print file identifier                                            *
!----------------------------------------------------------------------*

write(6,*)
write(6,'(6X,A)') 'Header of the ONEINT file:'
write(6,'(6X,A)') '--------------------------'
write(6,*)
write(Line,'(72A1)') (Header(i),i=1,72)
write(6,'(6X,A)') Line(:StrnLn(Line))
write(Line,'(72A1)') (Header(i),i=73,144)
write(6,'(6X,A)') Line(:StrnLn(Line))
write(6,*)

!----------------------------------------------------------------------*
!     Print the coordinates of the system                              *
!----------------------------------------------------------------------*

call PrCoor()

!----------------------------------------------------------------------*
!     Print comand list                                                *
!----------------------------------------------------------------------*

write(6,*)
write(6,'(6X,A)') 'The following tasks will be executed:'
write(6,'(6X,A)') '-------------------------------------'
write(6,*)
do iPr=1,120
  PrLine(iPr:iPr) = ' '
end do
Com = 'FFPT'
nSub1 = ComCtl(2,0,0)
clear = .false.
do iSub1=1,nSub1
  Sub1 = ComTab(2,iSub1,0,0)
  nSub2 = ComCtl(2,iSub1,0)
  nParm = 0
  if (nSub2 == 0) nParm = ComCtl(2,iSub1,1)
  ist = 25
  do iParm=0,nParm
    if (ComStk(2,iSub1,0,iParm)) then
      Parm = ComTab(2,iSub1,0,iParm)
      PrLine(1:4) = Com
      PrLine(9:12) = Sub1
      PrLine(ist:ist+3) = Parm
      ist = ist+4
      z = ComVal(2,iSub1,0,iParm)
      write(PrLine(ist:ist+7),'(F8.6)') z
      ist = ist+10
    end if
  end do
  if (PrLine(1:4) /= '    ') then
    if (clear) then
      write(6,'(14X,A)') PrLine(9:120)
    else
      write(6,'(6X,A)') PrLine
    end if
    do iPr=1,120
      PrLine(iPr:iPr) = ' '
    end do
    clear = .true.
  end if
  clear = .false.
  do iSub2=1,nSub2
    Sub2 = ComTab(2,iSub1,iSub2,0)
    nParm = ComCtl(2,iSub1,iSub2)
    ist = 25
    do iParm=0,nParm
      if (ComStk(2,iSub1,iSub2,iParm)) then
        Parm = ComTab(2,iSub1,iSub2,iParm)
        PrLine(1:4) = Com
        PrLine(9:12) = Sub1
        PrLine(17:20) = Sub2
        PrLine(ist:ist+3) = Parm
        ist = ist+4
        z = ComVal(2,iSub1,iSub2,iParm)
        write(data,'(F10.6)') z
        call LeftAd(data)
        iPoint = index(data,'.')
        nice = .true.
        do i=iPoint+1,Strnln(data)
          if (data(i:i) /= '0') nice = .false.
        end do
        if (nice) then
          PrLine(ist:ist+iPoint-1) = data(1:1+iPoint-2)
          ist = ist+iPoint+1
        else
          PrLine(ist:ist+iPoint+5) = data(1:1+iPoint+5)
          ist = ist+iPoint+7
        end if
      end if
    end do
    if (PrLine(1:4) /= '    ') then
      if (clear) then
        write(6,'(22X,A)') PrLine(17:120)
      else
        write(6,'(6X,A)') PrLine
      end if
      do iPr=1,120
        PrLine(iPr:iPr) = ' '
      end do
      clear = .true.
    end if
  end do
end do
do iTbl=1,mLbl
  write(6,'(6X,5A,I2,2A,F9.6)') 'GLBL    ','label="',gLblN(iTbl),'",','comp=',gLblC(iTbl),',','weight=',gLblW(iTbl)
end do
write(6,*)

!----------------------------------------------------------------------*
!     Terminate procedure                                              *
!----------------------------------------------------------------------*

return

end subroutine PrInp_FFPT
