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

subroutine PtEfld(H0,Ovlp,RR,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian                    *
!                i.e., add electric field perturbation operator        *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "input.fh"
real*8 H0(nSize), Ovlp(nSize), RR(nSize), Temp(nTemp)
character*8 Label
character*20 PriLbl
logical Debug, Exec, Orig, NoCntr
data Debug/.false./
dimension idum(1)

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Check if the command has been specified on input                 *
!                                                                      *
!----------------------------------------------------------------------*

Exec = .false.
Exec = Exec .or. ComStk(2,3,1,1)
Exec = Exec .or. ComStk(2,3,1,2)
Exec = Exec .or. ComStk(2,3,1,3)
if (.not. Exec) then
  return
end if

!----------------------------------------------------------------------*
!     Check if a origin has been specified                             *
!     The unspecified components of the origin are set to 0.0 !        *
!----------------------------------------------------------------------*

Orig = .false.
Orig = Orig .or. ComStk(2,3,2,1)
Orig = Orig .or. ComStk(2,3,2,2)
Orig = Orig .or. ComStk(2,3,2,3)
Orig = Orig .or. ComStk(2,3,2,4)
if (.not. Orig) then
  write(6,*) 'PtElfd: No matching center is found.'
  call Abend()
end if

XOrig = 0.0
YOrig = 0.0
ZOrig = 0.0
if (ComStk(2,3,2,1)) XOrig = ComVal(2,3,2,1)
if (ComStk(2,3,2,2)) YOrig = ComVal(2,3,2,2)
if (ComStk(2,3,2,3)) ZOrig = ComVal(2,3,2,3)
if (ComStk(2,3,2,4)) then
  iAtm = int(ComVal(2,3,2,4))
  if (iAtm < 0 .or. iAtm > nAtoms) then
    write(6,*) 'PtEfld: You specified a invalid atom number as the origin of the perturbation operator.'
    call Abend()
  end if
  XOrig = Coor(1,iAtm)
  YOrig = Coor(2,iAtm)
  ZOrig = Coor(3,iAtm)
end if
if (Debug) write(6,'(6X,A,3F12.6)') 'Origin of perturbation operator =',XOrig,YOrig,ZOrig

!----------------------------------------------------------------------*
!     Loop over the max possible number of centers and                 *
!     search for coincidence in the origin definitions                 *
!----------------------------------------------------------------------*

MxCntr = 9999
NoCntr = .true.
do iCntr=1,MxCntr
  if (NoCntr) then
    Label = 'EF1     '
    write(Label(4:8),'(I5)') iCntr
    iRc = -1
    iOpt1 = 1
    iOpt2 = 2
    iSyLbl = 0
    do iComp=1,3
      call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
      if (iRc == 0) then
        nInts = idum(1)
        call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
        X = Temp(nInts+1)
        Y = Temp(nInts+2)
        Z = Temp(nInts+3)
        if (X == XOrig .and. Y == YOrig .and. Z == ZOrig) NoCntr = .false.
      end if
    end do
  end if
end do
if (NoCntr) then
  write(6,*) 'PtEfld: You missed to specify the origin of the operator.'
  call Abend()
end if
if (Debug) write(6,'(6X,A,A)') 'Label of perturbation operator =',Label

!----------------------------------------------------------------------*
!     If centers match read the integrals and accumulate contribution  *
!----------------------------------------------------------------------*

do iComp=1,3
  if (ComStk(2,3,1,iComp)) then
    iRc = -1
    iOpt1 = 1
    iOpt2 = 2
    iSyLbl = 0
    Alpha = -ComVal(2,3,1,iComp)
    call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
    nInts = idum(1)
    if (iRc /= 0) goto 991
    call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
    if (iRc /= 0) goto 991
    call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
    call daxpy_(nInts,Alpha,Temp,1,H0,1)
    H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
    if (Debug) then
      write(6,'(6X,A,F8.6)') 'weight =',Alpha
      PriLbl = '        ; Comp =    '
      PriLbl(1:8) = Label
      write(PriLbl(19:20),'(I2)') iComp
      call PrDiOp(PriLbl,nSym,nBas,Temp)
    end if
  end if
end do

!----------------------------------------------------------------------*
!     Normal Exit                                                      *
!----------------------------------------------------------------------*

return

! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Ovlp)
  call Unused_real_array(RR)
end if

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*

991 write(6,*) 'PtEfld: Error reading ONEINT'
write(6,'(A,A)') 'Label=',Label
call Abend()

end subroutine PtEfld
