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

subroutine PtOkt0(H0,Ovlp,RR,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian,                   *
!                i.e., add quadrupole perturbation operator            *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "input.fh"
real*8 H0(nSize), Ovlp(nSize), RR(nSize), Temp(nTemp)
character*8 Label
character*20 PriLbl
dimension Cntr(3)
logical Debug, Exec, Orig
data Debug/.false./
dimension idum(1)

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Check if the command has been specified on input                 *
!                                                                      *
!----------------------------------------------------------------------*
!
Exec = .false.
Exec = Exec .or. ComStk(2,6,1,1)
Exec = Exec .or. ComStk(2,6,1,2)
Exec = Exec .or. ComStk(2,6,1,3)
Exec = Exec .or. ComStk(2,6,1,4)
Exec = Exec .or. ComStk(2,6,1,5)
Exec = Exec .or. ComStk(2,6,1,6)
Exec = Exec .or. ComStk(2,6,1,7)
Exec = Exec .or. ComStk(2,6,1,8)
Exec = Exec .or. ComStk(2,6,1,9)
Exec = Exec .or. ComStk(2,6,1,10)
if (.not. Exec) then
  return
end if

!----------------------------------------------------------------------*
!     Check if a origin has been specified                             *
!     The unspecified components of the origin are set to 0.0 !        *
!     If no origin has been given pick the center of mass!             *
!----------------------------------------------------------------------*

Orig = .false.
Orig = Orig .or. ComStk(2,6,2,1)
Orig = Orig .or. ComStk(2,6,2,2)
Orig = Orig .or. ComStk(2,6,2,3)
Orig = Orig .or. ComStk(2,6,2,4)

if (Orig) then
  XOrig = 0.0
  YOrig = 0.0
  ZOrig = 0.0
  if (ComStk(2,6,2,1)) XOrig = ComVal(2,6,2,1)
  if (ComStk(2,6,2,2)) YOrig = ComVal(2,6,2,2)
  if (ComStk(2,6,2,3)) ZOrig = ComVal(2,6,2,3)
  if (ComStk(2,6,2,4)) then
    iAtm = int(ComVal(2,6,2,4))
    if (iAtm < 0 .or. iAtm > nAtoms) then
      write(6,*) 'PtOkt0: You specified a invalid atom number as the origin of the perturbation operator.'
      call Abend()
    end if
    XOrig = Coor(1,iAtm)
    YOrig = Coor(2,iAtm)
    ZOrig = Coor(3,iAtm)
  end if
else
  call Get_dArray('Center of Mass',Cntr,3)
  XOrig = Cntr(1)
  YOrig = Cntr(2)
  ZOrig = Cntr(3)
end if
if (Debug) write(6,'(6X,A,3F12.6)') 'Origin of the perturbation operator =',XOrig,YOrig,ZOrig

!----------------------------------------------------------------------*
!     Loop over components                                             *
!----------------------------------------------------------------------*

do iComp=1,10
  if (ComStk(2,6,1,iComp)) then
    if (iComp == 1 .or. iComp == 4 .or. iComp == 6) then
      call PtOkt1('XR',Temp,RR)
    else if (iComp == 2 .or. iComp == 7 .or. iComp == 9) then
      call PtOkt1('YR',Temp,RR)
    else if (iComp == 3 .or. iComp == 8 .or. iComp == 10) then
      call PtOkt1('ZR',Temp,RR)
    end if
    Label = 'MltPl  3'
    PriLbl = 'MltPl  3; Comp =    '
    write(PriLbl(19:20),'(I2)') iComp
    iRc = -1
    iOpt1 = 1
    iOpt2 = 0
    iSyLbl = 0
    call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
    nInts = idum(1)
    if (iRc /= 0) goto 991
    call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
    if (iRc /= 0) goto 991
    call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
    X = Temp(nInts+1)
    Y = Temp(nInts+2)
    Z = Temp(nInts+3)
    if (X /= XOrig .or. Y /= YOrig .or. Z /= ZOrig) then
      write(6,*) 'PtOkt0: Input error, no matching center is found.'
      call Abend()
    end if
    Alpha = 5.0d0
    call DSCAL_(nInts+4,Alpha,Temp,1)
    if (iComp == 1 .or. iComp == 7 .or. iComp == 10) then
      Alpha = -3.0d0
      call daxpy_(nInts,Alpha,RR,1,Temp,1)
      Temp(nInts+4) = Temp(nInts+4)+Alpha*RR(nInts+4)
    else if (iComp /= 5) then
      Alpha = -1.0d0
      call daxpy_(nInts,Alpha,RR,1,Temp,1)
      Temp(nInts+4) = Temp(nInts+4)+Alpha*RR(nInts+4)
    end if
    Alpha = 0.5d0*ComVal(2,6,1,iComp)
    call daxpy_(nInts,Alpha,Temp,1,H0,1)
    H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
    if (Debug) then
      write(6,'(6X,A,F8.6)') 'weight =',ComVal(2,6,1,iComp)
      PriLbl = 'MltPl  3; Comp =    '
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
if (.false.) call Unused_real_array(Ovlp)

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*

991 write(6,*) 'PtOkt0: Error reading ONEINT'
write(6,'(A,A)') 'Label=',Label
call Abend()

end subroutine PtOkt0
