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

subroutine PtOkt0(H0,RR,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian,                   *
!                i.e., add quadrupole perturbation operator            *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nAtoms, nBas, nSym, ComStk, ComVal, Coor
use OneDat, only: sOpSiz
use Constants, only: Zero, One, Three, Five, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize, nTemp
real(kind=wp), intent(inout) :: H0(nSize), RR(nSize), Temp(nTemp)
character(len=8) :: Label
character(len=20) :: PriLbl
logical(kind=iwp) :: Exec, Orig
integer(kind=iwp) :: idum(1), iAtm, iCmp, iComp, iOpt1, iOpt2, iRc, iSyLbl, nInts
real(kind=wp) :: Alpha, Cntr(3), X, XOrig, Y, YOrig, Z, ZOrig
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Check if the command has been specified on input                 *
!                                                                      *
!----------------------------------------------------------------------*

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
  XOrig = Zero
  YOrig = Zero
  ZOrig = Zero
  if (ComStk(2,6,2,1)) XOrig = ComVal(2,6,2,1)
  if (ComStk(2,6,2,2)) YOrig = ComVal(2,6,2,2)
  if (ComStk(2,6,2,3)) ZOrig = ComVal(2,6,2,3)
  if (ComStk(2,6,2,4)) then
    iAtm = int(ComVal(2,6,2,4),kind=iwp)
    if (iAtm < 0 .or. iAtm > nAtoms) then
      write(u6,*) 'PtOkt0: You specified a invalid atom number as the origin of the perturbation operator.'
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
if (Debug) write(u6,'(6X,A,3F12.6)') 'Origin of the perturbation operator =',XOrig,YOrig,ZOrig

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
    iCmp = iComp
    Label = 'MltPl  3'
    PriLbl = 'MltPl  3; Comp =    '
    write(PriLbl(19:20),'(I2)') iComp
    iRc = -1
    iOpt1 = ibset(0,sOpSiz)
    iOpt2 = 0
    iSyLbl = 0
    call iRdOne(iRc,iOpt1,Label,iCmp,idum,iSyLbl)
    nInts = idum(1)
    if (iRc /= 0) call error()
    call RdOne(iRc,iOpt2,Label,iCmp,Temp,iSyLbl)
    if (iRc /= 0) call error()
    call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
    X = Temp(nInts+1)
    Y = Temp(nInts+2)
    Z = Temp(nInts+3)
    if (X /= XOrig .or. Y /= YOrig .or. Z /= ZOrig) then
      write(u6,*) 'PtOkt0: Input error, no matching center is found.'
      call Abend()
    end if
    Alpha = Five
    call DSCAL_(nInts+4,Alpha,Temp,1)
    if (iComp == 1 .or. iComp == 7 .or. iComp == 10) then
      Alpha = -Three
      call daxpy_(nInts,Alpha,RR,1,Temp,1)
      Temp(nInts+4) = Temp(nInts+4)+Alpha*RR(nInts+4)
    else if (iComp /= 5) then
      Alpha = -One
      call daxpy_(nInts,Alpha,RR,1,Temp,1)
      Temp(nInts+4) = Temp(nInts+4)+Alpha*RR(nInts+4)
    end if
    Alpha = Half*ComVal(2,6,1,iComp)
    call daxpy_(nInts,Alpha,Temp,1,H0,1)
    H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
    if (Debug) then
      write(u6,'(6X,A,F8.6)') 'weight =',ComVal(2,6,1,iComp)
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

contains

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*
subroutine error()

  write(u6,*) 'PtOkt0: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()

end subroutine error

end subroutine PtOkt0
