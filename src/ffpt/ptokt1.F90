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

subroutine PtOkt1(Oper,Temp1,Temp2)
!***********************************************************************
!                                                                      *
!     Objective: Construct the perturbation operator of the form       *
!                <X*R**2> = <XXX>+<XYY>+<XZZ>                          *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nAtoms, nBas, nSym, ComStk, ComVal, Coor
use OneDat, only: sOpSiz
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character(len=2), intent(in) :: Oper
real(kind=wp), intent(inout) :: Temp1(*), Temp2(*)
character(len=8) :: Label
character(len=20) :: PriLbl
logical(kind=iwp) :: Orig
integer(kind=iwp) :: idum(1), iAtm, iComp, iOpt1, iOpt2, iRc, iSyLbl, jComp, nInts
real(kind=wp) :: Alpha, Cntr(3), X, XOrig, Y, YOrig, Z, ZOrig
integer(kind=iwp), parameter :: xrComp(3) = [1,4,6], yrComp(3) = [2,7,9], zrComp(3) = [3,8,10]
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!     Start procedure                                                  *
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
      write(u6,*) 'PtOkt1: You specified a invalid atom number as the origin of the perturbation operator.'
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

do iComp=1,3
  if (Oper == 'XR') jComp = xrComp(iComp)
  if (Oper == 'YR') jComp = yrComp(iComp)
  if (Oper == 'ZR') jComp = zrComp(iComp)
  Label = 'MltPl  3'
  PriLbl = 'MltPl  3; Comp =    '
  write(PriLbl(19:20),'(I2)') jComp
  iRc = -1
  iOpt1 = ibset(0,sOpSiz)
  iOpt2 = 0
  iSyLbl = 0
  call iRdOne(iRc,iOpt1,Label,jComp,idum,iSyLbl)
  nInts = idum(1)
  if (iRc /= 0) call error()
  call RdOne(iRc,iOpt2,Label,jComp,Temp1,iSyLbl)
  call CmpInt(Temp1,nInts,nBas,nSym,iSyLbl)
  X = Temp1(nInts+1)
  Y = Temp1(nInts+2)
  Z = Temp1(nInts+3)
  if (X /= XOrig .or. Y /= YOrig .or. Z /= ZOrig) then
    write(u6,*) 'PtOkt1: Input error, no matching center is found.'
    call Abend()
  end if
  if (iComp == 1) then
    call dcopy_(nInts,Temp1,1,Temp2,1)
    Temp2(nInts+4) = Temp1(nInts+4)
  else
    Alpha = One
    call daxpy_(nInts,Alpha,Temp1,1,Temp2,1)
    Temp2(nInts+4) = Temp2(nInts+4)+Alpha*Temp1(nInts+4)
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

  write(u6,*) 'PtOkt1: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()

end subroutine error

end subroutine PtOkt1
