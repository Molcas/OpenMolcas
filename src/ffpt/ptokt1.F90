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
!                <X*R**2>=<XXX>+<XYY>+<XZZ>                            *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "input.fh"
dimension Temp1(*), Temp2(*)
character*2 Oper
character*8 Label
character*20 PriLbl
dimension Cntr(3)
integer xrComp(3)
data xrComp/1,4,6/
integer yrComp(3)
data yrComp/2,7,9/
integer zrComp(3)
data zrComp/3,8,10/
logical Debug, Orig
data Debug/.false./
dimension idum(1)

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
!
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
      write(6,*) 'PtOkt1: You specified a invalid atom number as the origin of the perturbation operator.'
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

do iComp=1,3
  if (Oper == 'XR') jComp = xrComp(iComp)
  if (Oper == 'YR') jComp = yrComp(iComp)
  if (Oper == 'ZR') jComp = zrComp(iComp)
  Label = 'MltPl  3'
  PriLbl = 'MltPl  3; Comp =    '
  write(PriLbl(19:20),'(I2)') jComp
  iRc = -1
  iOpt1 = 1
  iOpt2 = 0
  iSyLbl = 0
  call iRdOne(iRc,iOpt1,Label,jComp,idum,iSyLbl)
  nInts = idum(1)
  if (iRc /= 0) goto 991
  call RdOne(iRc,iOpt2,Label,jComp,Temp1,iSyLbl)
  call CmpInt(Temp1,nInts,nBas,nSym,iSyLbl)
  X = Temp1(nInts+1)
  Y = Temp1(nInts+2)
  Z = Temp1(nInts+3)
  if (X /= XOrig .or. Y /= YOrig .or. Z /= ZOrig) then
    write(6,*) 'PtOkt1: Input error, no matching center is found.'
    call Abend()
  end if
  if (iComp == 1) then
    call dcopy_(nInts,Temp1,1,Temp2,1)
    Temp2(nInts+4) = Temp1(nInts+4)
  else
    Alpha = 1.0d0
    call daxpy_(nInts,Alpha,Temp1,1,Temp2,1)
    Temp2(nInts+4) = Temp2(nInts+4)+Alpha*Temp1(nInts+4)
  end if
end do

!----------------------------------------------------------------------*
!     Normal Exit                                                      *
!----------------------------------------------------------------------*

return

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*

991 write(6,*) 'PtOkt1: Error reading ONEINT'
write(6,'(A,A)') 'Label=',Label
call Abend()

end subroutine PtOkt1
