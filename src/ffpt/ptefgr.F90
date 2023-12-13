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

subroutine PtEfGr(H0,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian,                   *
!                i.e., add electric field gradient pert. operator      *
!                Note: This operator is already in traceless form      *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nAtoms, nBas, nSym, ComStk, ComVal, Coor
use OneDat, only: sNoOri, sOpSiz
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize, nTemp
real(kind=wp), intent(inout) :: H0(nSize), Temp(nTemp)
character(len=8) :: Label
character(len=20) :: PriLbl
logical(kind=iwp) :: Exec, Orig, NoCntr
integer(kind=iwp) :: idum(1), iAtm, iCmp, iCntr, iComp, iOpt1, iOpt2, iRc, iSyLbl, jComp, MxCntr, nInts
real(kind=wp) :: Alpha, X, XOrig, Y, YOrig, Z, ZOrig
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Check if the command has been specified on input                 *
!                                                                      *
!----------------------------------------------------------------------*

Exec = .false.
Exec = Exec .or. ComStk(2,4,1,1)
Exec = Exec .or. ComStk(2,4,1,2)
Exec = Exec .or. ComStk(2,4,1,3)
Exec = Exec .or. ComStk(2,4,1,4)
Exec = Exec .or. ComStk(2,4,1,5)
Exec = Exec .or. ComStk(2,4,1,6)
if (.not. Exec) then
  return
end if

!----------------------------------------------------------------------*
!     Check if a origin has been specified                             *
!     The unspecified components of the origin are set to 0.0 !        *
!----------------------------------------------------------------------*

Orig = .false.
Orig = Orig .or. ComStk(2,4,2,1)
Orig = Orig .or. ComStk(2,4,2,2)
Orig = Orig .or. ComStk(2,4,2,3)
Orig = Orig .or. ComStk(2,4,2,4)
if (.not. Orig) then
  write(u6,*) 'PtEfGr: Input error, no matching center found.'
  call Abend()
end if

XOrig = Zero
YOrig = Zero
ZOrig = Zero
if (ComStk(2,4,2,1)) XOrig = ComVal(2,4,2,1)
if (ComStk(2,4,2,2)) YOrig = ComVal(2,4,2,2)
if (ComStk(2,4,2,3)) ZOrig = ComVal(2,4,2,3)
if (ComStk(2,4,2,4)) then
  iAtm = int(ComVal(2,4,2,4),kind=iwp)
  if (Debug) write(u6,'(6X,A,I2)') 'Origin of perturbation is centered at atom ',iAtm
  if (iAtm < 0 .or. iAtm > nAtoms) then
    write(u6,*) 'PrEfGr: iAtm < 0 .or. iAtm > nAtoms'
    write(u6,*) 'iAtm,nAtoms=',iAtm,nAtoms
    write(u6,*) 'You specified a invalid atom number as the origin of the perturbation operator.'
    call Abend()
  end if
  XOrig = Coor(1,iAtm)
  YOrig = Coor(2,iAtm)
  ZOrig = Coor(3,iAtm)
end if
if (Debug) write(u6,'(6X,A,3F12.6)') 'Origin of perturbation operator =',XOrig,YOrig,ZOrig

!----------------------------------------------------------------------*
!     Loop over the max possible number of centers and                 *
!     search for coincidence in the origin definitions                 *
!----------------------------------------------------------------------*

MxCntr = 9999
NoCntr = .true.
do iCntr=1,MxCntr
  if (NoCntr) then
    Label = 'EF2     '
    write(Label(4:8),'(I5)') iCntr
    iRc = -1
    iOpt1 = ibset(0,sOpSiz)
    iOpt2 = ibset(0,sNoOri)
    iSyLbl = 0
    do iComp=1,6
      iCmp = iComp
      call iRdOne(iRc,iOpt1,Label,iCmp,idum,iSyLbl)
      if (iRc == 0) then
        nInts = idum(1)
        call RdOne(iRc,iOpt2,Label,iCmp,Temp,iSyLbl)
        X = Temp(nInts+1)
        Y = Temp(nInts+2)
        Z = Temp(nInts+3)
        if (X == XOrig .and. Y == YOrig .and. Z == ZOrig) NoCntr = .false.
      end if
    end do
  end if
end do
if (NoCntr) then
  write(u6,*) 'PrEfGr: No center found!'
  write(u6,*) 'XOrig,YOrig,ZOrig=',XOrig,YOrig,ZOrig
  call Abend()
end if
if (Debug) write(u6,'(6X,A,A)') 'Label of perturbation operator =',Label

!----------------------------------------------------------------------*
!     Read all components and reduce the diagonal ones                 *
!----------------------------------------------------------------------*

do iComp=1,6
  if (ComStk(2,4,1,iComp)) then
    iRc = -1
    iOpt1 = ibset(0,sOpSiz)
    iOpt2 = ibset(0,sNoOri)
    iSyLbl = 0
    Alpha = -ComVal(2,4,1,iComp)

    if (iComp /= 6) then
      iCmp = iComp
      call iRdOne(iRc,iOpt1,Label,iCmp,idum,iSyLbl)
      nInts = idum(1)
      if (iRc /= 0) call error()
      call RdOne(iRc,iOpt2,Label,iCmp,Temp,iSyLbl)
      if (iRc /= 0) call error()
      call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
      call daxpy_(nInts,Alpha,Temp,1,H0,1)
      H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
    else
      ! Note that we have not explicit storage of the ZZ term!
      ! It has to be computed as -XX-YY.
      jComp = 1  ! the XX-term
      call iRdOne(iRc,iOpt1,Label,jComp,idum,iSyLbl)
      nInts = idum(1)
      if (iRc /= 0) call error()
      call RdOne(iRc,iOpt2,Label,jComp,Temp,iSyLbl)
      if (iRc /= 0) call error()
      call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
      call daxpy_(nInts,-Alpha,Temp,1,H0,1)
      H0(nInts+4) = H0(nInts+4)+Alpha*Temp(nInts+4)

      jComp = 4  ! the YY-term
      call iRdOne(iRc,iOpt1,Label,jComp,idum,iSyLbl)
      nInts = idum(1)
      if (iRc /= 0) call error()
      call RdOne(iRc,iOpt2,Label,jComp,Temp,iSyLbl)
      if (iRc /= 0) call error()
      call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
      call daxpy_(nInts,-Alpha,Temp,1,H0,1)
      H0(nInts+4) = H0(nInts+4)+Alpha*Temp(nInts+4)
    end if
    if (Debug) then
      write(u6,'(6X,A,F8.6)') 'weight =',ComVal(2,4,1,iComp)
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

contains

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*
subroutine error()

  write(u6,*) 'PtEfGr: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()

end subroutine error

end subroutine PtEfGr
