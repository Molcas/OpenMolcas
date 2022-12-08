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

subroutine PtGLbl(H0,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian,                   *
!                i.e., add any perturbation defined by the GLBL input  *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: mLbl, nBas, nSym, ComStk, gLblN, gLblC, gLblW
use OneDat, only: sNoOri, sOpSiz
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize, nTemp
real(kind=wp), intent(inout) :: H0(nSize), Temp(nTemp)
character(len=8) :: Label
character(len=20) :: PriLbl
logical(kind=iwp) :: Exec
integer(kind=iwp) :: idum(1), iComp, iLbl, iOpt1, iOpt2, iRc, iSyLbl, nInts
real(kind=wp) :: Alpha
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Check if the command has been specified on input                 *
!                                                                      *
!----------------------------------------------------------------------*

Exec = .false.
Exec = Exec .or. ComStk(3,0,0,0)
if (.not. Exec) then
  return
end if

!----------------------------------------------------------------------*
!     Loop over all components and add constributions to the           *
!     one-electron Hamiltonian                                         *
!----------------------------------------------------------------------*

do iLbl=1,mLbl
  Label = gLblN(iLbl)
  iComp = gLblC(iLbl)
  Alpha = gLblW(iLbl)
  if (Debug) then
    write(u6,'(6X,5A,I2,2A,G9.2)') 'GLBL    ','label  ="',gLblN(iLbl),'",','comp   =',gLblC(iLbl),',','weight =',gLblW(iLbl)
  end if
  iRc = -1
  iOpt1 = ibset(0,sOpSiz)
  iOpt2 = ibset(0,sNoOri)
  iSyLbl = 0
  call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
  nInts = idum(1)
  if (iRc /= 0) call error()
  call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
  if (iRc /= 0) call error()
  call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
  if (Debug) then
    PriLbl = Label//'; Comp =    '
    write(PriLbl(19:20),'(I2)') iComp
    call PrDiOp(PriLbl,nSym,nBas,Temp)
  end if
  call daxpy_(nInts,Alpha,Temp,1,H0,1)
  H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
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

  write(u6,*) 'PtGlbl: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()

end subroutine error

end subroutine PtGLbl
