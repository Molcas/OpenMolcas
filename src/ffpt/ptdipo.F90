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

subroutine PtDipo(H0,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian,                   *
!                i.e., add dipole perturbation operator                *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nBas, nSym, ComStk, ComVal
use OneDat, only: sNoOri, sOpSiz
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize, nTemp
real(kind=wp), intent(inout) :: H0(nSize), Temp(nTemp)
character(len=8) :: Label
character(len=20) :: PriLbl
logical(kind=iwp) :: Exec
integer(kind=iwp) :: idum(1), iCmp, iComp, iOpt1, iOpt2, iRc, iSyLbl, nInts
real(kind=wp) :: Alpha
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Check if the command has been specified on input                 *
!                                                                      *
!----------------------------------------------------------------------*

Exec = .false.
Exec = Exec .or. ComStk(2,1,1,1)
Exec = Exec .or. ComStk(2,1,1,2)
Exec = Exec .or. ComStk(2,1,1,3)
if (.not. Exec) then
  return
end if

!----------------------------------------------------------------------*
!     Loop over all components and add contributions to the            *
!     one-electron Hamiltonian                                         *
!----------------------------------------------------------------------*

do iComp=1,3
  if (ComStk(2,1,1,iComp)) then
    iCmp = iComp
    Label = 'MltPl  1'
    iRc = -1
    iOpt1 = ibset(0,sOpSiz)
    iOpt2 = ibset(0,sNoOri)
    iSyLbl = 0
    Alpha = -ComVal(2,1,1,iComp)
    call iRdOne(iRc,iOpt1,Label,iCmp,idum,iSyLbl)
    nInts = idum(1)
    if (iRc /= 0) call error()
    call RdOne(iRc,iOpt2,Label,iCmp,Temp,iSyLbl)
    if (iRc /= 0) call error()
    call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
    if (Debug) then
      write(u6,'(6X,A,F8.6)') 'weight =',Alpha
      PriLbl = 'MltPl  1; Comp =    '
      write(PriLbl(19:20),'(I2)') iComp
      call PrDiOp(PriLbl,nSym,nBas,Temp)
      write(u6,*) 'Nuclear contribution=',Temp(nInts+4)
    end if
    call daxpy_(nInts,Alpha,Temp,1,H0,1)
    H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
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

  write(u6,*) 'PtDipi: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()

end subroutine error

end subroutine PtDipo
