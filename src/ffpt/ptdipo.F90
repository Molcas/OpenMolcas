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

subroutine PtDipo(H0,Ovlp,RR,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian,                   *
!                i.e., add dipole perturbation operator                *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "input.fh"
real*8 H0(nSize), Ovlp(nSize), RR(nSize), Temp(nTemp)
character*8 Label
character*20 PriLbl
logical Debug, Exec
data Debug/.false./
dimension idum(1)

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
!     Loop over all components and add constributions to the           *
!     one-electron Hamiltonian                                         *
!----------------------------------------------------------------------*

do iComp=1,3
  if (ComStk(2,1,1,iComp)) then
    Label = 'MltPl  1'
    iRc = -1
    iOpt1 = 1
    iOpt2 = 2
    iSyLbl = 0
    Alpha = -ComVal(2,1,1,iComp)
    call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
    nInts = idum(1)
    if (iRc /= 0) goto 991
    call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
    if (iRc /= 0) goto 991
    call CmpInt(Temp,nInts,nBas,nSym,iSyLbl)
    if (Debug) then
      write(6,'(6X,A,F8.6)') 'weight =',Alpha
      PriLbl = 'MltPl  1; Comp =    '
      write(PriLbl(19:20),'(I2)') iComp
      call PrDiOp(PriLbl,nSym,nBas,Temp)
      write(6,*) 'Nuclear contribution=',Temp(nInts+4)
    end if
    call daxpy_(nInts,Alpha,Temp,1,H0,1)
    H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
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

991 write(6,*) 'PtDipi: Error reading ONEINT'
write(6,'(A,A)') 'Label=',Label
call Abend()

end subroutine PtDipo
