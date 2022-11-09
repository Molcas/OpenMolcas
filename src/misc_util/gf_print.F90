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

subroutine GF_Print(EVal,EVec,dDipM,iel,nDoF,nDim,ictl,IRInt,RedM,Lu_10,iOff)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, RF
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iel, nDoF, nDim, ictl, Lu_10, iOff
real(kind=wp), intent(in) :: EVal(nDim), EVec(2,nDoF,nDim), dDipM(nDim,iel), RedM(nDim)
real(kind=wp), intent(out) :: IRInt(nDim)
#include "Molcas.fh"
integer(kind=iwp), parameter :: Inc = 6
integer(kind=iwp) :: i, iHarm, iInt, iIRInt, j, Jnc, l, nChDisp
real(kind=wp) :: Tmp(Inc)
character(len=LenIn6) :: Label
character(len=120) :: Line
character(len=80) :: frmt
real(kind=wp), allocatable :: T(:,:)
character(len=LenIn6), allocatable :: ChDisp(:)

call Get_iScalar('nChDisp',nChDisp)
call mma_allocate(ChDisp,nChDisp,label='ChDisp')
call Get_cArray('ChDisp',ChDisp,LenIn6*nChDisp)

iIRInt = 0
do iHarm=1,nDim,Inc
  Jnc = min(Inc,nDim-iHarm+1)
  Label = ' '
  write(frmt,'(A,I3,A)') '(5X,A,1x,',Jnc,'(I7,3X))'
  write(u6,frmt) Label,(i,i=iHarm,iHarm+Jnc-1)
  write(u6,*)

  Label = 'Frequency:'
  write(frmt,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.2)'
  Line = ' '
  write(Line,frmt) Label,(EVal(i),i=iHarm,iHarm+Jnc-1)

  ! Replace minus signs with sign for imaginary unit.

  do i=1,120
    if (Line(i:i) == '-') Line(i:i) = 'i'
  end do
  write(u6,'(A)') Line
  write(u6,*)

  if (ictl /= 0) then
    Label = 'Intensity:'
    write(frmt,'(A,I3,A)') '(5X,A,1x,',Jnc,'ES10.3)'
    Tmp(1:Jnc) = Zero
    do l=1,iel
      Tmp(1:Jnc) = Tmp(1:Jnc)+dDipM(iHarm:iHarm+Jnc-1,l)**2
    end do
    write(u6,frmt) Label,(RF*Tmp(i),i=1,Jnc)
    IRInt(iIRInt+1:iIRInt+Jnc) = RF*Tmp(1:Jnc)
    iIRInt = iIRInt+Jnc
    Label = 'Red. mass:'
    write(frmt,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.5)'
    write(u6,frmt) Label,(RedM(i),i=iHarm,iHarm+Jnc-1)
    write(u6,*)
  else
    IRInt(iIRInt+1:iIRInt+Jnc) = Zero
    iIRInt = iIRInt+Jnc
  end if

  write(frmt,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.5)'
  do iInt=1,nDoF
    write(u6,frmt) ChDisp(iInt+iOff),(EVec(1,iInt,i),i=iHarm,iHarm+Jnc-1)
  end do
  write(u6,*)
  write(u6,*)
end do
call mma_deallocate(ChDisp)

call mma_allocate(T,nDoF,nDim,label='Temp')
T(:,:) = Evec(1,:,:)
Line = '*FREQUENCIES'
call WRH(Lu_10,1,[nDoF],[nDim],T,EVAL,1,Line)
call mma_deallocate(T)

if (ictl /= 0) then
  !write(Lu_10,*) '*BEGIN PROJECTED DIPOLE TRANSITIONS'
  Line = '*DIPOLE TRANSITIONS'
  do j=1,iel
    call WRH(Lu_10,1,[nDim],[nDim],[Zero],dDipM(1,j),2,Line)
  end do
  !write(Lu_10,*) '*END PROJECTED DIPOLE TRANSITIONS'
end if

return

end subroutine GF_Print
