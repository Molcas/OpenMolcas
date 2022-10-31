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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iel, nDoF, nDim, ictl, Lu_10, iOff
real(kind=wp) :: EVal(nDim), EVec(2,nDoF,nDim), dDipM(nDim,iel), IRInt(nDim), RedM(nDim)
#include "Molcas.fh"
#include "real.fh"
#include "constants2.fh"
#include "WrkSpc.fh"
integer(kind=iwp), parameter :: Inc = 6
integer(kind=iwp) :: i, iHarm, iInt, iIRInt, ij, ipT, j, Jnc, k, l, nChDisp
real(kind=wp) :: Tmp(Inc)
character(len=LenIn6) :: ChDisp(3*MxAtom), Label !IFG
character(len=120) :: Line
character(len=80) :: format

call Get_iScalar('nChDisp',nChDisp)
call Get_cArray('ChDisp',ChDisp,LenIn6*nChDisp)

iIRInt = 0
do iHarm=1,nDim,Inc
  Jnc = min(Inc,nDim-iHarm+1)
  Label = ' '
  write(format,'(A,I3,A)') '(5X,A,1x,',Jnc,'(I7,3X))'
  write(u6,format) Label,(i,i=iHarm,iHarm+Jnc-1)
  write(u6,*)

  Label = 'Frequency:'
  write(format,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.2)'
  Line = ' '
  write(Line,format) Label,(EVal(i),i=iHarm,iHarm+Jnc-1)

  ! Replace minus signs with sign for imaginary unit.

  do i=1,120
    if (Line(i:i) == '-') Line(i:i) = 'i'
  end do
  write(u6,'(A)') Line
  write(u6,*)

  if (ictl /= 0) then
    Label = 'Intensity:'
    write(format,'(A,I3,A)') '(5X,A,1x,',Jnc,'ES10.3)'
    call dcopy_(Jnc,[Zero],0,Tmp,1)
    do k=1,Jnc
      do l=1,iel
        Tmp(k) = tmp(k)+dDipM(k+iHarm-1,l)**2
      end do
    end do
    write(u6,format) Label,(RF*tmp(i),i=1,Jnc)
    do i=1,Jnc
      iIRInt = iIRInt+1
      IRInt(iIRInt) = RF*tmp(i)
    end do
    Label = 'Red. mass:'
    write(format,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.5)'
    write(u6,format) Label,(RedM(i),i=iHarm,iHarm+Jnc-1)
    write(u6,*)
  else
    do i=1,Jnc
      iIRInt = iIRInt+1
      IRInt(iIRInt) = Zero
    end do
  end if

  write(format,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.5)'
  do iInt=1,nDoF
    write(u6,format) ChDisp(iInt+iOff),(EVec(1,iInt,i),i=iHarm,iHarm+Jnc-1)
  end do
  write(u6,*)
  write(u6,*)
end do

call GetMem('Temp','ALLO','REAL',ipT,nDim*nDoF)
ij = -1
do i=1,nDim
  do j=1,nDoF
    ij = ij+1
    Work(ipT+ij) = Evec(1,j,i)
  end do
end do
call WRH(Lu_10,1,[nDoF],[nDim],Work(ipT),EVAL,1,'*FREQUENCIES')
call GetMem('Temp','FREE','REAL',ipT,nDim*nDoF)

if (ictl /= 0) then
  !write(Lu_10,*) '*BEGIN PROJECTED DIPOLE TRANSITIONS'
  do j=1,iel
    call WRH(Lu_10,1,[ndim],[ndim],[Zero],dDipM(1,j),2,'*DIPOLE TRANSITIONS')
  end do
  !write(Lu_10,*) '*END PROJECTED DIPOLE TRANSITIONS'
end if

return

end subroutine GF_Print
