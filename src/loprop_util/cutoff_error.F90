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

subroutine CutOff_Error(l,lMax,rMP,xrMP,nij,EC,C_o_C,nElem,Scratch_New,Scratch_Org,nAtoms,iPrint,Cut_Off_Error)

use Real_Spherical

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "status.fh"
real*8 rMP(nij,nElem), xrMP(nij,nElem), EC(3,nij), C_o_C(3)
real*8 Scratch_New(nij*(2+lMax+1)), Scratch_Org(nij*(2+lMax+1))
character*80 Banner_Line
! Statement function
mElem(i) = (i+1)*(i+2)*(i+3)/6

iEnd = mElem(lMax)
iStrt = mElem(l)+1
ij = 0
do iAtom=1,nAtoms
  do jAtom=1,iAtom
    ij = ij+1
    call ReExpand(xrMP,nij,nElem,C_o_C,EC(1,ij),ij,lMax)

    ! Set all elements corresponding to l+1 = 0
    do iElem=iStrt,iEnd
      xrMP(ij,iElem) = Zero
    end do

    call ReExpand(xrMP,nij,nElem,EC(1,ij),C_o_C,ij,lMax)
  end do
end do

if (iPrint >= 1) then
  write(6,*)
  write(Banner_Line,'(A,I2)') 'Errors introduced by zeroing multipole moments greater than l = ',l
  call Banner(Banner_Line,1,80)
end if
Sum = 0.0d0
iElem = mElem(l)+1
do k=l+1,lMax
  if (iPrint >= 1) then
    write(6,*)
    write(6,'(A,I1)') 'l=',k
    write(6,*)
    write(6,*) ' m     Original       New            Error            Percent'
    write(6,*)
  end if

  kDim = (k+1)*(k+2)/2
  call DGEMM_('N','N',nij,2*k+1,kDim,1.0d0,xrMP(1,iElem),nij,RSph(ipSph(k)),kDim,0.0d0,Scratch_New,nij)
  call DGEMM_('N','N',nij,2*k+1,kDim,1.0d0,rMP(1,iElem),nij,RSph(ipSph(k)),kDim,0.0d0,Scratch_Org,nij)

  iOff = 1
  rms = 0.0d0
  do m=-k,k
    Original = DDot_(nij,[One],0,Scratch_Org(iOff),1)
    Estimated = DDot_(nij,[One],0,Scratch_New(iOff),1)
    Error = Original-Estimated

    Sum = Sum+Error*Error
    rms = rms+Error*Error
    Percent = 0.0d0
    if (abs(Error) < 1.0D-8) then
      Percent = 0.0d0
    else if (abs(Original) > 1.0D-13) then
      Percent = abs(Error/Original)*100.0d0
    else
      Percent = -999.99d0
    end if

    if (iPrint >= 1) then
      if (Percent >= Zero) then
        write(6,'(I3,3F16.8,4X,F6.2)') m,Original,Estimated,Error,Percent
      else
        write(6,'(I3,3F16.8,4X,A)') m,Original,Estimated,Error,'Infinite'
      end if
    end if

    iOff = iOff+nij
  end do
  if (iPrint >= 1) then
    rms = sqrt(rms/dble(2*k+1))
    write(6,*)
    write(6,'(A,F16.8)') 'Root mean square = ',rms
  end if

  iElem = iElem+(k+1)*(k+2)/2
end do

Cut_Off_Error = Sum

return

end subroutine CutOff_Error
