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

use Real_Spherical, only: ipSph, RSph
use Index_Functions, only: nTri3_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: l, lMax, nij, nElem, nAtoms, iPrint
real(kind=wp), intent(in) :: rMP(nij,nElem), EC(3,nij), C_o_C(3)
real(kind=wp), intent(inout) :: xrMP(nij,nElem)
real(kind=wp), intent(out) :: Scratch_New(nij*(2+lMax+1)), Scratch_Org(nij*(2+lMax+1)), Cut_Off_Error
integer(kind=iwp) :: iAtom, iElem, iEnd, ij, iOff, iStrt, jAtom, k, kDim, m
real(kind=wp) :: Error, Estimated, Original, Percent, rms, rSum
character(len=80) :: Banner_Line
real(kind=wp), external :: DDot_

iEnd = nTri3_Elem1(lMax)
iStrt = nTri3_Elem1(l)+1
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
  write(u6,*)
  write(Banner_Line,'(A,I2)') 'Errors introduced by zeroing multipole moments greater than l = ',l
  call Banner(Banner_Line,1,80)
end if
rSum = Zero
iElem = nTri3_Elem1(l)+1
do k=l+1,lMax
  if (iPrint >= 1) then
    write(u6,*)
    write(u6,'(A,I1)') 'l=',k
    write(u6,*)
    write(u6,*) ' m     Original       New            Error            Percent'
    write(u6,*)
  end if

  kDim = (k+1)*(k+2)/2
  call DGEMM_('N','N',nij,2*k+1,kDim,One,xrMP(1,iElem),nij,RSph(ipSph(k)),kDim,Zero,Scratch_New,nij)
  call DGEMM_('N','N',nij,2*k+1,kDim,One,rMP(1,iElem),nij,RSph(ipSph(k)),kDim,Zero,Scratch_Org,nij)

  iOff = 1
  rms = Zero
  do m=-k,k
    Original = DDot_(nij,[One],0,Scratch_Org(iOff),1)
    Estimated = DDot_(nij,[One],0,Scratch_New(iOff),1)
    Error = Original-Estimated

    rSum = rSum+Error*Error
    rms = rms+Error*Error
    Percent = Zero
    if (abs(Error) < 1.0e-8_wp) then
      Percent = Zero
    else if (abs(Original) > 1.0e-13_wp) then
      Percent = abs(Error/Original)*100.0_wp
    else
      Percent = -huge(Percent)
    end if

    if (iPrint >= 1) then
      if (Percent >= Zero) then
        write(u6,'(I3,3F16.8,4X,F6.2)') m,Original,Estimated,Error,Percent
      else
        write(u6,'(I3,3F16.8,4X,A)') m,Original,Estimated,Error,'Infinite'
      end if
    end if

    iOff = iOff+nij
  end do
  if (iPrint >= 1) then
    rms = sqrt(rms/real(2*k+1,kind=wp))
    write(u6,*)
    write(u6,'(A,F16.8)') 'Root mean square = ',rms
  end if

  iElem = iElem+(k+1)*(k+2)/2
end do

Cut_Off_Error = rSum

return

end subroutine CutOff_Error
