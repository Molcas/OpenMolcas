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

subroutine SG_Setup_MCLR()

use molcas, only: MxLev
use sguga, only: SG_Init_Simple
use general_data, only: iSpin, nActEl, nElec3, nHole1, nRS1, nRS2, nRS3, nSym
use general_data, only: nRas, nRasEl, nRsPrt
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iBas, iq, ISM(1:MxLev), iSym, Level(MxLev), nLev, nRs1T
integer(kind=iwp), parameter :: iState=1

nLev = 0
do iSym=1,nSym
  do iBas=1,nRs1(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
do iSym=1,nSym
  do iBas=1,nRs2(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
do iSym=1,nSym
  do iBas=1,nRs3(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do

if (nHole1+nElec3 /= 0) then
   nRsPrt=3
   nRas(:,1)=nRs1(:)
   nRas(:,2)=nRs2(:)
   nRas(:,3)=nRs3(:)
  nRs1T = sum(nRs1(1:nSym))
   nRasEl(1)=2*nRs1T-nHole1
   nRasEl(2)=nActel-nElec3
   nRasEl(3)=nActel
else
   nRsPrt=1
   nRas(:,1)=nRs2(:)
   nRasEl(1)=nActel
end if

Level(1:MxLev)=[(iq,iq=1,MxLev)]

Call SG_Init_Simple(istate,nSym,nActEl,iSpin,     &
                    nRas,nRasEl,nRsPrt,            &
                    xLevel=Level, xL2Act=Level,    &
                    xNLEV=nLev, xNSM=ISM)

end subroutine SG_Setup_MCLR
