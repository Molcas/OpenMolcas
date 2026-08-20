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

subroutine Setup_RASSCF()
use Molcas, only: MxLev
use rasscf_global, only: NSM
use general_data, only: nActel, nElec3, nHole1, nRs1, nRs2, nRs3, nSym, NLEV, Level, &
                        NGAS, NGSSH, nRas,nRasEl,nRsPrt
use Definitions, only: iwp
implicit none
integer(kind=iwp) :: IGAS, iq, ISYM, nRs1T, NSTA

NLEV = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    NSTA = NLEV+1
    NLEV = NLEV+NGSSH(IGAS,ISYM)
    NSM(NSTA:NLEV) = ISYM
  end do
end do

Level(1:MxLev)=[(iq,iq=1,MxLev)]

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

end subroutine Setup_RASSCF
