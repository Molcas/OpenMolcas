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

!----------------------------------------------------------------------*
! A routine inspired by the wr_motra_info utility.                     *
!----------------------------------------------------------------------*
subroutine WrRdSim(iLu,iOpt,iDisk,iTcSim,nTcSim,Etot,Radie,nPart,Gmma,Gam,Esav)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iLu, iOpt, nTcSim
integer(kind=iwp), intent(inout) :: iDisk, iTcSim(nTcSim), nPart
real(kind=wp), intent(inout) :: Etot, Radie, Gmma, Gam, Esav
integer(kind=iwp) :: iDum(1)
real(kind=wp) :: Dum(1)

call iDaFile(iLu,iOpt,iTcSim,nTcSim,iDisk)
Dum(1) = Etot
call dDaFile(iLu,iOpt,Dum,1,iDisk)
Etot = Dum(1)
Dum(1) = Radie
call dDaFile(iLu,iOpt,Dum,1,iDisk)
Radie = Dum(1)
iDum(1) = nPart
call iDaFile(iLu,iOpt,iDum,1,iDisk)
nPart = iDum(1)
Dum(1) = Gmma
call dDaFile(iLu,iOpt,Dum,1,iDisk)
Gmma = Dum(1)
Dum(1) = Gam
call dDaFile(iLu,iOpt,Dum,1,iDisk)
Gam = Dum(1)
Dum(1) = Esav
call dDaFile(iLu,iOpt,Dum,1,iDisk)
Esav = Dum(1)

return

end subroutine WrRdSim
