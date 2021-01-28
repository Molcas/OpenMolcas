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

subroutine WR_VibRot_Info1(Lu,iOpt,iDisk,ntit1,J1A,J2A,lambda,n0,nvib1,Redm,Umax,Umin,ngrid,isn1,isn2,Req,xMass1,xMass2)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt
integer(kind=iwp), intent(inout) :: iDisk, ntit1, J1A, J2A, lambda, n0, nvib1, ngrid, isn1, isn2
real(kind=wp), intent(inout) :: Redm, Umax, Umin, Req, xMass1, xMass2
integer(kind=iwp) :: idum(1)
real(kind=wp) :: dum(1)

idum(1) = ntit1
call iDaFile(Lu,iOpt,idum,1,iDisk)
ntit1 = idum(1)

idum(1) = J1A
call iDaFile(Lu,iOpt,idum,1,iDisk)
J1A = idum(1)

idum(1) = J2A
call iDaFile(Lu,iOpt,idum,1,iDisk)
J2A = idum(1)

idum(1) = lambda
call iDaFile(Lu,iOpt,idum,1,iDisk)
lambda = idum(1)

idum(1) = n0
call iDaFile(Lu,iOpt,idum,1,iDisk)
n0 = idum(1)

idum(1) = nvib1
call iDaFile(Lu,iOpt,idum,1,iDisk)
nvib1 = idum(1)

dum(1) = Redm
call dDaFile(Lu,iOpt,dum,1,iDisk)
Redm = dum(1)

dum(1) = Umax
call dDaFile(Lu,iOpt,dum,1,iDisk)
Umax = dum(1)

dum(1) = Umin
call dDaFile(Lu,iOpt,dum,1,iDisk)
Umin = dum(1)

idum(1) = ngrid
call iDaFile(Lu,iOpt,idum,1,iDisk)
ngrid = idum(1)

idum(1) = isn1
call iDaFile(Lu,iOpt,idum,1,iDisk)
isn1 = idum(1)

idum(1) = isn2
call iDaFile(Lu,iOpt,idum,1,iDisk)
isn2 = idum(1)

dum(1) = Req
call dDaFile(Lu,iOpt,dum,1,iDisk)
Req = dum(1)

dum(1) = xMass1
call dDaFile(Lu,iOpt,dum,1,iDisk)
xMass1 = dum(1)

dum(1) = xMass2
call dDaFile(Lu,iOpt,dum,1,iDisk)
xMass2 = dum(1)

return

end subroutine WR_VibRot_Info1
