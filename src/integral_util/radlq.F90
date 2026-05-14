!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993, Kurt Pfingst                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Radlq(Zeta,nZeta,lsum,Rnr,icop)
!***********************************************************************
!                                                                      *
! Object: to compute the radial part of the continuum Coulomb          *
!         integrals outside the R-matrix sphere                        *
!                                                                      *
! Called from: KneInt                                                  *
!                                                                      *
! Author: K.Pfingst 21/5/93                                            *
!***********************************************************************

use fx, only: f_interface
use rmat, only: EpsAbs, EpsRel, ExpSum, l, RMatR
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, lSum, icop
real(kind=wp), intent(in) :: Zeta(nZeta)
real(kind=wp), intent(out) :: Rnr(nZeta,0:lsum)
integer(kind=iwp), parameter :: limit = 200, lenw = 4*limit
integer(kind=iwp) :: ier, ir, iScrt(limit), iZeta, Last, nEval
real(kind=wp) :: AbsEr, reslt, Scrt(lenw)
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif
procedure(f_interface) :: fradf

call Untested('Radlq')

!                                                                      *
!***********************************************************************
!                                                                      *
do ir=0,lsum
  do iZeta=1,nZeta
    expsum = Zeta(iZeta)
    ier = 0
    l = ir-icop
    call dqagi(fradf,Rmatr,1,Epsabs,Epsrel,reslt,abser,neval,ier,limit,lenw,last,iScrt,Scrt)
    if (ier > 0) then
      call WarningMessage(1,' WARNING in Radlq; Consult output for details!')
      write(u6,*) ' ier=',ier,' Error in Dqagi called from Radlq.'
      write(u6,*) ' result=',reslt
      write(u6,*) ' abser =',abser
      write(u6,*) ' neval =',neval
      write(u6,*) ' WARNING in Radlq'
    end if
    Rnr(iZeta,ir) = reslt
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) ' Result in Radlq'
write(Label,'(A)') ' Rnr'
call RecPrt(Label,' ',Rnr(1,0),nZeta,lsum+1)
#endif

end subroutine Radlq
