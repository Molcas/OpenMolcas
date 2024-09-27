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
subroutine Radlc(Zeta,nZeta,lsum,Rnr)
!***********************************************************************
!                                                                      *
! Object: to compute the radial part of the continuum multipole        *
!         integrals within an R-matrix run                             *
!                                                                      *
! Called from: Mltint                                                  *
!                                                                      *
! Author: K.Pfingst 21/5/93                                            *
!***********************************************************************

use fx, only: f_interface
use rmat, only: EpsAbs, EpsRel, ExpSum, KeyR, l, NagInt, QuadPack, RMatR, TestInt
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, lSum
real(kind=wp), intent(in) :: Zeta(nZeta)
real(kind=wp), intent(out) :: Rnr(nZeta,0:lsum)
integer(kind=iwp), parameter :: limit = 200, lenw = 4*limit
integer(kind=iwp) :: ier2, ir, iScrt(limit), iZeta, Last, nEval
real(kind=wp) :: AbsEr, reslt, Result2, Scrt(lenw)
#ifdef _DEBUGPRINT_
character(len=80) ::Label
#endif
procedure(f_interface) :: fradf

reslt = Zero
call Untested('Radlc')

!***********************************************************************

do ir=0,lsum
  do iZeta=1,nZeta
    expsum = Zeta(iZeta)
    l = ir
    !if (quadpack) then
    if (quadpack .and. (.not. testint)) then
      ier2 = -1
      call dqag(fradf,Zero,Rmatr,Epsabs,Epsrel,Keyr,result2,abser,neval,ier2,limit,lenw,last,iScrt,Scrt)
      if (ier2 /= 0) then
        call WarningMessage(1,' WARNING in Radlc; Consult the output for details!')
        write(u6,*)
        write(u6,*) ' WARNING in Radlc'
        write(u6,*)
        write(u6,*) ' ier=',ier2,' Error in Dqag called from Radlc.'
        write(u6,*) ' result=',result2
        write(u6,*) ' abser=',abser
        write(u6,*) ' neval=',neval
        write(u6,*)
      end if
      reslt = result2
    !else if (Nagint) then
    else if (Nagint .and. (.not. testint)) then
      call WarningMessage(2,'Radlc: Nagint option not implemented!')
      call Abend()
      !ier1 = -1
      !call d01ajf(fradf,Zero,Rmatr,Epsabs,Epsrel,result1,abser,wrk1,4*INtparm,iwrk1,INtparm,ier1)
      !if (ier1 /= 0) then
      !   write(u6,*)
      !   write(u6,*) ' WARNING in Radlc'
      !   write(u6,*)
      !   write(u6,*) ' ier=',ier1,' Error in D01ajf called from Radlc.'
      !   write(u6,*) ' result=',result1
      !   write(u6,*) ' abser=',abser
      !   write(u6,*) ' intparm=',intparm
      !   write(u6,*)
      !end if
      !reslt = result1
    end if
    if (testint) then
      call WarningMessage(2,'Radlc: testint option not implemented!')
      call Abend()
      !ier2 = -1
      !call dqag(fradf,Zero,Rmatr,Epsabs,Epsrel,Keyr,result2,abser,neval,ier2,limit,lenw,last,iScrt,Scrt)
      !if (ier2 /= 0) then
      !   write(u6,*)
      !   write(u6,*) ' WARNING in Radlc'
      !   write(u6,*)
      !   write(u6,*) ' ier=',ier2,' Error in Dqag called from Radlc.'
      !   write(u6,*) ' result=',result2
      !   write(u6,*) ' abser=',abser
      !   write(u6,*) ' neval=',neval
      !   write(u6,*)
      !end if

      !ier1 = -1
      !call d01ajf(fradf,Zero,Rmatr,Epsabs,Epsrel,result1,abser,wrk1,4*INtparm,iwrk1,INtparm,ier1)
      !if ((ier1 /= 0) .or. (ier2 /= 0)) then
      !   write(u6,*)
      !   write(u6,*) ' WARNING in Radlc'
      !   write(u6,*)
      !   write(u6,*) ' ier=',ier1,' Error in D01ajf called from Radlc.'
      !   write(u6,*) ' result=',result1
      !   write(u6,*) ' abser=',abser
      !   write(u6,*) ' intparm=',intparm
      !   write(u6,*)
      !end if
      !reslt = result1

      !diff = abs(result2-result1)
      !diff1 = abs((result2-result1)/result2)
      !if ((diff > epsabs) .or. (diff1 > epsrel)) then
      !   write(u6,*)
      !   write(u6,*) ' WARNING in Radlc'
      !   write(u6,*)
      !   write(u6,*) ' DIFFabs =',diff
      !   write(u6,*)
      !   write(u6,*) ' DIFFrel =',diff1
      !   write(u6,*)
      !   write(u6,*) ' NAG result=',result1
      !   write(u6,*)
      !   write(u6,*) ' NAG error =',ier1
      !   write(u6,*)
      !   write(u6,*) ' QUAD result=',result2
      !   write(u6,*)
      !   write(u6,*) ' QUAD error =',ier2
      !   write(u6,*)
      !end if
    end if
    Rnr(iZeta,ir) = reslt
  end do
end do

!***********************************************************************

#ifdef _DEBUGPRINT_
write(u6,*) ' Result in Radlc'
write(Label,'(A)') ' Rnr'
call RecPrt(Label,' ',Rnr(1,0),nZeta,lsum+1)
#endif

end subroutine Radlc
