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
! Copyright (C) 1986, Per-Olof Widmark                                 *
!               1986, Bjorn O. Roos                                    *
!***********************************************************************
!***********************************************************************
!                                                                      *
! WRITTEN IN 1986 BY                                                   *
! PER-OLOF WIDMARK AND BJOERN O. ROOS                                  *
! DEPARTMENT OF THEORETICAL CHEMISTRY                                  *
! UNIVERSITY OF LUND                                                   *
! SWEDEN                                                               *
!                                                                      *
!***********************************************************************

subroutine vibrotmain(ireturn)

use Vibrot_globals, only: EoutO, ifPrWf, iobs, npoint, Titobs, Vibwvs, Vibwvs1, Vibwvs2
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: ncase, ngrid, nvib, i
real(kind=wp) :: Umin, Umax, E0, dE0, Redm, Teas, Req, sc, temp
real(kind=wp) :: R(npoint+4), PotR(npoint+4)

! Logical units
! Vibwvs can be saved for later Transition moment calculations
! where it will be redefined as Vibwvs1 or Vibwvs2:
Vibwvs = 12
Vibwvs1 = Vibwvs+1
Vibwvs2 = Vibwvs+2
call Daname(Vibwvs,'VIBWVS')
call Daname(Vibwvs1,'VIBWVS1')
call Daname(Vibwvs2,'VIBWVS2')
ncase = 0
call Vibinp(ncase,ngrid,nvib,Umin,Umax,R,PotR,E0,dE0,Redm,Teas,Req,sc,temp)
select case(ncase)
  case(1)
    call Vibrot(ngrid,nvib,Umin,Umax,R,PotR,E0,dE0,Redm,Req,sc,temp)
    if (IfPrWf > 0) call PrWf_VibRot(ngrid,R)
  case(2)
    do i=1,iobs
      call Vibtrm(ngrid,Umin,Umax,Teas,R,EoutO(1,i),Titobs(i))
    end do
  case default
    write(u6,1000)
    call Quit_OnUserError()
end select
call Daclos(Vibwvs)
call Daclos(Vibwvs1)
call Daclos(Vibwvs2)
ireturn = 0

return

1000 format(/1x,'No ROVIbrational or TRANsition keywordsspecified in input.')

end subroutine Vibrotmain
