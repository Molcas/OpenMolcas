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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine makedx_cvb(dx,nparm,ioptc,heigvec,heigval,dxp,gradp,w2,close2conv,nposeig,scalesmall,wrongstat,nnegeig,opth,alfastart, &
                      alfa)

use casvb_global, only: alftol, cnrm, cnrmtol, eigwrngtol, exp12tol, expct, form2AF, formAD, formAF, grdwrngtol, hh, ip
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nparm, nposeig, nnegeig
real(kind=wp), intent(out) :: dx(nparm), w2(nparm), alfa
integer(kind=iwp), intent(inout) :: ioptc
real(kind=wp), intent(in) :: heigvec(nparm,nparm), heigval(nparm), gradp(nparm), alfastart
real(kind=wp), intent(inout) :: dxp(nparm)
logical(kind=iwp), intent(in) :: close2conv, scalesmall, wrongstat, opth
integer(kind=iwp) :: i
real(kind=wp) :: exp1, exp2
real(kind=wp), external :: dnrm2_

alfa = alfastart
! << Update is small >>
if ((cnrm < hh) .and. scalesmall) then
  ! Close to wrong stationary point:
  if (wrongstat) then
    if (dnrm2_(nparm,gradp,1) < grdwrngtol) then
      write(u6,*) ' Gradient too small - not using information!'
      w2(:) = Zero
      do i=1,nnegeig
        if (heigval(i) >= eigwrngtol) w2(i) = sign(One,gradp(i))
      end do
      do i=nnegeig+1,nparm
        if (heigval(i) <= -eigwrngtol) w2(i) = sign(One,gradp(i))
      end do
      call getdxp_cvb(dxp,w2,heigval,nnegeig,nparm,alfa)
      cnrm = dnrm2_(nparm,dxp,1)
    end if
  else
    if ((.not. opth) .and. (ip >= 2)) write(u6,form2AF) ' Scaling update from :',cnrm,' to :',hh
  end if
  dxp(:) = hh/cnrm*dxp(:)
  cnrm = hh
else if (cnrm >= hh) then
  call optalf_cvb(heigval,gradp,nparm,hh,alfa,nnegeig,alfastart,alftol)
  call getdxp_cvb(dxp,gradp,heigval,nnegeig,nparm,alfa)
  call expec_cvb(dxp,gradp,heigval,nnegeig,nparm,expct,exp1,exp2)
  cnrm = dnrm2_(nparm,dxp,1)
  if ((.not. opth) .and. (ip >= 2)) write(u6,formAF) ' Alpha and norm of update :',alfa,cnrm
end if

if ((ioptc > 0) .and. ((.not. opth) .and. (cnrm < cnrmtol))) then
  if (ip >= 0) then
    write(u6,'(a)') ' '
    write(u6,formAD) ' WARNING - predicted update too small :',cnrm,cnrmtol
  end if
  ioptc = -2
  return
end if
do
  call expec_cvb(dxp,gradp,heigval,nnegeig,nparm,expct,exp1,exp2)
  if ((exp1 >= -exp12tol) .and. (exp2 <= exp12tol)) exit
  dxp(:) = 0.9_wp*dxp(:)
  cnrm = dnrm2_(nparm,dxp,1)
  if (cnrm < cnrmtol) then
    write(u6,formAD) ' Norm of update too small :',cnrm,cnrmtol
    call abend_cvb()
  end if
end do
if ((ip >= 2) .and. close2conv .and. ((exp1 < Zero) .or. (exp2 > Zero))) then
  write(u6,*) ' Warning - not a max/min direction !'
  if (nnegeig > 0) write(u6,*) ' Expected change for maximized variables :',exp1
  if (nposeig > 0) write(u6,*) ' Expected change for minimized variables :',exp2
end if
call mxatb_cvb(heigvec,dxp,nparm,nparm,1,dx)

return

end subroutine makedx_cvb
