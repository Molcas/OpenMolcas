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

subroutine setmocom_cvb()

use casvb_global, only: iact_mo, nact_mo, nbas_mo, nbasf_mo, nbasi_mo, nbasisq_mo, nbassqf_mo, nbassqi_mo, nsym_mo
use Definitions, only: iwp, u6

implicit none
#include "rasdim.fh"
#include "jobiph_j.fh"
integer(kind=iwp) :: i, ii, j
logical(kind=iwp), parameter :: debug = .false.

nsym_mo = nsym_j
nbasi_mo(:) = nbas_j(:)

nbas_mo = 0
nbasisq_mo = 0
do i=1,8
  nbassqi_mo(i) = nbasi_mo(i)*nbasi_mo(i)
  nbasf_mo(i) = nbas_mo
  nbassqf_mo(i) = nbasisq_mo
  nbas_mo = nbas_mo+nbasi_mo(i)
  nbasisq_mo = nbasisq_mo+nbassqi_mo(i)
end do

nact_mo = 0
do i=1,8
  do j=1,nash_j(i)
    nact_mo = nact_mo+1
    iact_mo(nact_mo) = nbasf_mo(i)+nfro_j(i)+nish_j(i)+j
  end do
end do

if (debug) then
  write(u6,*) ' MO interface'
  write(u6,*) ' ------------'
  write(u6,*) ' nsym    :',nsym_mo
  write(u6,*) ' nbas    :',nbas_mo
  write(u6,*) ' nbasisq :',nbasisq_mo
  write(u6,*) ' nbasi   :',nbasi_mo
  write(u6,*) ' nbassqi :',nbassqi_mo
  write(u6,*) ' nbasf   :',nbasf_mo
  write(u6,*) ' nbassqf :',nbassqf_mo
  write(u6,*) ' nact    :',nact_mo
  write(u6,*) ' iact    :',(iact_mo(ii),ii=1,nact_mo)
end if

return

end subroutine setmocom_cvb
