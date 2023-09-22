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

subroutine cnfprint_cvb()

use casvb_global, only: nconf_fr, ndetvb_fr, nel_fr, nfrag, nvbr_fr
use Definitions, only: iwp, u6

implicit none
#include "main_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1, idum(1), ifrag, ioffs, nconf_off
integer(kind=iwp), external :: mstacki_cvb
logical(kind=iwp), external :: recinpcmp_cvb, up2date_cvb ! ... Make: up to date? ...

if (recinpcmp_cvb(4)) call touch_cvb('CNFPRINT')

if ((ip(1) >= 0) .and. (.not. up2date_cvb('CNFPRINT'))) then
  i1 = mstacki_cvb(max(noe,noe*nconf))
  call rdioff_cvb(1,recinp,ioffs)
  call rdis_cvb(idum,1,recinp,ioffs)
  !noe1 = idum(1)
  call rdis_cvb(idum,1,recinp,ioffs)
  !nconf1 = idum(1)
  call rdis_cvb(idum,1,recinp,ioffs)
  !kbasiscvb1 = idum(1)
  call rdis_cvb(iwork(i1),noe*nconf,recinp,ioffs)
  if (nconf == 0) then
    do i=1,min(nel,norb)
      iwork(i+i1-1) = 1
    end do
    do i=1,nel-norb
      iwork(i+i1-1) = 2
    end do
  end if
  nconf_off = 0
  do ifrag=1,nfrag
    if (nfrag > 1) write(u6,'(/,a,i3)') ' Configuration list for wavefunction fragment',ifrag
    write(u6,'(/,a)') ' Spatial VB configurations'
    write(u6,'(a)') ' -------------------------'
    write(u6,'(a)') '     Conf. =>   Orbitals'
    call cnfprt_cvb(iwork(noe*nconf_off+i1),nconf_fr(ifrag),nel_fr(ifrag))
    write(u6,'(/,a,i6)') ' Number of VB configurations :',nconf_fr(ifrag)
    write(u6,'(a,i6)') '           VB structures     :',nvbr_fr(ifrag)
    write(u6,'(a,i6)') '           VB determinants   :',ndetvb_fr(ifrag)
    nconf_off = nconf_off+nconf_fr(ifrag)
  end do
  call mfreei_cvb(i1)
  call make_cvb('CNFPRINT')
end if

return

end subroutine cnfprint_cvb
