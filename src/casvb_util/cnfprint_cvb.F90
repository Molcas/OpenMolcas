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

use casvb_global, only: ipr, nconf, nconf_fr, ndetvb_fr, nel, nel_fr, nfrag, noe, norb, nvbr_fr, recinp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: idum(1), ifrag, ioffs, nconf_off
integer(kind=iwp), allocatable :: scr(:)
logical(kind=iwp), external :: recinpcmp_cvb, up2date_cvb ! ... Make: up to date? ...

if (recinpcmp_cvb(4)) call touch_cvb('CNFPRINT')

if ((ipr(1) >= 0) .and. (.not. up2date_cvb('CNFPRINT'))) then
  call mma_allocate(scr,max(noe,noe*nconf),label='scr')
  call rdioff_cvb(1,recinp,ioffs)
  call rdis_cvb(idum,1,recinp,ioffs)
  !noe1 = idum(1)
  call rdis_cvb(idum,1,recinp,ioffs)
  !nconf1 = idum(1)
  call rdis_cvb(idum,1,recinp,ioffs)
  !kbasiscvb1 = idum(1)
  call rdis_cvb(scr,noe*nconf,recinp,ioffs)
  if (nconf == 0) then
    scr(1:min(nel,norb)) = 1
    scr(1:nel-norb) = 2
  end if
  nconf_off = 0
  do ifrag=1,nfrag
    if (nfrag > 1) write(u6,'(/,a,i3)') ' Configuration list for wavefunction fragment',ifrag
    write(u6,'(/,a)') ' Spatial VB configurations'
    write(u6,'(a)') ' -------------------------'
    write(u6,'(a)') '     Conf. =>   Orbitals'
    call cnfprt_cvb(scr(noe*nconf_off+1),nconf_fr(ifrag),nel_fr(ifrag))
    write(u6,'(/,a,i6)') ' Number of VB configurations :',nconf_fr(ifrag)
    write(u6,'(a,i6)') '           VB structures     :',nvbr_fr(ifrag)
    write(u6,'(a,i6)') '           VB determinants   :',ndetvb_fr(ifrag)
    nconf_off = nconf_off+nconf_fr(ifrag)
  end do
  call mma_deallocate(scr)
  call make_cvb('CNFPRINT')
end if

return

end subroutine cnfprint_cvb
