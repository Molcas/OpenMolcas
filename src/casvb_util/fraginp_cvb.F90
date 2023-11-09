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

subroutine fraginp_cvb()

use casvb_global, only: confsinp, i2s_fr, nalf_fr, nbet_fr, nconf, nconf_fr, nel_fr, nfrag, nMs_fr, noe, nS_fr
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: istr, istr2, mavaili, mxconf, nread
real(kind=wp) :: dum(1), Scurr
integer(kind=iwp), allocatable :: tmp(:,:)
integer(kind=iwp), parameter :: ncmp = 4, nstrin = 2
character(len=*), parameter :: string(nstrin) = ['WAVE    ','CON     ']

do
  call fstring_cvb(string,nstrin,istr,ncmp,2)
  if (istr == 1) then
    ! 'WAVE'
    nfrag = nfrag+1
    nel_fr(nfrag) = 0
    call int_cvb(nel_fr(nfrag),1,nread,1)
    nMs_fr(nfrag) = 0
    nS_fr(nfrag) = 0
    nalf_fr(1,nfrag) = 0
    nbet_fr(1,nfrag) = 0
    i2s_fr(1,nfrag) = -1
    do
      Scurr = -One
      call real_cvb(dum,1,nread,1)
      Scurr = dum(1)
      if (Scurr == -One) exit
      nS_fr(nfrag) = nS_fr(nfrag)+1
      i2s_fr(nS_fr(nfrag),nfrag) = nint(Two*Scurr)
    end do
  else if (istr == 2) then
    ! 'CON'
    if (nfrag == 0) then
      nfrag = 1
      nel_fr(nfrag) = 0
      nMs_fr(nfrag) = 0
      nS_fr(nfrag) = 0
      nalf_fr(1,nfrag) = 0
      nbet_fr(1,nfrag) = 0
      i2s_fr(1,nfrag) = -1
    end if

    call mma_maxINT(mavaili)
    mxconf = max(mavaili/2,0)/noe
    call mma_allocate(tmp,noe,mxconf,label='confsinp')
    if (allocated(confsinp)) then
      tmp(:,1:size(confsinp,2)) = confsinp(:,:)
      call mma_deallocate(confsinp)
    end if
    call move_alloc(tmp,confsinp)
    nconf_fr(nfrag) = 1
    nconf = nconf+1
    do
      if (mxconf < nconf) then
        write(u6,*) ' Insufficient memory for configuration read',mavaili,mxconf,nconf
        call abend_cvb()
      end if
      confsinp(:,nconf) = 0
      call int_cvb(confsinp(:,nconf),noe,nread,1)
      call fstring_cvb('CON',1,istr2,3,2)
      if (istr2 == 0) exit
      nconf_fr(nfrag) = nconf_fr(nfrag)+1
      nconf = nconf+1
    end do
    call mma_allocate(tmp,noe,nconf,label='confsinp')
    tmp(:,:) = confsinp(:,1:nconf)
    call mma_deallocate(confsinp)
    call move_alloc(tmp,confsinp)
  end if
  if (istr == 0) exit
end do

return

end subroutine fraginp_cvb
