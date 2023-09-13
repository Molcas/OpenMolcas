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

subroutine rdcivec_cvb(detvec,fn,reord)
!***********************************************************************
!                                                                      *
!     Read the contents of the JOBIPH file.                            *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "rasdim.fh"
#include "jobiph_j.fh"
character*(*) fn
logical debug
data debug/.false./
dimension detvec(*)
dimension ncix(8)
logical reord
dimension rdum(1)

iwr = 0

call getnci_cvb(ncix,nactel_j,ispin_j-1,lsym_j)
ndet_j = ncix(1)

lujob = 15
call daname_cvb(lujob,fn)
! Allocate at least NDET words for each vector, since this is
! required by csdtvc:
!      Call GetMem('OCIvec','Allo','Real',ipCI,nConf_j*nroots_j)
call GetMem('OCIvec','Allo','Real',ipCI,nConf_j*nroots_j+ndet_j-nconf_j)
if (iwr == 0) then
  do i=1,nroots_j
    j = iroot_j(i)
    iDisk = iadr15_j(4)
    do k=1,j-1
      call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
    end do
    call dDaFile(LuJob,2,Work(ipCI+(i-1)*nconf_j),nConf_j,iDisk)
  end do

  if (reord) then
    call GetMem('ipci2','Allo','Real',ipCI2,nConf_j)
    call reord2_cvb(work(ipci),work(ipci2),1)
    call fmove_cvb(work(ipci2),work(ipci),nconf_j)
    call GetMem('ipci2','Free','Real',ipCI2,idum)
  end if

  call csf2det_cvb(work(ipci),detvec,lsym_j,1)
else if (iwr == 1) then
  call csf2det_cvb(work(ipci),detvec,lsym_j,2)

  if (reord) then
    call GetMem('ipci2','Allo','Real',ipCI2,nConf_j)
    call reord2_cvb(work(ipci),work(ipci2),0)
    call fmove_cvb(work(ipci2),work(ipci),nconf_j)
    call GetMem('ipci2','Free','Real',ipCI2,idum)
  end if

  do i=1,nroots_j
    j = iroot_j(i)
    iDisk = iadr15_j(4)
    do k=1,j-1
      call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
    end do
    call dDaFile(LuJob,1,Work(ipCI+(i-1)*nconf_j),nConf_j,iDisk)
  end do
end if
if (debug) then
  do i=0,nroots_j-1
    write(6,'(a,i3,a)') ' (CSF) CI vector ',i+1,' :'
    write(6,'(a)') ' ---------------------'
    call mxprint_cvb(work(ipci+nconf_j*i),1,nconf_j,0)
  end do
end if
call GetMem('OCIvec','Free','Real',ipCI,idum)
call daclos_cvb(lujob)

return

end subroutine rdcivec_cvb
