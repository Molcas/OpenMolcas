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

subroutine check_create_klvab_t3_mem(vblock)
! this routine finds the upper estimate of the memory
! requirements of the most demanding step in create_klvab_t3

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: vblock
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
integer(kind=iwp) :: mem, mem_avail, mem_trial, vblock_my
real(kind=wp), parameter :: kb = 1024.0_wp

!.0 - calculate vblock_my

call my_block(vblock,vblock_my)
!
if (printkey >= 10) then
  write(u6,*)
  write(u6,*) 'check_create_klvab_t3_mem '
  write(u6,*)
  write(u6,'(A,3(i5,1x))') 'nc,no,nv',nc,no,nv
  write(u6,'(A,3(i5,1x))') 'maxdim,vblock,vblock_my',maxdim,vblock,vblock_my
end if

!.1 !create
mem = vblock*vblock*(no+nv)+nv*((nv*(nv+1))/2)+nv*nv+nc*maxdim+nc*maxdim*maxdim+ &
      max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim)
!.2 !klvaa_vvvo
mem_trial = vblock*vblock*(no+nv)+(nv*(nv*(nv+1))/2)+nv*nv+vblock_my*vblock_my*no*no+2*maxdim*maxdim*no*no

if (mem_trial > mem) mem = mem_trial
!.3 !create
mem_trial = vblock*vblock*(no+nv)+(nv*(nv*(nv+1))/2)+nv*nv+vblock_my*vblock_my*no*no+2*maxdim*maxdim*no*no

if (mem_trial > mem) mem = mem_trial
!.4 !create
mem_trial = no*no*vblock*(no+nv)+no*nv*(no*(no+1)/2)+vblock*no*no+nc*(no*(no+1)/2)+nc*no*nv+ &
            max(nc*((no*(no+1))/2),nc*no*maxdim,nc*no*nv)

if (mem_trial > mem) mem = mem_trial
!.5 !create
mem_trial = no*no*vblock*(no+nv)+no*nv*(no*(no+1)/2)+vblock*no*no+nv*vblock_my*no*no+2*maxdim*maxdim*no*no

if (mem_trial > mem) mem = mem_trial
!.6 !klvaa_oovo
mem_trial = no*no*vblock*(no+nv)+no*nv*(no*(no+1)/2)+vblock*no*no+nv*vblock_my*(((no-1)*no)/2)+2*maxdim*maxdim*no*no

if (mem_trial > mem) mem = mem_trial
!.7 !klvaa_oovo
mem_trial = (((no-1)*no)/2)*vblock*vblock+vblock_my*vblock_my*no*no+nc*no*maxdim+2*max(nc*no*maxdim,maxdim*maxdim*no*no)

if (mem_trial > mem) mem = mem_trial
!.8 !klvaa_oovo
mem_trial = no*no*vblock*vblock+vblock_my*vblock_my*no*no+nc*no*maxdim+2*max(nc*no*maxdim,maxdim*maxdim*no*no)

if (mem_trial > mem) mem = mem_trial

if (printkey >= 10) then
  write(u6,*)
  write(u6,'(A,f10.1,A,f7.1,A,f3.1,A)') 'Memory required for the reorg. step = ',real(8*mem,kind=wp)/kb,' kb ', &
                                        real(8*mem,kind=wp)/kb**2,' Mb ',real(8*mem,kind=wp)/kb**3,' Gb'
end if

! - calculate available free memory

call GetMem('(T)','Max','Real',mem_avail,mem_avail)
!
if (printkey >= 10) then
  write(u6,'(A,f10.1,A,f7.1,A,f3.1,A)') 'Available memory                    = ',real(8*mem_avail,kind=wp)/kb,' kb ', &
                                        real(8*mem_avail,kind=wp)/kb**2,' Mb ',real(8*mem_avail,kind=wp)/kb**3,' Gb'
  write(u6,*)
end if

! - check, if mem fits

if (mem_avail < mem) then
  write(u6,*) 'Not enough memory for the transformation step '
  call Abend()
end if

return

end subroutine check_create_klvab_t3_mem
