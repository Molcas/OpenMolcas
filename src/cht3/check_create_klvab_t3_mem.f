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
        subroutine check_create_klvab_t3_mem (vblock)
!
! this routine finds the upper estimate of the memory
! requirements of the most demanding step in create_klvab_t3
!
        implicit none
!
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
        integer vblock,vblock_my
        integer mem,mem_trial,mem_avail
!
!.0 - calculate vblock_my
!
        call my_block (vblock,vblock_my)
!
        if (printkey.ge.10) then
        write (6,*)
        write (6,*) 'check_create_klvab_t3_mem '
        write (6,*)
        write (6,'(A,3(i5,1x))') 'nc,no,nv',nc,no,nv
        write (6,'(A,3(i5,1x))') 'maxdim,vblock,vblock_my',             &
     & maxdim,vblock,vblock_my
        end if
!
!.1 !create
        mem=vblock*vblock*(no+nv)+                                      &
     & nv*((nv*(nv+1))/2)+nv*nv+nc*maxdim+nc*maxdim*maxdim+             &
     & max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim)
!.2 !klvaa_vvvo
        mem_trial=vblock*vblock*(no+nv)+                                &
     & (nv*(nv*(nv+1))/2)+nv*nv+vblock_my*vblock_my*no*no+              &
     & 2*maxdim*maxdim*no*no
!
        if (mem_trial.gt.mem) mem=mem_trial
!.3 !create
        mem_trial=vblock*vblock*(no+nv)+                                &
     & (nv*(nv*(nv+1))/2)+nv*nv+vblock_my*vblock_my*no*no+              &
     & 2*maxdim*maxdim*no*no
!
        if (mem_trial.gt.mem) mem=mem_trial
!.4 !create
        mem_trial=no*no*vblock*(no+nv)+                                 &
     & no*nv*(no*(no+1)/2)+vblock*no*no+nc*(no*(no+1)/2)+               &
     & nc*no*nv+max(nc*((no*(no+1))/2),nc*no*maxdim,nc*no*nv)
!
        if (mem_trial.gt.mem) mem=mem_trial
!.5 !create
        mem_trial=no*no*vblock*(no+nv)+                                 &
     & no*nv*(no*(no+1)/2)+vblock*no*no+nv*vblock_my*no*no+             &
     & 2*maxdim*maxdim*no*no
!
        if (mem_trial.gt.mem) mem=mem_trial
!.6 !klvaa_oovo
        mem_trial=no*no*vblock*(no+nv)+                                 &
     & no*nv*(no*(no+1)/2)+vblock*no*no+                                &
     & nv*vblock_my*(((no-1)*no)/2)+2*maxdim*maxdim*no*no
!
        if (mem_trial.gt.mem) mem=mem_trial
!.7 !klvaa_oovo
        mem_trial=(((no-1)*no)/2)*vblock*vblock+                        &
     & vblock_my*vblock_my*no*no+nc*no*maxdim+                          &
     & 2*max(nc*no*maxdim,maxdim*maxdim*no*no)
!
        if (mem_trial.gt.mem) mem=mem_trial
!.8 !klvaa_oovo
        mem_trial=no*no*vblock*vblock+                                  &
     & vblock_my*vblock_my*no*no+nc*no*maxdim+                          &
     & 2*max(nc*no*maxdim,maxdim*maxdim*no*no)
!
        if (mem_trial.gt.mem) mem=mem_trial
!
           if (printkey.ge.10) then
              write (6,*)
              write (6,'(A,f10.1,A,f7.1,A,f3.1,A)')                     &
     &       'Memory required for the reorg. step = ',                  &
     &        (8.0d0*mem)/(1024),' kb ',                                &
     &        (8.0d0*mem)/(1024*1024),' Mb ',                           &
     &        (8.0d0*mem)/(1024*1024*1024),' Gb '
           end if
!
! - calculate available free memory
!
        Call GetMem('(T)','Max','Real',mem_avail,mem_avail)
!
           if (printkey.ge.10) then
              write (6,'(A,f10.1,A,f7.1,A,f3.1,A)')                     &
     &       'Available memory                    = ',                  &
     &        (8.0d0*mem_avail)/(1024),' kb ',                          &
     &        (8.0d0*mem_avail)/(1024*1024),' Mb ',                     &
     &        (8.0d0*mem_avail)/(1024*1024*1024),' Gb '
              write (6,*)
           end if
!
! - check, if mem fits
!
        if (mem_avail.lt.mem) then
          write (6,*) 'Not enough memory for the transformation step '
          call Abend()
        end if
!
        return
        end
