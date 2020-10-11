************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2015, Steven Vancoillie                                *
************************************************************************
* Molcas task manager
* Steven Vancoillie, june 2015
*   Reimplementation of the task manager routines
*   based solely on atomic read+increment operations on
*   a global task counter.
*
*   For GA, we use the ga_read_inc function, and for DGA
*   we use the gtsk_nxtval function. The latter has a special
*   "gtsk_nxtval_even" sibling for better performance. It is
*   only used when both MPI and DGA are active. The routines
*   that use this special flavour can be found in the second
*   part of the file.
*
*   The task lists are stacked, and have to be freed in
*   the order of last initialized = freed first!

      block data block_tsk
        implicit none
#include "tsk.fh"
        data list_counter/0/
      end

      subroutine init_tsk(id,n)
      implicit none
#include "tsk.fh"
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
#endif
      integer :: id, n

#ifdef _debug_trace_
#endif

      if (list_counter.eq.mxtsklst) then
        call sysabendmsg ('init_tsk',
     &    'no free task lists available',' ')
      end if
      list_counter = list_counter + 1

      id = list_counter
      ntasks(id) = n
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        if (.not.ga_create(MT_INT,1,1,'gltskl',0,0,task_counter(id)))
     &    call sysabendmsg ('init_tsk',
     &      'failed to create global task list',' ')
#  if !defined (_GA_)
        call ga_zero(task_counter(id))
        call gtsk_setup(id,task_counter(id))
#  else
        call ga_fill(task_counter(id),1)
#  endif
      else
        task_counter(id) = 1
      endif
#else
      task_counter(id) = 1
#endif

#ifdef _debug_trace_
#endif
      end

      subroutine free_tsk(id)
      implicit none
#include "tsk.fh"
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#endif
      integer :: id

#ifdef _debug_trace_
#endif

      if (list_counter.eq.0) then
        call sysabendmsg ('free_tsk',
     &    'attempting to free a non-existent task list.',' ')
      end if
      if (id .ne. list_counter) then
        call sysabendmsg ('free_tsk',
     &    'only stack-based task lists are supported.',' ')
      end if
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        if (.not.ga_destroy(task_counter(id)))
     &    call sysabendmsg ('free_tsk',
     &      'failed to destroy global task list.',' ')
#  if !defined (_GA_)
        call gtsk_reset(id)
#  endif
      end if
#endif
      list_counter = list_counter - 1

#ifdef _debug_trace_
#endif
      end

      logical function rsv_tsk(id,task)
      implicit none
#include "tsk.fh"
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#endif
      integer :: id, task

#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
* (atomically) read+increment next task number
#  if !defined (_GA_)
        task = gtsk_nxtval(id,1)
#  else
        task = ga_read_inc(task_counter(id),1,1,1)
#  endif
      else
        task = task_counter(id)
        task_counter(id) = task_counter(id) + 1
      endif
#else
      task = task_counter(id)
      task_counter(id) = task_counter(id) + 1
#endif
      rsv_tsk = task.le.ntasks(id)
      end

*****************************************************************
* The "even" flavour of the task routines just give each process
* a spread of numbers just as they would do with a loop stride.
*****************************************************************
      subroutine init_tsk_even(id,n)
      implicit none
#include "tsk.fh"
#include "para_info.fh"
      integer :: id, n

#ifdef _debug_trace_
#endif

      if (list_counter.eq.mxtsklst) then
        call sysabendmsg ('init_tsk_even',
     &    'no free task lists available',' ')
      end if
      list_counter = list_counter + 1

      id = list_counter
      ntasks(id) = n
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        task_counter(id) = myrank + 1
      else
        task_counter(id) = 1
      endif
#else
      task_counter(id) = 1
#endif

#ifdef _debug_trace_
#endif
      end

      subroutine free_tsk_even(id)
      implicit none
#include "tsk.fh"
      integer :: id

#ifdef _debug_trace_
#endif

      if (list_counter.eq.0) then
        call sysabendmsg ('free_tsk_even',
     &    'attempting to free a non-existent task list.',' ')
      end if
      if (id .ne. list_counter) then
        call sysabendmsg ('free_tsk_even',
     &    'only stack-based task lists are supported.',' ')
      end if
      list_counter = list_counter - 1

#ifdef _debug_trace_
#endif
      end

      logical function rsv_tsk_even(id,task)
      implicit none
#include "tsk.fh"
#include "para_info.fh"
      integer :: id, task

#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        task = task_counter(id)
* the next task is a stride of <nprocs> away
        task_counter(id) = task_counter(id) + nprocs
      else
        task = task_counter(id)
        task_counter(id) = task_counter(id) + 1
      endif
#else
      task = task_counter(id)
      task_counter(id) = task_counter(id) + 1
#endif
      rsv_tsk_even = task.le.ntasks(id)
      end
