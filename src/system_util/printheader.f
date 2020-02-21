************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
CSVC: print a banner with module name and runtime information
      subroutine print_module_header(modulename)
#ifdef _OPENMP
*$    use omp_lib
#endif
      implicit none
      character(*) :: modulename

      character(100) :: line
      integer :: i

#include "para_info.fh"
#include "unixinfo.fh"
#ifdef _MOLCAS_MPP_
      integer :: nprocs_global
      integer, external :: GAnNodes
      character(16) :: proc
#endif
#include "WrkSpc.fh"
      real*8 :: bytes
      integer :: order, group
      character(3) :: unit(0:8) =
     & ['  B',' kB',' MB',' GB',' TB', ' PB', ' EB', ' ZB', ' YB']
      character(16) :: memory
      integer :: nthreads
      character(16) :: threads

      write (6,*)
      write (6,'(50a)')('()',i=1,50)

      write(6,'(a)')
      line = '&' // trim(modulename)
      call upcase(line)
      call center_text(line)
      write(6,'(a)') trim(line)
      write(6,'(a)')

#ifdef _MOLCAS_MPP_
      nprocs_global = GAnNodes()
      write(proc,'(I16)') nprocs_global
      if (nprocs_global.gt.1) then
        if (nprocs.gt.1) then
          line = 'launched '//trim(adjustl(proc))//' MPI processes, '//
     &           'running in PARALLEL mode (work-sharing enabled)'
        else
          line = 'launched '//trim(adjustl(proc))//' MPI processes, '//
     &           'running in SERIAL mode (work-sharing disabled)'
        end if
      else
        line = 'only a single process is used, running in SERIAL mode'
      end if
#else
      line = 'only a single process is used'
#endif
      call center_text(line)
      write(6,'(a)') trim(line)

#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif

* mxmem is the number of 8-byte words (real*8) we have available.
      bytes = 8*mxmem
      order = FLOOR(LOG10(bytes))
      group = MIN(order/3,8)
      if (MOD(order,3).eq.0) then
        write (memory,'(F3.1,A)') bytes/10**(3*group),unit(group)
      else
        write (memory,'(I3,A)') INT(bytes/10**(3*group)),unit(group)
      end if
* report the maximum number of threads available
      if (nthreads.eq.1) then
        write (threads,'(A)') '1 thread'
      else
        write (threads,'(I8,A8)') nthreads, ' threads'
      end if
* if OPENMP is not compiled in, we don't know how many threads
* could be available in linear algebra libraries...
#ifndef _OPENMP
      threads=trim(threads)//'?'
#endif

      line = 'available to each process: '//
     &        trim(adjustl(memory))//' of memory, '//
     &        trim(adjustl(threads))
      call center_text(line)
      write(6,'(a)') trim(line)

      line = 'pid:'
#ifdef _MOLCAS_MPP_
      nprocs_global = GAnNodes()
      write(proc,'(I16)') nprocs_global
      if (nprocs_global.gt.1) then
        if (nprocs.gt.1) line = 'master pid:'
      end if
#endif
      write(line,'(a,1x,i0)') trim(line),pid
      call center_text(line)
      write(6,'(a)') trim(line)

      write (6,'(50a)')('()',i=1,50)
      write (6,*)
      end

      subroutine center_text(line)
CSVC: centers the text of a line
      implicit none
      character(*) :: line
      character(100) :: text
      integer :: linewidth, textwidth, textoffset

      text = adjustl(line)
      linewidth = len(line)
      textwidth = len_trim(text)
      textoffset = (linewidth-textwidth)/2
      line = ' '
      line(textoffset+1:textoffset+textwidth)=text
      end
