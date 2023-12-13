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

subroutine print_module_header(modulename)
!SVC: print a banner with module name and runtime information

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs
#endif
#ifdef _OPENMP
use omp_lib, only: omp_get_max_threads
#endif
use UnixInfo, only: pid
use stdalloc, only: mxMem
use Definitions, only: wp, iwp, u6, RtoB

implicit none
character(len=*) :: modulename
character(len=100) :: line
integer(kind=iwp) :: order, group, nthreads
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: nprocs_global
character(len=16) :: proc
integer(kind=iwp), external :: GAnNodes
#endif
real(kind=wp) :: bytes
character(len=16) :: memory, threads
character(len=*), parameter :: unit(0:8) = ['  B',' kB',' MB',' GB',' TB',' PB',' EB',' ZB',' YB']
logical(kind=iwp), external :: Reduce_Prt

if (Reduce_Prt()) return
write(u6,*)
write(u6,'(a)') repeat('()',50)

write(u6,'(a)')
line = '&'//trim(modulename)
call upcase(line)
call center_text(line)
write(u6,'(a)') trim(line)
write(u6,'(a)')

#ifdef _MOLCAS_MPP_
nprocs_global = GAnNodes()
write(proc,'(I16)') nprocs_global
if (nprocs_global > 1) then
  if (nprocs > 1) then
    line = 'launched '//trim(adjustl(proc))//' MPI processes, running in PARALLEL mode (work-sharing enabled)'
  else
    line = 'launched '//trim(adjustl(proc))//' MPI processes, running in SERIAL mode (work-sharing disabled)'
  end if
else
  line = 'only a single process is used, running in SERIAL mode'
end if
#else
line = 'only a single process is used'
#endif
call center_text(line)
write(u6,'(a)') trim(line)

#ifdef _OPENMP
nthreads = omp_get_max_threads()
#else
nthreads = 1
#endif

! mxmem is the number of 8-byte words (real*8) we have available.
bytes = RtoB*mxMem
order = floor(log10(bytes))
group = min(order/3,8)
if (mod(order,3) == 0) then
  write(memory,'(F3.1,A)') bytes/10**(3*group),unit(group)
else
  write(memory,'(I3,A)') int(bytes/10**(3*group)),unit(group)
end if
! report the maximum number of threads available
if (nthreads == 1) then
  write(threads,'(A)') '1 thread'
else
  write(threads,'(I8,A8)') nthreads,' threads'
end if
! if OPENMP is not compiled in, we don't know how many threads
! could be available in linear algebra libraries...
#ifndef _OPENMP
threads = trim(threads)//'?'
#endif

line = 'available to each process: '//trim(adjustl(memory))//' of memory, '//trim(adjustl(threads))
call center_text(line)
write(u6,'(a)') trim(line)

line = 'pid:'
#ifdef _MOLCAS_MPP_
nprocs_global = GAnNodes()
write(proc,'(I16)') nprocs_global
if (nprocs_global > 1) then
  if (nprocs > 1) line = 'master pid:'
end if
#endif
write(line,'(a,1x,i0)') trim(line),pid
call center_text(line)
write(u6,'(a)') trim(line)

write(u6,'(a)') repeat('()',50)
write(u6,*)

end subroutine print_module_header
