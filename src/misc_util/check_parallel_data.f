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
* Copyright (C) 2018, Denis Jelovina                                   *
************************************************************************


c if action="C" checks if array x(N) is identical across processes
c               returns status=.true. if data are idential
c               oterwise returns status=.false.
c                value of status is rank-independent
c if action="S" copy data from masterto all processes
#ifdef _MOLCAS_MPP_
      subroutine check_parallel_data(x,n,status,action)
      implicit none
      integer :: n
      real *8 :: x(n)
      logical :: status
      integer :: itemp,ns,j,i,k,irank
      character :: action
#include "para_info.fh"
#include "WrkSpc.fh"
      status=.true.
      if (nProcs.eq.1) return

      ns=n*nProcs
      call getmem('x_prll','ALLO','REAL',itemp, ns)
      do i=1,ns
         Work(itemp+i-1)=0
      end do
c
      j=MyRank*n
      do i=1,ns
         Work(j+itemp+i-1)=x(i)
      end do
      call GADsum(Work(itemp),ns)
c
      if(action.eq."C")then
         do irank=0,nProcs-1
c     comparing MyRank vs rank j
            if(irank.eq.MyRank) cycle
            i=MyRank*n
            j=irank*n
            do k=1,n
               if(Work(i+itemp+k-1).ne.Work(j+itemp+k-1))then
                  status=.false.
                  exit
               end if
            end do
         end do

      else if(action.eq."S")then
c     copy data from master to MyRank
         if(MyRank.ne.0)then
            j=0
            do i=1,ns
               x(i)=Work(j+itemp+i-1)
            end do
         end if
c
      else
         write (6,*) "check_parallel_data(), illegal value:"
         write (6,*) "action=",action
         write (6,*) "correct function call!!"
         call abort()
      end if

      call getmem('x_prll','FREE','REAL',itemp,ns)
      return
      end
#elif defined (NAGFOR)
      subroutine empty_check_parallel_data()
c     NAG compiler doesn't line empty files 
      implicit none
      return
      end
#endif
