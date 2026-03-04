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
! Copyright (C) 1995, Martin Schuetz                                   *
!               1998, Roland Lindh                                     *
!               2000-2015, Steven Vancoillie                           *
!***********************************************************************
!***********************************************************************
! This Module contains subroutines and functions which interface calls *
! to the Global Array Tools (GA)                                       *
!  DISTRIBUTED DATA PARALLEL VERSION for SCF                           *
!***********************************************************************
! SubRoutine GAInit                                                    *
!  ->     initialize GA, returns rank of process and # of processes    *
! SubRoutine GATerminate                                               *
!  ->     finalize GA                                                  *
! SubRoutine GASync                                                    *
!  -> synchronize processes...                                         *
! Subroutine GABrdcst(dType,Buf,nByte,Root)                            *
!  -> broadcast message Buf of size nByte from Root                    *
! SubRoutine GAStp(msg,ierr)                                           *
!  -> terminate parallel application...                                *
!     msg:      message to be printed...                               *
!     ierr:     error code...                                          *
! SubRoutine GAIGOP(k,n,op)                                            *
!  -> integer global operation; stub routine to ga_igop...             *
!     k(n):     global vector                                          *
!     op:       global operation '+','*','max','min','absmax','absmin' *
! SubRoutine GADGOP(x,n,op)                                            *
!  -> double global operation; stub routine to ga_dgop...              *
!     x(n):     global vector                                          *
!     op:       global operation '+','*','max','min','absmax','absmin' *
! SubRoutine GAAccP(iGA,ilo,ihi,jlo,jhi,buf,ld,alpha)                  *
!  -> accumulate to GA patch; stub routine to ga_acc...                *
!     iGA:      GA handle                                              *
!     ilo,ihi,jlo,jhi: defines GA patch...                             *
!     buf:      local buffer, containing data to accumulate...         *
!     ld:       leading dimension of buf...                            *
!     alpha:    scaling factor...                                      *
! SubRoutine GADupl(iGA1,iGA2)                                         *
!  -> duplicate & copy a global array...                               *
!     iGA1,iGA2  :       GA handles...                                 *
! SubRoutine GAAdd(alpha,iGA1,beta,iGA2,iGA3)                          *
!  -> add to global arrays; stub routine to ga_dadd...                 *
!     iGA3 = alpha * iGA1 + beta * iGA2                                *
!----------------------------------------------------------------------*
!     written by:                                                      *
!     M. Schuetz, University of Lund, Sweden, 1995                     *
!                                                                      *
!     modified by:                                                     *
!     R. Lindh, University of Lund, Sweden, 1998                       *
!     S. Vancoillie, University of Lund, Sweden, 2010-2015             *
!***********************************************************************

subroutine GAInit()
!***********************************************************************
!     purpose: initialize DGA and set the global rank and number of    *
!              processes in mpp_procid and mpp_nprocs. Then also set   *
!              the (initial) local myRank and nProcs variables.        *
!     called from: DPMP2 (distributed parallel MP2)                    *
!     calls to: MPI-2/DGA routines                                     *
!***********************************************************************

use Para_Info, only: MyRank, nProcs
#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_nprocs, mpp_procid, mpp_workshare
use Definitions, only: iwp
#endif

implicit none
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: molcas_nprocs, iRC
character(len=8) :: molcas_nprocs_env
#include "global.fh"

! SVC: bypass MPI initialization if only 1 process, this is needed for a
! specific version of GEO (so that the serial tasks which are run by MPI
! do not try to re-initialize MPI). This has the consequence that for any
! calculation where the number of processes is 1, calls to MPI/GA will fail
! at runtime (even though it will compile when inside _MOLCAS_MPP_!)
call getenvf('MOLCAS_NPROCS',molcas_nprocs_env)
if (molcas_nprocs_env(1:1) == ' ') then
  molcas_nprocs = -1
else
  read(molcas_nprocs_env,*) molcas_nprocs
end if
if (molcas_nprocs /= 1) then
# ifdef _GA_
  call mpi_init(iRC)
  call ga_initialize()
  call ga_replace_ma()
# else
  call ga_initialize()
# endif
  mpp_procid = ga_nodeid()
  mpp_nprocs = ga_nnodes()
  mpp_workshare = .true.
  ! make each slave process go to its proper work directory
  call slaveschdir(mpp_procid,iRC)
  if (iRC /= 0) call Abend()
else
  mpp_procid = 0
  mpp_nprocs = 1
  mpp_workshare = .false.
end if
MyRank = mpp_procid
nProcs = mpp_nprocs
#else
MyRank = 0
nProcs = 1
#endif

end subroutine GAInit
!=!=
subroutine GATerminate()
! SVC: This terminates the parallel runtime after which the processes
! are no longer allowed to make ga/mpi calls. When called more than
! once, the routine does nothing. This is to support early
! termination of the parallel runtime without actually exiting.
! In such a situation, when the process eventually finishes, it
! will call this routine again, but then doing nothing. Such a use
! case is e.g. when we want to terminate slave processes and only
! continue to run the master process in serial mode.

#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_nprocs
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: FirstCall = .true.
#ifdef _GA_
integer(kind=iwp) :: iErr
#endif

if (FirstCall) then
  FirstCall = .false.
  if (mpp_nprocs > 1) then
    call ga_terminate()
#   ifdef _GA_
    call mpi_finalize(iErr)
#   endif
  end if
end if
#endif

end subroutine GATerminate
!=!=
subroutine GASync()

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par

implicit none

if (Is_Real_Par()) call ga_sync()
#endif

end subroutine GASync
!=!=
subroutine GABrdcst(dType,Buf,nByte,Root)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: dType, nByte, Root
character(len=*) :: Buf
#ifdef _MOLCAS_MPP_
interface
  subroutine GA_Brdcst(tp,buf,lenbuf,root)
    import :: iwp
    integer(kind=iwp) :: tp, lenbuf, root
    type(*) :: buf
  end subroutine GA_Brdcst
end interface

if (Is_Real_Par()) call GA_Brdcst(dType,Buf,nByte,Root)
#else
#include "macros.fh"
unused_var(dType)
unused_var(Buf)
unused_var(nByte)
unused_var(Root)
#endif

end subroutine GABrdcst
!=!=
subroutine GAStp(msg,ierr)

use Definitions, only: iwp

implicit none
character(len=*) :: msg
integer(kind=iwp) :: ierr
#ifdef _MOLCAS_MPP_

call ga_error(msg,ierr)
#else
#include "macros.fh"
unused_var(msg)
unused_var(ierr)
#endif

end subroutine GAStp
!=!=
subroutine GADGOP(x,n,op)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: x(n)
character(len=*) :: op
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"

if (Is_Real_Par()) call ga_dgop(MT_DBL,x,n,op)
#else
#include "macros.fh"
unused_var(x)
unused_var(op)
#endif

end subroutine GADGOP
!=!=
subroutine GAdGOp_Scal(x,op)

use Definitions, only: wp

implicit none
real(kind=wp) :: x
character(len=*) :: op
real(kind=wp) :: x_arr(1)

x_arr(1) = x
call GAdGOp(x_arr,1,op)
x = x_arr(1)

end subroutine GAdGOp_Scal
!=!=
subroutine GADSUM(x,n)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: x(n)
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"

if (Is_Real_Par()) call ga_dgop(MT_DBL,x,n,'+')
#else
#include "macros.fh"
unused_var(x)
#endif

end subroutine GADSUM
!=!=
subroutine GAdSum_Scal(x)

use Definitions, only: wp

implicit none
real(kind=wp) :: x
real(kind=wp) :: x_arr(1)

x_arr(1) = x
call GAdSum(x_arr,1)
x = x_arr(1)

end subroutine GAdSum_Scal
!=!=
subroutine GAIGOP(k,n,op)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: n, k(n)
character(len=*) :: op
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"

if (Is_Real_Par()) call ga_igop(MT_INT,k,n,op)
#else
#include "macros.fh"
unused_var(k)
unused_var(op)
#endif

end subroutine GAIGOP
!=!=
subroutine GAiGOp_Scal(k,op)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: k
character(len=*) :: op
integer(kind=iwp) :: k_arr(1)

k_arr(1) = k
call GAiGOp(k_arr,1,op)
k = k_arr(1)

end subroutine GAiGOp_Scal
!=!=
subroutine GAAccP(iGA,ilo,ihi,jlo,jhi,buf,ld,alpha)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iGA, ilo, ihi, jlo, jhi, ld
real(kind=wp) :: buf(ld,*), alpha
#ifdef _MOLCAS_MPP_

if (Is_Real_Par()) call ga_acc(iGA,ilo,ihi,jlo,jhi,buf,ld,alpha)
#else
#include "macros.fh"
unused_var(iGa)
unused_var(ilo)
unused_var(ihi)
unused_var(jlo)
unused_var(jhi)
unused_var(buf(1,1))
unused_var(alpha)
#endif

end subroutine GAAccP
!=!=
subroutine GADupl(iGA1,iGA2)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use Definitions, only: u6
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iGA1, iGA2
#ifdef _MOLCAS_MPP_
logical(kind=iwp) :: ok
character(len=6) :: gaLbl2
character(len=5) :: gaLbl
#include "global.fh"

if (.not. Is_Real_Par()) return
if (iGA1 >= 0) return
call ga_inquire_name(iGA1,gaLbl)
write(gaLbl2,'(A,I1)') gaLbl,2

ok = ga_duplicate(iGA1,iGA2,gaLbl2)
if (.not. ok) then
  write(u6,*) 'GADupl: ga_duplicate not OK!'
  call GAStp('GADupl',42)
  call Abend()
end if

call ga_copy(iGA1,iGA2)
#else
#include "macros.fh"
unused_var(iGA1)
unused_var(iGA2)
#endif

end subroutine GADupl
!=!=
subroutine GAAdd(alpha,iGA1,beta,iGA2,iGA3)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: alpha, beta
integer(kind=iwp) :: iGA1, iGA2, iGA3
#ifdef _MOLCAS_MPP_

if (.not. Is_Real_Par()) return
if ((iGA1 >= 0) .or. (iGA2 >= 0) .or. (iGA3 >= 0)) return
call ga_dadd(alpha,iGA1,beta,iGA2,iGA3)
#else
#include "macros.fh"
unused_var(alpha)
unused_var(iGA1)
unused_var(beta)
unused_var(iGA2)
unused_var(iGA3)
#endif

end subroutine GAAdd
!=!=
function GANodeID()

#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_procid
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: GANodeID

#ifdef _MOLCAS_MPP_
GANodeID = mpp_procid
#else
GANodeID = 0
#endif

end function GANodeID
!=!=
function GAnNodes()

#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_nprocs
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: GAnNodes

#ifdef _MOLCAS_MPP_
GAnNodes = mpp_nprocs
#else
GAnNodes = 1
#endif

end function GAnNodes
