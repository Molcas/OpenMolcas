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
* Copyright (C) 1995, Martin Schuetz                                   *
*               1998, Roland Lindh                                     *
*               2000-2015, Steven Vancoillie                           *
************************************************************************
************************************************************************
* This Module contains subroutines and functions which interface calls *
* to the Global Array Tools (GA)                                       *
*  DISTRIBUTED DATA PARALLEL VERSION for SCF                           *
************************************************************************
* SubRoutine GAInit                                                    *
*  ->     initialize GA, returns rank of process and # of processes    *
* SubRoutine GATerminate                                               *
*  ->     finalize GA                                                  *
* SubRoutine GASync                                                    *
*  -> synchronize processes...                                         *
* Subroutine GABrdcst(dType,Buf,nByte,Root)                            *
*  -> broadcast message Buf of size nByte from Root                    *
* SubRoutine GAStp(msg,ierr)                                           *
*  -> terminate parallel application...                                *
*     msg:      message to be printed...                               *
*     ierr:     error code...                                          *
* SubRoutine GAIGOP(k,n,op)                                            *
*  -> integer global operation; stub routine to ga_igop...             *
*     k(n):     global vector                                          *
*     op:       global operation '+','*','max','min','absmax','absmin' *
* SubRoutine GADGOP(x,n,op)                                            *
*  -> double global operation; stub routine to ga_dgop...              *
*     x(n):     global vector                                          *
*     op:       global operation '+','*','max','min','absmax','absmin' *
* SubRoutine GAAccP(iGA,ilo,ihi,jlo,jhi,buf,ld,alpha)                  *
*  -> accumulate to GA patch; stub routine to ga_acc...                *
*     iGA:      GA handle                                              *
*     ilo,ihi,jlo,jhi: defines GA patch...                             *
*     buf:      local buffer, containing data to accumulate...         *
*     ld:       leading dimension of buf...                            *
*     alpha:    scaling factor...                                      *
* SubRoutine GADupl(iGA1,iGA2)                                         *
*  -> duplicate & copy a global array...                               *
*     iGA1,iGA2  :       GA handles...                                 *
* SubRoutine GAAdd(alpha,iGA1,beta,iGA2,iGA3)                          *
*  -> add to global arrays; stub routine to ga_dadd...                 *
*     iGA3 = alpha * iGA1 + beta * iGA2                                *
*----------------------------------------------------------------------*
*     written by:                                                      *
*     M. Schuetz, University of Lund, Sweden, 1995                     *
*                                                                      *
*     modified by:                                                     *
*     R. Lindh, University of Lund, Sweden, 1998                       *
*     S. Vancoillie, University of Lund, Sweden, 2010-2015             *
************************************************************************
      SubRoutine GAInit
*     purpose: initialize DGA and set the global rank and number of    *
*              processes in mpp_procid and mpp_nprocs. Then also set   *
*              the (initial) local myRank and nProcs varibles.         *
*     called from: DPMP2 (distributed parallel MP2)                    *
*     calls to: MPI-2/DGA routines                                     *
      Implicit None
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#  include "mpp_info.fh"
#  include "global.fh"
      Character(8) :: molcas_nprocs_env
      Integer :: molcas_nprocs, iRC

C     SVC: bypass MPI initialization if only 1 process, this is needed for a
C     specific version of GEO (so that the serial tasks which are run by MPI
C     do not try to re-initialize MPI). This has the consequence that for any
C     calculation where the number of processes is 1, calls to MPI/GA will fail
C     at runtime (even though it will compile when inside _MOLCAS_MPP_!)
      Call getenvf('MOLCAS_NPROCS',molcas_nprocs_env)
      If(molcas_nprocs_env(1:1).eq.' ') Then
         molcas_nprocs=-1
      Else
         Read(molcas_nprocs_env,*) molcas_nprocs
      End If
      If(molcas_nprocs.ne.1) Then
#  ifdef _GA_
         Call mpi_init(iRC)
         Call ga_initialize()
         Call ga_replace_ma()
#  else
         Call ga_initialize()
#  endif
         mpp_procid=ga_nodeid()
         mpp_nprocs=ga_nnodes()
         mpp_workshare = .True.
C make each slave process go to its proper work directory
         Call slaveschdir (mpp_procid, iRC)
         IF (iRC.NE.0) CALL Abend()
      Else
         mpp_procid=0
         mpp_nprocs=1
         mpp_workshare = .False.
      End If
      MyRank = mpp_procid
      nProcs = mpp_nprocs
#else
      MyRank=0
      nProcs=1
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GATerminate
CSVC: This terminates the parallel runtime after which the processes
C     are no longer allowed to make ga/mpi calls. When called more than
C     once, the routine does nothing. This is to support early
C     termination of the parallel runtime without actually exiting.
C     In such a situation, when the process eventually finishes, it
C     will call this routine again, but then doing nothing. Such a use
C     case is e.g. when we want to terminate slave processes and only
C     continue to run the master process in serial mode.
      Implicit None
#ifdef _MOLCAS_MPP_
#  include "mpp_info.fh"
#  include "global.fh"
      Logical, Save :: FirstCall = .true.
#ifdef _GA_
      Integer iErr
#endif
      Logical, External :: Is_Real_Par

      if (FirstCall) then
        FirstCall=.false.
        if(mpp_nprocs.gt.1) then
          Call ga_terminate()
#  ifdef _GA_
          Call mpi_finalize(iErr)
#  endif
        endif
      endif
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GASync
      Implicit None
#ifdef _MOLCAS_MPP_
#  include "global.fh"
      Logical, External :: Is_Real_Par

      If (Is_Real_Par()) Then
         Call ga_sync()
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      Subroutine GABrdcst(dType,Buf,nByte,Root)
      Implicit None
      Integer       dType,nByte,Root
      Character(*)  buf
#ifdef _MOLCAS_MPP_
#  include "global.fh"

      Logical, External :: Is_Real_Par
      If (Is_Real_Par()) CALL GA_Brdcst(dType,Buf,nByte,Root)
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(dType)
         Call Unused_character(Buf)
         Call Unused_integer(nByte)
         Call Unused_integer(Root)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GAStp(msg,ierr)
      Implicit None
      Character*(*) msg
      Integer ierr
#ifdef _MOLCAS_MPP_
#  include "global.fh"

      Call ga_error(msg,ierr)
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_character(msg)
         Call Unused_integer(ierr)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GADGOP(x,n,op)
      Implicit None
      Integer n
      Real*8 x(n)
      Character*(*) op
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
      Logical, External :: Is_Real_Par

      If (Is_Real_Par()) Then
         Call ga_dgop(MT_DBL,x,n,op)
      End If
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(x)
         Call Unused_character(op)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GAdGOp_Scal(x,op)
      Implicit None
      Real*8 x
      Character(*) op
      Real*8 x_arr(1)

      x_arr(1)=x
      Call GAdGOp(x_arr,1,op)
      x=x_arr(1)
      End
*----------------------------------------------------------------------*
      SubRoutine GADSUM(x,n)
      Implicit None
      Integer n
      Real*8 x(n)
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
      Logical, External ::  Is_Real_Par

      If (Is_Real_Par()) Then
         Call ga_dgop(MT_DBL,x,n,'+')
      End If
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(x)
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GAdSum_Scal(x)
      Implicit None
      Real*8 x
      Real*8 x_arr(1)
      x_arr(1)=x
      Call GAdSum(x_arr,1)
      x=x_arr(1)
      End
*----------------------------------------------------------------------*
      SubRoutine GAIGOP(k,n,op)
      Implicit None
      Integer n
      Integer k(n)
      Character*(*) op
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
      Logical, External :: Is_Real_Par

      If (Is_Real_Par()) Then
         Call ga_igop(MT_INT,k,n,op)
      End If
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(k)
         Call Unused_character(op)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GAiGOp_Scal(k,op)
      Implicit None
      Integer k
      Character(*) op
      Integer k_arr(1)

      k_arr(1)=k
      Call GAiGOp(k_arr,1,op)
      k=k_arr(1)
      End
*----------------------------------------------------------------------*
      SubRoutine GAAccP(iGA,ilo,ihi,jlo,jhi,buf,ld,alpha)
      Implicit None
      Integer iGA,ilo,ihi,jlo,jhi,ld
      Real*8 buf(1:ld,1:*),alpha
#ifdef _MOLCAS_MPP_
#  include "global.fh"
      Logical, External :: Is_Real_Par

      If (Is_Real_Par()) Then
         Call ga_acc(iGA,ilo,ihi,jlo,jhi,buf,ld,alpha)
      End IF
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iGA)
         Call Unused_integer(ilo)
         Call Unused_integer(ihi)
         Call Unused_integer(jlo)
         Call Unused_integer(jhi)
         Call Unused_real_array(buf)
         Call Unused_real(alpha)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GADupl(iGA1,iGA2)
      Implicit None
      Integer iGA1,iGA2
#ifdef _MOLCAS_MPP_
#  include "global.fh"
      Logical, External :: Is_Real_Par
      Character gaLbl*5,gaLbl2*6
      Logical ok

      If (.Not. Is_Real_Par()) Return
      If (iGA1.ge.0) Return
      Call ga_inquire_name(iGA1,gaLbl)
      Write(gaLbl2,'(A,I1)') gaLbl,2

      ok=ga_duplicate(iGA1,iGA2,gaLbl2)
      If (.NOT.ok) Then
        Write (6,*) 'GADupl: ga_duplicate not OK!'
        Call GAStp('GADupl',42)
        Call QTrace
        Call Abend()
      End If

      Call ga_copy(iGA1,iGA2)
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iGA1)
         Call Unused_integer(iGA2)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      SubRoutine GAAdd(alpha,iGA1,beta,iGA2,iGA3)
      Implicit None
      Integer iGA1,iGA2,iGA3
      Real*8 alpha,beta
#ifdef _MOLCAS_MPP_
#  include "global.fh"
      Logical, External :: Is_Real_Par

      If (.Not. Is_Real_Par()) Return
      If ((iGA1.ge.0).OR.(iGA2.ge.0).OR.(iGA3.ge.0)) Return
      Call ga_dadd(alpha,iGA1,beta,iGA2,iGA3)
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(alpha)
         Call Unused_integer(iGA1)
         Call Unused_real(beta)
         Call Unused_integer(iGA2)
         Call Unused_integer(iGA3)
      End If
#endif
      Return
      End
*----------------------------------------------------------------------*
      Integer Function GANodeID()
      Implicit None
#ifdef _MOLCAS_MPP_
#  include "mpp_info.fh"
      GANodeID = mpp_procid
#else
      GANodeID = 0
#endif
      Return
      End
*----------------------------------------------------------------------*
      Integer Function GAnNodes()
      Implicit None
#ifdef _MOLCAS_MPP_
#  include "mpp_info.fh"
      GAnNodes = mpp_nprocs
#else
      GAnNodes = 1
#endif
      Return
      End
