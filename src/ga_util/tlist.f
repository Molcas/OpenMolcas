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
      Subroutine Init_TList(Triangular,P_Eff)
      implicit real*8 (a-h,o-z)
c      real*8  distrib,PQpTsk,TskLw,TskHi,MinPQ1,a,fint,tskmin,tskmax
      Logical Triangular,Alloc
#include "real.fh"
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
c      fint(a)=a-dmod(a,one)
c     Write (*,*) 'T_Status=',T_Status
      If (T_Status.eq.Active) return
      T_Status=Active
      Alloc=.false.
      Call IA_TList(Triangular,P_Eff, Alloc)
      return
      end
*
      Subroutine Alloc_TList(Triangular,P_Eff)
      implicit real*8 (a-h,o-z)
c      real*8  distrib,PQpTsk,TskLw,TskHi,MinPQ1,a,fint,tskmin,tskmax
      Logical Triangular,Alloc
#include "real.fh"
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
      If (T_Status.eq.Active) return
      Alloc=.true.
      Call IA_TList(Triangular,P_Eff, Alloc)
      return
      end
*
* removed _entry_
*
      Subroutine IA_TList(Triangular,P_Eff, Alloc)
      implicit real*8 (a-h,o-z)
      real*8  distrib,PQpTsk,TskLw,TskHi,MinPQ1,a,fint,tskmin,tskmax
      Logical Triangular,Alloc
* parameters concerning task distribution...
      Integer iDen_PQ, iDen_Tsk, MinPQ, MxnTsk
      Parameter ( iDen_PQ  = 2 )
      Parameter ( iDen_Tsk = 4 )
      Parameter ( MinPQ    = 4 )
* max number of tasks in tasklist per node...
      Parameter ( MxnTsk = 100 )
#include "real.fh"
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
      save ntasks_alloc
      fint(a)=dble(int(a))
      MinPQ1= MinPQ
      MxnTsk1=MxnTsk
      iDen_PQ1=iDen_PQ
      iDen_Tsk1=iDen_Tsk
c
C     Call Nr_Shells(nSkal)
C     P  = nSkal*(nSkal+1)/2
      P = P_Eff
      If (Triangular) Then
         PQ = P*(P+1d0)/2d0
      Else
         PQ = P*P
      End If
      if(Alloc) then
        ipTskQ = 0
        ntasks_alloc = 0
      end if
      nTasks = nint(Min(PQ,dble(MxnTsk1*nProcs)))
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
*
      if(Alloc) then
        ntasks_alloc=ntasks
        Call GetMem('TskMap','ALLO','Real',ipTskM,2*nTasks)
        call dcopy_(2*nTasks,Zero,0,Work(ipTskM),1)
        Call GetMem('TskQ','ALLO','Real',ipTskQ,2*nTasks)
c       Write (*,*) 'init_tlist ipTskQ @ ',ipTskQ,' nTasks=',nTasks
        call dcopy_(2*nTasks,Not_Used,0,Work(ipTskQ),1)
        Call GetMem('TskLst','ALLO','INTE',ipTskL,nTasks*2)
        Return
      end if
*
      tskmin=1.d14
      tskmax=zero
      TskLw=one
      TskHi=zero
      iTsk=0

*     REPEAT
  100 Continue
        distrib = fint(PQ/dble(iDen_PQ1))
        nTaskpP=nTasks/nProcs
        nTaskpP_seg=nTaskpP/iDen_Tsk1
        PQpTsk=MinPQ1
        If (nTaskpP_seg.ge.1)
     >    PQpTsk=fint(distrib/dble(nTaskpP_seg*nProcs))
        If (PQpTsk.gt.MinPQ1) Then
          PQpTsk  = fint(distrib/dble(nTaskpP_seg*nProcs))
          distrib = fint(PQpTsk*dble(nTaskpP_seg*nProcs))
          PQ=PQ-distrib
          kTskHi=nTaskpP_seg*nProcs
        Else If (PQ.gt.MinPQ1 .and. nTasks.gt.1) Then
          PQpTsk=Max(MinPQ1,fint((PQ+dble(nTasks)-2d0)/dble(nTasks-1)))
          kTskHi=nint(fint(PQ/PQpTsk))
          PQ=PQ-dble(kTskHi)*PQpTsk
        Else If (PQ.gt.0d0) Then
          PQpTsk=PQ
          kTskHi=1
          PQ=0d0
        Else
          kTskHi=0
          Write (6,*) 'Init_TList: you should not be here!'
          Call Abend()
        End If
        Do kTsk = 1, kTskHi
          TskHi=TskHi+PQpTsk
          Work(ipTskM+2*iTsk)=TskLw
          Work(ipTskM+2*iTsk+1)=TskHi
          tskmin=min(tskmin,(TskHi-TskLw+one))
          tskmax=max(tskmax,(TskHi-TskLw+one))
          TskLw=TskHi+one
          iTsk=iTsk+1
        End Do
        nTasks=nTasks-kTskHi
      If (abs(PQ).gt.1.d-10) Go To 100
*     UNTIL (PQ == 0)
      If (nTasks.lt.0) Then
          Write (6,*) 'nTasks.lt.0'
          Write (6,*) 'MyRank=',MyRank
          Call Abend
      End If
      nTasks=iTsk
*
c     Call RecPrt('TskM',' ',Work(ipTskM),2,nTasks)
      Return
      End
*
      Subroutine Free_TList
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
*
c     Write (*,*) 'Free_TList, T_Status=',T_Status
      If (T_Status.ne.Active) Return
      T_Status=Inactive
*
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
C     Call GetMem('TskQ  ','Free','Real',ipTskQ,2*nTasks_alloc)
C     Call GetMem('TskMap','Free','Real',ipTskM,2*nTasks_alloc)
      Call Free_Work(ipTskQ)
      Call Free_Work(ipTskM)
*
c     Write (*,*) 'Free_TList, T_Status=',T_Status
      Return
      End
