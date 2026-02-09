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
      use definitions, only: iwp, wp, u6
      Use Para_Info, Only: MyRank, nProcs, Is_Real_Par
      use TList_Mod
      use stdalloc, only: mma_allocate
      use Constants, only: Zero, One, Two
      implicit real*8 (a-h,o-z)
      Logical(kind=iwp), intent(in):: Triangular
      real(kind=wp), intent(in):: P_Eff

      real(kind=wp)  distrib,PQpTsk,TskLw,TskHi,MinPQ1,a,fint,
     &               tskmin,tskmax
* parameters concerning task distribution...
      Integer(kind=iwp), Parameter:: iDen_PQ  = 2
      Integer(kind=iwp), Parameter:: iDen_Tsk = 4
      Integer(kind=iwp), Parameter:: MinPQ    = 4
* max number of tasks in tasklist per node...
      Integer(kind=iwp), Parameter:: MxnTsk = 100

      fint(a)=dble(int(a))

      If (Allocated(TskL)) Return
      MinPQ1= MinPQ
      MxnTsk1=MxnTsk
      iDen_PQ1=iDen_PQ
      iDen_Tsk1=iDen_Tsk
c
      P = P_Eff
      If (Triangular) Then
         PQ = P*(P+One)/Two
      Else
         PQ = P*P
      End If
      nTasks = nint(Min(PQ,dble(MxnTsk1*nProcs)))
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
*
      Call mma_allocate(TskM,2,nTasks,Label='TskM')
      TskM(:,:)=0
      Call mma_allocate(TskQ,2,nTasks,Label='TskQ')
      TskQ(:,:)=Not_Used
      Call mma_allocate(TskL,nTasks*2,Label='TskL')
*
      tskmin=1.d14
      tskmax=Zero
      TskLw=One
      TskHi=Zero
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
          PQpTsk=Max(MinPQ1,fint((PQ+dble(nTasks)-Two)/dble(nTasks-1)))
          kTskHi=nint(fint(PQ/PQpTsk))
          PQ=PQ-dble(kTskHi)*PQpTsk
        Else If (PQ.gt.Zero) Then
          PQpTsk=PQ
          kTskHi=1
          PQ=0d0
        Else
          kTskHi=0
          Write (u6,*) 'Init_TList: you should not be here!'
          Call Abend()
        End If
        Do kTsk = 1, kTskHi
          TskHi=TskHi+PQpTsk
          iTsk=iTsk+1
          TskM(1,iTsk)=TskLw
          TskM(2,iTsk)=TskHi
          tskmin=min(tskmin,(TskHi-TskLw+One))
          tskmax=max(tskmax,(TskHi-TskLw+One))
          TskLw=TskHi+One
        End Do
        nTasks=nTasks-kTskHi
      If (abs(PQ).gt.1.d-10) Go To 100
*     UNTIL (PQ == 0)
      If (nTasks.lt.0) Then
          Write (u6,*) 'nTasks.lt.0'
          Write (u6,*) 'MyRank=',MyRank
          Call Abend()
      End If
      nTasks=iTsk
*
      End Subroutine Init_TList
*
      Subroutine Free_TList()
      use TList_Mod, only: TskQ,TskM
      Use Para_Info, Only: nProcs, Is_Real_Par
      use stdalloc, only: mma_deallocate
      implicit None
*
      If (.Not.Allocated(TskQ)) Return
*
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
      Call mma_deallocate(TskQ)
      Call mma_deallocate(TskM)
*
      End Subroutine Free_TList
