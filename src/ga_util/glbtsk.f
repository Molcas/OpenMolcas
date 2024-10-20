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
************************************************************************
* SubRoutine GATskL(create,nTsk,igaTsk)                                *
*  -> initialize or kill global task list                              *
*     create:   Logical: .TRUE. -> create / .FALSE. -> kill            *
*     nTsk:     # of tasks                                             *
*     igaTsk:   global array handle to global task list (on return)    *
************************************************************************
      SubRoutine GATskL(create,nTsk,igaTsk)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: Is_Real_Par
#endif
      Implicit None
      Logical create
      Integer nTsk,igaTsk
#ifdef _MOLCAS_MPP_
#  include "global.fh"
#  include "mafdecls.fh"
      Logical ok
      Integer nProcs,Chk
*
      If (.Not. Is_Real_Par()) Return
      If (create) Then
        nProcs=ga_nnodes()
        Chk=(nTsk+nProcs-1)/nProcs
        ok=ga_create(MT_INT,nTsk,1,'GlTskL',Chk,1,igaTsk)
        If (.NOT.ok) Then
          Write (6,*) 'GATskL: ga_create not OK!'
          Call GAStp('GATskL',42)
          Call Abend()
        End If
        Call ga_zero(igaTsk)
      Else
        ok=ga_destroy(igaTsk)
      End If
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(create)
         Call Unused_integer(nTsk)
         Call Unused_integer(igaTsk)
      End If
#endif
      End SubRoutine GATskL
*----------------------------------------------------------------------*
      SubRoutine GATskL_Zero(igaTsk)
      Implicit None
      Integer igaTsk
#ifdef _MOLCAS_MPP_
#  include "global.fh"

      Call ga_zero(igaTsk)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(igaTsk)
#endif
      End SubRoutine GATskL_Zero
************************************************************************
* Integer Function RsvTsk(igaTsk,iTskLs,nTsk,iStart)                   *
*  -> reserve a task for my node and mark it on the global task list   *
*     as reserved.                                                     *
*     igaTsk:   global array handle to global task list                *
*     iTskLs:   my private task list (my favourite sequence of tasks)  *
*     nTsk:     # of tasks                                             *
*     iStart:   starting value of index of private task list           *
************************************************************************
      Integer Function RsvTsk(igaTsk,iTskLs,nTsk,mTsk,iStart,iS,iE)
#ifdef _MOLCAS_MPP_
#  ifdef _GA_
      Implicit None
      Integer nTsk,mTsk,igaTsk,iTskLs(nTsk,2),iStart,iS,iE
#  include "global.fh"
      Logical Reserved
      Integer iCnt,iTsk
*
      If (iStart.gt.mTsk) Then
        iTsk=0
        iCnt=nTsk
      Else
        Do iCnt = iStart, mTsk
          iTsk = iTskLs(iCnt,1)
*         try to reserve iTsk on global task list...
          Reserved = ga_read_inc(igaTsk,iTsk,1,1) .ne. 0
          If (Reserved) Then
             iE = iE - 1
             iTskLs(iE,2)=iTsk
          Else
             iS = iS + 1
             iTskLs(iS,2)=iTsk
             Go To 100
          End If
        End Do
        iTsk=0
  100   Continue
*
      End If
      RsvTsk=iTsk
      iStart=iCnt
#  else
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nTsk,mTsk,igaTsk,iTskLs(nTsk,2),iStart,iS,iE
#  include "global.fh"
      Logical Reserved
      Integer iCnt,iTsk
      Integer :: One=1
      Integer, allocatable:: TSKR(:)
*
      If (iStart.gt.mTsk) Then
        iTsk=0
        iCnt=nTsk
      Else
        Call mma_allocate(TSKR,mTsk,Label='TSKR')
        Call ga_readb_inc(igaTsk,mTsk,iTskLs(iStart,1),TSKR)
        Do iCnt = iStart, mTsk
          iTsk = iTskLs(iCnt,1)
          Reserved = TSKR(1+iCnt-iStart) .ne. 0
          If (Reserved) Then
             iE = iE - 1
             iTskLs(iE,2)=iTsk
          Else
             iS = iS + 1
             iTskLs(iS,2)=iTsk
             Go To 100
          End If
        End Do
        iTsk=0
  100   Continue
        Call mma_deallocate(TSKR)
      End If
      RsvTsk=iTsk
      iStart=iCnt
#  endif
#else
      Implicit None
      Integer nTsk,mTsk,igaTsk,iTskLs(nTsk,2),iStart,iS,iE
      RsvTsk=0
      iStart=nTsk
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(igaTsk)
         Call Unused_integer_array(iTskLs)
         Call Unused_integer(mTsk)
         Call Unused_integer(iS)
         Call Unused_integer(iE)
      End If
#endif
      Return
      End

