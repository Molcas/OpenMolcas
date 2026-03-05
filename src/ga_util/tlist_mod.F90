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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module TList_Mod

use Para_Info, only: Is_Real_Par, MyRank, nProcs
use stdalloc, only: mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp) :: iStrt_TList, iTCnSt, iTskCan, mTasks, nTasks
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iEnd_TList, igaTsk
#endif
real(kind=wp) :: P, PQ, QLast(2)
logical(kind=iwp) :: GT_Status = .false., PP_Status = .false.
integer(kind=iwp), allocatable, target :: TskL(:)
real(kind=wp), allocatable :: TskM(:,:), TskQ(:,:)
real(kind=wp), parameter :: Not_Used = -One

public :: Free_GTList, Free_PPList, Free_TList, Init_GTList, Init_PPList, Init_TList, Pos_QLast, Put_QLast, QLast, ReInit_GTList, &
          ReInit_PPList, Rsv_GTList

contains

subroutine Init_TList(Triangular,P_Eff)

  use stdalloc, only: mma_allocate
  use Constants, only: Zero, Two

  logical(kind=iwp), intent(in) :: Triangular
  real(kind=wp), intent(in) :: P_Eff
  integer(kind=iwp) :: iDen_PQ1, iDen_Tsk1, iTsk, kTsk, kTskHi, MxnTsk1, nTaskpP, nTaskpP_seg
  real(kind=wp) distrib, PQpTsk, TskLw, TskHi, MinPQ1, tskmin, tskmax
  integer(kind=iwp), parameter :: iDen_PQ = 2, iDen_Tsk = 4, MinPQ = 4, MxnTsk = 100

  if (allocated(TskL)) return
  MinPQ1 = MinPQ
  MxnTsk1 = MxnTsk
  iDen_PQ1 = iDen_PQ
  iDen_Tsk1 = iDen_Tsk

  P = P_Eff
  if (Triangular) then
    PQ = P*(P+One)/Two
  else
    PQ = P*P
  end if
  nTasks = nint(min(PQ,real(MxnTsk1*nProcs,kind=wp)))
  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return

  call mma_allocate(TskM,2,nTasks,Label='TskM')
  TskM(:,:) = 0
  call mma_allocate(TskQ,2,nTasks,Label='TskQ')
  TskQ(:,:) = Not_Used
  call mma_allocate(TskL,nTasks*2,Label='TskL')

  tskmin = 1.0e14_wp
  tskmax = Zero
  TskLw = One
  TskHi = Zero
  iTsk = 0

  ! REPEAT
  do
    distrib = fint(PQ/real(iDen_PQ1,kind=wp))
    nTaskpP = nTasks/nProcs
    nTaskpP_seg = nTaskpP/iDen_Tsk1
    PQpTsk = MinPQ1
    if (nTaskpP_seg >= 1) PQpTsk = fint(distrib/real(nTaskpP_seg*nProcs,kind=wp))
    if (PQpTsk > MinPQ1) then
      PQpTsk = fint(distrib/real(nTaskpP_seg*nProcs,kind=wp))
      distrib = fint(PQpTsk*real(nTaskpP_seg*nProcs,kind=wp))
      PQ = PQ-distrib
      kTskHi = nTaskpP_seg*nProcs
    else if ((PQ > MinPQ1) .and. (nTasks > 1)) then
      PQpTsk = max(MinPQ1,fint((PQ+real(nTasks,kind=wp)-Two)/real(nTasks-1,kind=wp)))
      kTskHi = int(PQ/PQpTsk)
      PQ = PQ-real(kTskHi,kind=wp)*PQpTsk
    else if (PQ > Zero) then
      PQpTsk = PQ
      kTskHi = 1
      PQ = Zero
    else
      kTskHi = 0
      write(u6,*) 'Init_TList: you should not be here!'
      call Abend()
    end if

    do kTsk=1,kTskHi
      TskHi = TskHi+PQpTsk
      iTsk = iTsk+1
      TskM(1,iTsk) = TskLw
      TskM(2,iTsk) = TskHi
      tskmin = min(tskmin,(TskHi-TskLw+One))
      tskmax = max(tskmax,(TskHi-TskLw+One))
      TskLw = TskHi+One
    end do
    nTasks = nTasks-kTskHi

    ! until (PQ == 0)
    ! if (abs(PQ) > 1.0e-10_wp) cycle
    ! exit
    if (abs(PQ) <= 1.0e-10_wp) exit

  end do

  if (nTasks < 0) then
    write(u6,*) 'nTasks < 0'
    write(u6,*) 'MyRank=',MyRank
    call Abend()
  end if
  nTasks = iTsk

end subroutine Init_TList

subroutine Free_TList()

  if (.not. allocated(TskQ)) return

  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
  call mma_deallocate(TskQ)
  call mma_deallocate(TskM)

end subroutine Free_TList

subroutine Init_GTList()

  if (GT_Status) return
  GT_Status = .true.

  iTCnSt = 1
  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
# ifdef _MOLCAS_MPP_
  ! create global tasklist...
  call GATskL(.true.,nTasks,igaTsk)
# endif

end subroutine Init_GTList

subroutine ReInit_GTList()

  if (.not. GT_Status) then
    write(u6,*) 'ReInit_GTList: List not active!'
    call Abend()
  end if
  iTCnSt = 1
  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
# ifdef _MOLCAS_MPP_
  ! initialize global tasklist...
  call GATskL_Zero(igaTsk)
# endif

end subroutine ReInit_GTList

subroutine Free_GTList()

  if (.not. GT_Status) return
  GT_Status = .false.

  iTCnSt = 1
  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
# ifdef _MOLCAS_MPP_
  ! create global tasklist...
  call GATskL(.false.,nTasks,igaTsk)
# endif

end subroutine Free_GTList

function Rsv_GTList(TskLw,TskHi,iOpt,NewBatch)

  logical(kind=iwp) :: Rsv_GTList
  real(kind=wp), intent(out) :: TskLw, TskHi
  integer(kind=iwp), intent(in) :: iOpt
  logical(kind=iwp), intent(out) :: NewBatch
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: MyTask
  integer(kind=iwp), external :: RsvTsk
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Rsv_GTList = .false.
  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) then
    if (iTCnSt == 1) then
      iTCnSt = iTCnSt+1
      iTskCan = iTskCan+1
      Rsv_GTList = .true.
      TskLw = One
      TskHi = PQ
      iStrt_TList = 1
      if (iOpt == 0) then
        NewBatch = .true.
      else
        NewBatch = .false.
      end if
    end if
# ifdef _MOLCAS_MPP_
  else
    !
    !---- NewBatch is true if
    !      1) Batch never processed by the node before
    !      2) The sequence of executed batches is broken

    if (iOpt == 0) then
      MyTask = RsvTsk(igaTsk,TskL,nTasks,nTasks,iTCnST,iStrt_TList,iEnd_TList)
      NewBatch = .true.
    else if (iOpt == 1) then
      MyTask = RsvTsk(igaTsk,TskL,nTasks,mTasks,iTCnST,iStrt_TList,iEnd_TList)
      NewBatch = iStrt_TList > mTasks
    else if (iOpt == 2) then
      MyTask = RsvTsk(igaTsk,TskL,nTasks,nTasks,iTCnST,iStrt_TList,iEnd_TList)
      NewBatch = (iStrt_TList > mTasks) .or. (iStrt_TList /= iTCnST)
    else
      MyTask = 0
      write(u6,*) 'Rsv_GTList: Invalid option:',iOpt
      call Abend()
    end if
    if (MyTask >= 1) then
      Rsv_GTList = .true.
      TskLw = TskM(1,MyTask)
      TskHi = TskM(2,MyTask)
      iTCnSt = iTCnSt+1
      iTskCan = iTskCan+1
    end if
# endif
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *

end function Rsv_GTList

subroutine Init_PPList()

  integer(kind=iwp) :: i, iE, iTsk
  integer(kind=iwp), pointer :: TskList(:,:)
  logical(kind=iwp), parameter :: Debug = .false.

  if (Debug) then
    if (PP_Status) then
      write(u6,*) 'Init_PPList: Active'
    else
      write(u6,*) 'Init_PPList: InActive'
    end if
  end if

  if (PP_Status) return
  PP_Status = .true.

  iTskCan = 0
  mTasks = 0
  iStrt_TList = 0
# ifdef _MOLCAS_MPP_
  iEnd_TList = nTasks+1
# endif
  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return

  TskList(1:nTasks,1:2) => TskL(1:2*nTasks)
  TskList(:,1) = 0
  do iTsk=0,nTasks-1
    TskList(1+iTsk,1) = mod(iTsk+MyRank,nTasks)+1
  end do

  ! Copy list in inverse order to ensure the proper order if
  ! reinitiated at once.

  iE = nTasks-1
  TskList(:,2) = 0
  do i=0,nTasks-1
    TskList(1+iE,2) = TskList(i+1,1)
    iE = iE-1
  end do
  nullify(TskList)

  QLast(1) = Not_Used
  QLast(2) = Not_Used

end subroutine Init_PPList

subroutine ReInit_PPList(Semi_Direct)

  logical(kind=iwp), intent(in) :: Semi_Direct
  integer(kind=iwp) i, iCount, iE
  integer(kind=iwp), pointer :: TskList(:,:)
  logical(kind=iwp), parameter :: Debug = .false.

  if (Debug) then
    if (PP_Status) then
      write(u6,*) 'ReInit_PPList: Active'
    else
      write(u6,*) 'ReInit_PPList: InActive'
    end if
  end if
  if (.not. PP_Status) then
    write(u6,*) 'ReInit_PPList: List is not active!'
    call Abend()
  end if
  iTskCan = 0
  mTasks = iStrt_TList
  if (nProcs == 1) then
    iStrt_TList = 0
#   ifdef _MOLCAS_MPP_
    iEnd_TList = nTasks+1
#   endif
  else
    if (Semi_Direct) then
      TskList(1:nTasks,1:2) => TskL(1:2*nTasks)

      ! Copy first the task indices of tasks that were exectuted
      call ICopy(mTasks,TskList(:,2),1,TskList(:,1),1)

      ! Now copy task indices of tasks which were not executed by this node.
      ! Change the order so that the first task it the largest in the list.

      iE = mTasks
      iCount = 1
      do i=mTasks,nTasks-1
        if (iCount > myRank) then
          TskList(1+i,1) = TskList(i+1,2)
        else
          TskList(1+i,1) = TskList(iE+1,2)
          iE = iE-1
          iCount = iCount+1
        end if
      end do

      nullify(TskList)
    end if

    iStrt_TList = 0
#   ifdef _MOLCAS_MPP_
    iEnd_TList = nTasks+1
#   endif

    QLast(1) = Not_Used
    QLast(2) = Not_Used
  end if

end subroutine ReInit_PPList

subroutine Free_PPList()

  if (.not. allocated(TskL)) return
  PP_Status = .false.

  if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
  call mma_deallocate(TskL)

end subroutine Free_PPList

subroutine Put_QLast()

  if (.not. allocated(TskQ)) return
  TskQ(1,iTskCan) = QLast(1)
  TskQ(2,iTskCan) = QLast(2)

  QLast(1) = Not_Used
  QLast(2) = Not_Used

end subroutine Put_QLast

subroutine Pos_QLast(Disc)

  use Definitions, only: RtoI

  real(kind=wp), intent(inout) :: Disc
  integer(kind=iwp) :: iWR(2), mInts
  real(kind=wp) :: Dummy(1), Quad_ijkl, RST_triplet
  logical(kind=iwp), parameter :: Copy = .true., NoCopy = .false.

  if (.not. allocated(TskQ)) return

  Quad_ijkl = TskQ(1,iTskCan)
  RST_triplet = TskQ(2,iTskCan)
  if (Quad_ijkl == Not_Used) return

  ! If already at the right position return

  if ((Quad_ijkl == QLast(1)) .and. (RST_triplet == QLast(2))) return

  do

    call iRBuf(iWR,2,Copy)
    call dRBuf(QLast,2,Copy)
    mInts = iWR(2)
    if ((QLast(1) == Quad_ijkl) .and. (QLast(2) == RST_triplet)) then
      if (mInts > 0) call dRBuf(Dummy,mInts,NoCopy)
      Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
      return
    else if (QLast(1) <= Quad_ijkl) then
      if (mInts > 0) call dRBuf(Dummy,mInts,NoCopy)
      Disc = Disc+real(2/RtoI+2+mInts,kind=wp)
      cycle
    else
      write(u6,*) 'Pos_QLast: batch is lost!'
      write(u6,'(A,2F10.1)') 'Index,1.0:  ',QLast(1),QLast(2)
      write(u6,'(A,2F10.1)') 'Looking for ',Quad_ijkl,RST_triplet
      write(u6,*) ' iTskCan,=',iTskCan
      call RecPrt('TskQ',' ',TskQ,2,iTskCan)
      write(u6,*)
      call Abend()
    end if

  end do

  write(u6,*) 'Pos_QLast: Fatal problem!'
  call Abend()

end subroutine Pos_QLast

elemental function fint(a)

  real(kind=wp) :: fint
  real(kind=wp), intent(in) :: a

  fint = real(int(a),kind=wp)

end function fint

end module TList_Mod
