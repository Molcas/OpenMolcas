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

subroutine Init_TList(Triangular,P_Eff)

use definitions, only: iwp, wp, u6
use Para_Info, only: MyRank, nProcs, Is_Real_Par
use TList_Mod, only: nTasks, P, PQ, TskL, TskM, TskQ, Not_Used
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Two

implicit none
logical(kind=iwp), intent(in) :: Triangular
real(kind=wp), intent(in) :: P_Eff
real(kind=wp) distrib, PQpTsk, TskLw, TskHi, MinPQ1, a, fint, tskmin, tskmax
! parameters concerning task distribution...
integer(kind=iwp), parameter :: iDen_PQ = 2
integer(kind=iwp), parameter :: iDen_Tsk = 4
integer(kind=iwp), parameter :: MinPQ = 4
! max number of tasks in tasklist per node...
integer(kind=iwp), parameter :: MxnTsk = 100
integer(kind=iwp) iDen_PQ1, iDen_Tsk1, iTsk, kTsk, kTskHi, MxnTsk1, nTaskpP, nTaskpP_seg

fint(a) = dble(int(a))

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
nTasks = nint(min(PQ,dble(MxnTsk1*nProcs)))
if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return

call mma_allocate(TskM,2,nTasks,Label='TskM')
TskM(:,:) = 0
call mma_allocate(TskQ,2,nTasks,Label='TskQ')
TskQ(:,:) = Not_Used
call mma_allocate(TskL,nTasks*2,Label='TskL')

tskmin = 1.d14
tskmax = Zero
TskLw = One
TskHi = Zero
iTsk = 0

! REPEAT
do
  distrib = fint(PQ/dble(iDen_PQ1))
  nTaskpP = nTasks/nProcs
  nTaskpP_seg = nTaskpP/iDen_Tsk1
  PQpTsk = MinPQ1
  if (nTaskpP_seg >= 1) PQpTsk = fint(distrib/dble(nTaskpP_seg*nProcs))
  if (PQpTsk > MinPQ1) then
    PQpTsk = fint(distrib/dble(nTaskpP_seg*nProcs))
    distrib = fint(PQpTsk*dble(nTaskpP_seg*nProcs))
    PQ = PQ-distrib
    kTskHi = nTaskpP_seg*nProcs
  else if ((PQ > MinPQ1) .and. (nTasks > 1)) then
    PQpTsk = max(MinPQ1,fint((PQ+dble(nTasks)-Two)/dble(nTasks-1)))
    kTskHi = nint(fint(PQ/PQpTsk))
    PQ = PQ-dble(kTskHi)*PQpTsk
  else if (PQ > Zero) then
    PQpTsk = PQ
    kTskHi = 1
    PQ = 0d0
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
  if (abs(PQ) <= 1.d-10) exit

end do

if (nTasks < 0) then
  write(u6,*) 'nTasks < 0'
  write(u6,*) 'MyRank=',MyRank
  call Abend()
end if
nTasks = iTsk

end subroutine Init_TList
