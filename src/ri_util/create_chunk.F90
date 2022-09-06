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

subroutine Create_Chunk(LenVec,NumVec,IncVec)

use Chunk_mod
#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs, Is_Real_Par
#endif

implicit real*8(A-H,O-Z)
#include "stdalloc.fh"
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
external ga_create_irreg
logical ga_create_irreg, ok
integer myMap(2)
#endif

#ifdef _MOLCAS_MPP_
if (NumVec <= 0) then
  call WarningMessage(2,'Create_Chunk: Failure NumVec <= 0')
  write(6,*) 'NumVec=',NumVec
  call Abend()
end if

if (Is_Real_Par()) then
  call mma_allocate(iMap,nProcs+1,Label='iMap')
  iMap(:) = 0

  FullSize = dble(LenVec*NumVec)
  call mma_maxDBLE(MaxMem)
  iMap(1+MyRank) = MaxMem
  call GAIGOP(iMap,nProcs,'+')
  TotalMemory = 0.0d0
  itmp = iMap(1)

  ! Find the smallest possible memory allocation!

  do i=1,nProcs-1
    itmp = min(itmp,iMap(1+i))
  end do
  TotalMemory = dble(itmp)*dble(nProcs)

  ! Compute the number of vectors to handle at the time

  if (TotalMemory > FullSize) then

    IncVec = NumVec

  else

    IncVec = int(dble(NumVec)*(TotalMemory/FullSize))

  end if
  if (IncVec <= 0) then
    call WarningMessage(2,'Create_Chunk: Failure IncVec <= 0')
    write(6,*) 'FullSize=',FullSize
    write(6,*) 'NumVec=',NumVec
    write(6,*) 'LenVec=',LenVec
    write(6,*) 'TotalMemory=',TotalMemory
    write(6,*) 'Local size of memory'
    write(6,*) (iMap(i),i=1,nProcs)
    write(6,*) 'iTmp=',iTmp
    call Abend()
  end if

  ! Compute the number of vectors per node, This also defines the Map array.

  iNode0 = 0
  iStart = 1
  do iNode=0,nProcs-1
    if (iStart == 1) iNode0 = iNode
    iMap(1+iNode) = iStart
    iStart = iStart+max((IncVec+iNode)/nProcs,1)
  end do
  iMap(1+nProcs) = iStart
  IncVec0 = iStart

  !call Put_iArray('DistVec',iMap,nProcs+1)

  nBlocks = nProcs-iNode0
  myMap(1) = 1
  Ok = GA_Create_Irreg(mt_dbl,LenVec,IncVec0,'Chunk',myMap,1,iMap(1+iNode0),nBlocks,ip_Chunk)
  if (.not. Ok) then
    call WarningMessage(2,'Error in GA_Create_Irreg')
    call Abend()
  end if
else
  call mma_maxDBLE(MaxMem)
  IncVec = min(NumVec,MaxMem/LenVec)
  call mma_allocate(Chunk,LenVec*IncVec,Label='Chunk')
end if

!                                                                      *
!***********************************************************************
!                                                                      *
#else
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_maxDBLE(MaxMem)
IncVec = min(NumVec,MaxMem/LenVec)
call mma_allocate(Chunk,LenVec*IncVec,Label='Chunk')

#endif

return

end subroutine Create_Chunk
