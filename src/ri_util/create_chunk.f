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
      Subroutine Create_Chunk(LenVec,NumVec,IncVec)
      use Chunk_mod
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: MyRank, nProcs, Is_Real_Par
#endif
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
      External  ga_create_irreg
      Logical   ga_create_irreg, ok
      Integer myMap(2)


*
      If (NumVec.le.0) Then
         Call WarningMessage(2,'Create_Chunk: Failure NumVec.le.0')
         Write (6,*) 'NumVec=',NumVec
         Call Abend()
      End If

      If (Is_Real_Par()) Then
         Call mma_allocate(iMap,nProcs+1,Label='iMap')
         iMap(:)=0
*
         FullSize=DBLE(LenVec*NumVec)
         Call mma_maxDBLE(MaxMem)
         iMap(1+MyRank) = MaxMem
         Call GAIGOP(iMap,nProcs,'+')
         TotalMemory=0.0D0
         itmp=iMap(1)
*
*        Find the smallest possible memory allocation!
*
         Do i = 1, nProcs-1
            itmp = Min(itmp,iMap(1+i))
         End Do
         TotalMemory=DBLE(itmp)*DBLE(nProcs)
*
*        Compute the number of vectors to handle at the time
*
         If (TotalMemory.gt.FullSize) Then
*
            IncVec = NumVec
*
         Else
*
            IncVec = INT ( DBLE(NumVec) *(TotalMemory/FullSize))
*
         End If
         If (IncVec.le.0) Then
            Call WarningMessage(2,
     &                  'Create_Chunk: Failure IncVec.le.0')
            Write (6,*) 'FullSize=',FullSize
            Write (6,*) 'NumVec=',NumVec
            Write (6,*) 'LenVec=',LenVec
            Write (6,*) 'TotalMemory=',TotalMemory
            Write (6,*) 'Local size of memory'
            Write (6,*) (iMap(i),i=1,nProcs)
            Write (6,*) 'iTmp=',iTmp
            Call Abend()
         End If
*
*        Compute the number of vectors per node, This also defines
*        the Map array.
*
         iNode0=0
         iStart = 1
         Do iNode = 0, nProcs-1
            If (iStart.eq.1) iNode0=iNode
            iMap(1+iNode) = iStart
            iStart = iStart + Max((IncVec+iNode)/nProcs,1)
         End Do
         iMap(1+nProcs)=iStart
         IncVec0=iStart
*
C        Call Put_iArray('DistVec',iMap,nProcs+1)
*
         nBlocks=nProcs-iNode0
         myMap(1)=1
         Ok = GA_Create_Irreg(mt_dbl,LenVec,IncVec0,'Chunk',
     &                        myMap,1,
     &                        iMap(1+iNode0),nBlocks,ip_Chunk)
         If (.Not. Ok) Then
            Call WarningMessage(2,'Error in GA_Create_Irreg')
            Call Abend()
         End If
      Else
         Call mma_maxDBLE(MaxMem)
         IncVec=Min(NumVec,MaxMem/LenVec)
         Call mma_allocate(Chunk,LenVec*IncVec,Label='Chunk')
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
#else
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_maxDBLE(MaxMem)
      IncVec=Min(NumVec,MaxMem/LenVec)
      Call mma_allocate(Chunk,LenVec*IncVec,Label='Chunk')
*
#endif
      Return
      End
