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
      Subroutine Create_Chunk(ip_iMap,ip_Chunk,LenVec,NumVec,IncVec)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: MyRank, nProcs, Is_Real_Par
#endif
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
      External  ga_create_irreg
      Logical   ga_create_irreg, ok
      Integer myMap(2)
*
      If (Is_Real_Par()) Then
         Call GetMem('iMap','Allo','Inte',ip_iMap,nProcs+1)
         Call IZero(iWork(ip_iMap),nProcs+1)
*
         FullSize=DBLE(LenVec*NumVec)
         Call GetMem('MemMax','Max','Real',iDummy,MaxMem)
         iWork(ip_iMap+MyRank) = MaxMem
         Call GAIGOP(iWork(ip_iMap),nProcs,'+')
         TotalMemory=0.0D0
         itmp=iWork(ip_iMap)
*
*        Find the smallest possible memory allocation!
*
         Do i = 1, nProcs-1
            itmp = Min(itmp,iWork(ip_iMap+i))
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
         If (NumVec.le.0) Then
            Call WarningMessage(2,
     &                  'Create_Chunk: Failure NumVec.le.0')
            Write (6,*) 'NumVec=',NumVec
            Call Abend()
         End If
         If (IncVec.le.0) Then
            Call WarningMessage(2,
     &                  'Create_Chunk: Failure IncVec.le.0')
            Write (6,*) 'FullSize=',FullSize
            Write (6,*) 'NumVec=',NumVec
            Write (6,*) 'LenVec=',LenVec
            Write (6,*) 'TotalMemory=',TotalMemory
            Write (6,*) 'Local size of memory'
            Write (6,*) (iWork(ip_iMap+i),i=0,nProcs-1)
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
            iWork(ip_iMap+iNode) = iStart
            iStart = iStart + Max((IncVec+iNode)/nProcs,1)
         End Do
         iWork(ip_iMap+nProcs)=iStart
*
C        Call Put_iArray('DistVec',iWork(ip_iMap),nProcs+1)
*
         nBlocks=nProcs-iNode0
         myMap(1)=1
         Ok = GA_Create_Irreg(mt_dbl,LenVec,IncVec,'Chunk',
     &                        myMap,1,
     &                        iWork(ip_iMap+iNode0),nBlocks,
     &                        ip_Chunk)
         If (.Not. Ok) Then
            Call WarningMessage(2,'Error in GA_Create_Irreg')
            Call Abend()
         End If
      Else
         ip_iMap=ip_iDummy
         Call GetMem('MemMax','Max','Real',iDummy,MaxMem)
         IncVec=Min(NumVec,MaxMem/LenVec)
         Call GetMem('Chunk','Allo','Real',ip_Chunk,LenVec*IncVec)
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
#else
*                                                                      *
************************************************************************
*                                                                      *
      ip_iMap=ip_iDummy
      Call GetMem('MemMax','Max','Real',iDummy,MaxMem)
      IncVec=Min(NumVec,MaxMem/LenVec)
      Call GetMem('Chunk','Allo','Real',ip_Chunk,LenVec*IncVec)
*
#endif
      Return
      End
