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
! Copyright (C) 1994-1996, Martin Schuetz                              *
!***********************************************************************
!***********************************************************************
! This Module comprises 8 subroutines and 5 functions, which are all   *
! used to store, recall and administrate vectors from e.g. subsequent  *
! iterations, which can be assigned to an ordered keylist, e.g. to a   *
! list of iteration numbers.                                           *
! The internal data structure used to control the vectors is a linked  *
! list of nodes, each containing a pointer to the memory or disk       *
! location of the vector, together with some other information.        *
! The head of the list is a special node called CNOD which contains    *
! the following info:                                                  *
!     CNOD:    addr(CNOD)   -> error code, if 0: op performed OK       *
!              addr(CNOD)+1 -> ptr to first NODE in list (iroot)       *
!              addr(CNOD)+2 -> actual length of list (lislen)          *
!              addr(CNOD)+3 -> # of vectors wished to keep in memory   *
!                              (incore)                                *
!              addr(CNOD)+4 -> unused yet                              *
!              addr(CNOD)+5 -> unused yet                              *
! One NODE of the linked list itself contains six integers and is      *
! organized as follows:                                                *
!     NODE:    addr(NODE)   -> ptr to next node                        *
!              addr(NODE)+1 -> ptr to stored vector                    *
!              addr(NODE)+2 -> next free position, if on disk          *
!                              unused otherwise                        *
!              addr(NODE)+3 -> length of vector                        *
!              addr(NODE)+4 -> iteration number                        *
!              addr(NODE)+5 -> 1, if in core; 0, if on disk            *
! The 11 subroutines and functions operating on the linked lists have  *
! the following purposes:                                              *
! IniLst(LList,incore)                                                 *
!  -> initialize list                                                  *
! PutVec(vec,lvec,iterat,opcode,LList)                                 *
!  -> store vector on list, ev. move the tailnode vector on disk       *
! GetVec(iterat,LList,inode,vec,lvec)                                  *
!  -> fetch vector from list, which corresponds to iterat              *
! GetNod(iterat,LList,inode)                                           *
!  -> does not read out vector, just returns address of node in inode  *
! InfNod(inode,iterat,ipnext,ipvec,lvec)                               *
!  -> returns info of node indicated by inode                          *
! Logical InCore(inode)                                                *
!  -> returns true, if corresponding vector is incore, false otherwise *
! Integer iVPtr(vptr1,nvtr1,inode)                                     *
!  -> uses InfNod or GetVec to read out vector of inode, depending     *
!     if InCore or not. If the result is ivptr1, the vector was read   *
!     from disk. ivptr1 has to be allocated before (same vector length *
!     as stored in inode).                                             *
! Integer LstPtr(iterat,LList)                                         *
!  -> uses GetNod and InfNod to get the pointer to the vector, which   *
!     corresponds to iterat. The pointer is the return value of the    *
!     function. If the vector is not InCore, the function terminates   *
!     with an Error message. This function is used only in cases,      *
!     where the vector for sure is stored in core (e.g. last entries   *
!     in LList).                                                       *
! Logical LLErr(LList)                                                 *
!  -> checks, if ErrCode was set in previous LL Operation              *
! Integer LLLen(LList)                                                 *
!  -> returns the actual length of the LL                              *
! KilLst(LList)                                                        *
!  -> kill list and free all memory                                    *
! DmpLst(LList,LUnit,lDskPt)                                           *
!  -> flush list and dump it on corresponding DA listfile              *
!     lDskPt: DA pointer where listnodes begin (on output)             *
! RclLst(LList,LUnit,lDskPt,NoAllo)                                    *
!  -> recall and restore list from corresponding DA listfile           *
!     lDskPt: DA pointer where listnodes begin  (on input)             *
!     NoAllo: amount of memory not to allocate...                      *
!----------------------------------------------------------------------*
!     written by:                                                      *
!     M. Schuetz                                                       *
!     University of Lund, Sweden, 1994/96                              *
!***********************************************************************
      SubRoutine IniLst(iLList,incore)
      use LnkLst, only: Debug_LnkLst, lLList, nLList
      Implicit None
      Integer iLList,incore

      Debug_LnkLst=.False.
!
!     allocate list header CNOD
      lLList=lLList+1
      iLList=lLList
      nLList(iLList,0)=0
      nLList(iLList,1)=0
      nLList(iLList,2)=0
      nLList(iLList,3)=incore

!
      If (Debug_LnkLst) Then
         Write (6,*) 'IniLst'
         Call StlLst(iLList)
      End If
!
      Return
      End SubRoutine IniLst
!----------------------------------------------------------------------*

      SubRoutine PutVec(vec,lvec,iterat,opcode,iLList)
!     NoAllo is the amount of memory (in DWords) one wants to keep
!     for other purposes.
!     opcode is a 4 character string:
!     'NOOP' -no operation, if same iterat is found on head node
!     'OVWR' -overwrite vector, if same iterat is found in any node
!             If lvec is different from lvec stored on node, old
!             vector is not overwritten and iWork(LList) is set to
!             ErrCode 1
!
      use LnkLst, only: Debug_LnkLst, lLList, nLList, SCF_V, MaxNodes
      use stdalloc, only: mma_allocate
      Implicit None
!
!     declaration subroutine parameters
      Integer lvec,iterat,iLList,iroot,lislen
      Real*8 vec(lvec)
      Character opcode*4
!
!     declaration local variables
      Integer iPtr2,MaxMem
!     Integer iDskPt,len
!
#include "SysDef.fh"
!
      If (Debug_LnkLst) Then
         Write (6,*) 'PutVec'
         Call StlLst(iLList)
      End If
!
!     clear ErrCode
      nLList(iLList,0)=0
!     read listhead
      iroot =nLList(iLList,1)
      lislen=nLList(iLList,2)

      Select Case(opcode)

      Case('NOOP')

         If ((iroot.gt.0).AND.(nLList(iroot,4).eq.iterat)) Return

      Case('OVWR')
         Do While ((iroot>0))
            If (nLList(iroot,3).ne.lvec) Then
! Set error code: inconsistency in vector lengths
               nLList(iLList,0)=1
            Else If (nLList(iroot,4)==iterat) Then
               SCF_V(iroot)%A(1:lVec)=vec(1:lVec)
               Return
            End If
            iroot=nLList(iroot,0)
         End Do

      Case default

!         opcode unknown
          Write (6,*) 'PutVec: opcode unknown'
          Write (6,'(A,A)') 'opcode=',opcode
          Call Abend()

      End Select

      iroot =nLList(iLList,1)
!     check if there is still enough memory to store vec
      Call mma_maxDBLE(MaxMem)
!     let's allocate some memory
!     allocate new node
      lLList=lLList+1
      iPtr2=lLList
      If (iPtr2.gt.Maxnodes) Then
         Write (6,*) 'PutVec: iPtr2.gt.Maxnodes'
         Call Abend()
      End If
      If (Allocated(SCF_V(iPtr2)%A)) Then
         Write (6,*) 'Node already allocated'
         Write (6,*) 'iPtr2=',iPtr2
         Call Abend()
      End If
      Call mma_allocate(SCF_V(iPtr2)%A,lVec,Label='LVec')
      nLList(iPtr2,0)=iroot
      nLList(iPtr2,1)=iPtr2
      nLList(iPtr2,2)=0
      nLList(iPtr2,3)=lvec
      nLList(iPtr2,4)=iterat
      nLList(iPtr2,5)=1


      SCF_V(iPtr2)%A(:)=Vec(:)
      iroot=iPtr2
      lislen=lislen+1
      nLList(iLList,1)=iroot
      nLList(iLList,2)=lislen

!
      Return
      End SubRoutine PutVec
!----------------------------------------------------------------------*


      SubRoutine GetVec(iterat,iLList,inode,vec,lvec)
!     searches linked list for node corresponding to iterat, starting
!     from iroot=iWork(LList+1).
!     inode points to the node found after searching.
!     the vector stored at inode then is copied to vec, but only
!     if lvec is equal to the vector length stored in inode. Otherwise
!     the entry is considered as inconsistent.
!     inode is set to zero, if no entry was found, or to its negative
!     address, if an inconsistent entry was found.
!     if LList<0, then -LList is interpreted as a direct node address,
!     and not the address of the listhead (faster access).
      use LnkLst, only: nLList, SCF_V
      Implicit None
!
!     declaration subroutine parameters
      Integer lvec,iterat,iLList,inode
      Real*8 vec(lvec)
!
#include "SysDef.fh"
!
      inode=nLList(iLList,1)
      If (inode<=0) Then
         Write (6,*) 'GetVec: iNode<=0'
         Call Abend()
      End If

      Do While ((nLList(inode,4).ne.iterat).and.(nLList(inode,0).ne.0))
         inode=nLList(inode,0)
      End Do

      If (nLList(inode,4).eq.iterat) Then
!       we've found matching entry, so check if consistent
        If (nLList(inode,3).eq.lvec) Then
!          everything's alright, we made it, let's copy to vec
           vec(1:lVec)=SCF_V(iNode)%A(1:lVec)
        Else
! inconsistency
          write(6,*)' Found inconsistency.'
!          inode=-inode
          inode=0
        End If
      Else
        inode=0
      End If
!
      End SubRoutine GetVec
!----------------------------------------------------------------------*


      SubRoutine GetNod(iterat,iLList,inode)
!     searches linked list for node corresponding to iterat, starting
!     from iroot. inode points to the node found after searching.
!     inode is set to zero and iWork(LList)=0 is set to ErrCode 1,
!     if no correspondance was found.
      use LnkLst, only: Debug_LnkLst, nLList
      Implicit None
!
!     declaration subroutine parameters
      Integer iterat,iLList,inode
!
      If (Debug_LnkLst) Then
         Write (6,*) 'GetNod'
         Call StlLst(iLList)
      End If
!
!     clear ErrCode
      nLList(iLList,0)=0
!     set inode to iroot
      inode=nLList(iLList,1)
      If (inode<=0) Then
         Write (6,*) 'GetNod: iNode<=0'
         Write (6,*) 'iLList=',iLList
         Call Abend()
      End If

      Do While ((nLList(inode,4).ne.iterat).and.(nLList(inode,0).ne.0))
        inode=nLList(inode,0)
      End Do
      If (nLList(inode,4).eq.iterat) Then
!       we've found matching entry
      Else
        Write (6,*) 'GetNod: Warning!'
        inode=0
        nLList(iLList,0)=1
      End If
!
      End SubRoutine GetNod
!----------------------------------------------------------------------*


      SubRoutine InfNod(inode,iterat,ipnext,ipvec,lvec)
!     returns info of node indicated by inode. iterat,ipnext,ipvec,lvec
!     are overwritten with the corresponding info on the node
      use LnkLst, Only: nLList
      Implicit None
!
!     declaration of procedure parameters
      Integer inode,iterat,ipnext,ipvec,lvec
!

      iterat=nLList(inode,4)
      ipnext=nLList(inode,0)
      ipvec= nLList(inode,1)
      lvec=  nLList(inode,3)
!
      Return
      End SubRoutine InfNod
!----------------------------------------------------------------------*


      Logical Function InCore(inode)
!     returns true, if corresponding vector is incore, false otherwise
      use LnkLst, Only: nLList
      Implicit None
      Integer inode
!

      If (nLList(inode,5).eq.1) Then
        InCore=.TRUE.
      Else
        InCore=.FALSE.
      End If
      Return
      End Function InCore
!----------------------------------------------------------------------*


      Logical Function LLErr(iLList)
!     checks, if ErrCode was set in previous LL Operation
      use LnkLst, Only: nLList
      Implicit None
      Integer iLList
!

      If (nLList(iLList,0).eq.0) Then
        LLErr=.FALSE.
      Else
        LLErr=.TRUE.
      End If
      Return
      End Function LLErr
!----------------------------------------------------------------------*


      Integer Function LLLen(iLList)
!     returns the actual length of the LL
      use LnkLst, Only: nLList
      Implicit None
      Integer iLList
!

      LLLen=nLList(iLList,2)
      Return
      End Function LLLen
!----------------------------------------------------------------------*


      Subroutine iVPtr(vptr1,nvptr1,inode)
!     uses InfNod or GetVec to read out vector of inode, depending
!     if InCore or not. If the result is ivptr1, the vector was read
!     from disk. ivptr1 has to be allocated before (same vector length
!     as stored in inode). The return value points to the memory
!     location of the read vector. If the function fails (in GetVec),
!     the inode value of GetVec is returned.
!
!     2017-03-15:Converted to return the array in vptr1.
      use LnkLst, only: SCF_V, nLLIst
      Implicit None
      Integer nvptr1,ivptr2,inode,idum1,idum2,idum3
      Logical, External:: InCore
      Real*8  vptr1(nvptr1)
!
      If (InCore(inode)) Then
        Call InfNod(inode,idum1,idum2,ivptr2,idum3)
        vPtr1(1:nvptr1)=SCF_V(inode)%A(1:nvptr1)
      Else
        Call GetVec(nLList(inode,4),inode,inode,vptr1,nLList(inode,3))
      End If
!
      Return
      End Subroutine iVPtr
!----------------------------------------------------------------------*


      Integer Function LstPtr(iterat,iLList)
!     uses GetNod and InfNod to obtain the pointer to the vector, which
!     corresponds to iterat. The pointer is the return value of the
!     function. If the vector is not InCore, the function terminates
!     with an Error message. This function is used only in cases,
!     where the vector for sure is stored in core (e.g. last entries
!     in LList).
      Implicit None
!
!     declaration subroutine parameters
      Integer iterat,iLList
!
!     declaration local variables
      Integer inode,idum1,idum2,idum3,ivptr
!     and functions
      Logical, External::  InCore
!
      LstPtr=-999999
      Call GetNod(iterat,iLList,inode)
      If (inode.eq.0) Then
! Hmmm, no entry found in LList, that's strange
        Write (6,*) 'LstPtr: inode.le.0'
        Write (6,*) 'inode=',inode
        Call Abend()
      Else If (InCore(inode)) Then
        Call InfNod(inode,idum1,idum2,ivptr,idum3)
        LstPtr=ivptr
      Else
! Hmmm, no incore hit for this entry, that's strange
        Write (6,*) 'LstPtr: no incore hit for this entry'
        Write (6,*) 'inode=',inode
        Call Abend()
      End If
      End Function LstPtr
!----------------------------------------------------------------------*


      SubRoutine KilLst(iLList)
      use LnkLst, only: Debug_LnkLst, nLList, SCF_V
      use stdalloc, only: mma_deallocate
      Implicit None
!     Free all memory of linked list LList
!     local vars
      Integer iLList,iroot, iFlag
!
!
      If (Debug_LnkLst) Then
         Write (6,*) 'KilLst'
         Call StlLst(iLList)
      End If
!
      iroot=nLList(iLList,1)
      Do While (iroot.ne.0)
        iFlag=nLList(iroot,5)

        If (iFlag.eq.1) Then
          Call mma_deallocate(SCF_V(iroot)%A)
        End If
        iroot=nLList(iroot,0)
      End Do
!
      End SubRoutine KilLst
!----------------------------------------------------------------------*


      SubRoutine DmpLst(iLList,LUnit,lDskPt)
      use LnkLst, only: SCF_V, nLList, NodSiz
      use stdalloc, only: mma_deallocate
      Implicit None
      Integer iLList,LUnit,lDskPt
!
#include "SysDef.fh"
      Integer iDskPt, iPtr1, iPtr2, iRoot, Len
!
!     clear ErrCode
      nLList(iLList,0)=0
!     read listhead
      iroot=nLList(iLList,1)

      If (iroot.le.0) Then
! linked list has zero length, that's strange
! save in any case listhead...
        lDskPt=0
        iDskPt=lDskPt
        Call iDaFile(LUnit,1,nLList(iLList,0),NodSiz,iDskPt)
        Return
      End If

      Do While (nLList(iroot,5).eq.1)
         iPtr1=iroot
         iPtr2=iPtr1
! go either to last element or list or to last element in core
! and flush vector
         Do While ((nLList(iPtr1,0).ne.0).and.(nLList(iPtr1,5).eq.1))
          iPtr2=iPtr1
          iPtr1=nLList(iPtr1,0)
        End Do
        If (nLList(iPtr1,5).eq.1) Then
! nothing written on disk yet
          iDskPt=0
          iPtr2=iPtr1
          iPtr1=nLList(iPtr1,1)
        Else
          iDskPt=nLList(iPtr1,2)
          iPtr1=nLList(iPtr2,1)
        End If
! now set up node again and write vec to disk
        nLList(iPtr2,1)=iDskPt
        nLList(iPtr2,5)=0
        len=nLList(iPtr2,3)

        Call dDaFile(LUnit,1,SCF_V(iPtr2)%A,len,iDskPt)
        Call mma_deallocate(SCF_V(iPtr2)%A)
      End Do
      lDskPt=iDskPt
!     now all vectors are flushed... so dump linked list...
      Call iDaFile(LUnit,1,nLList(iLList,0),NodSiz,iDskPt)
!
      Do While (iroot.ne.0)
        iPtr1=iroot
        iroot=nLList(iroot,0)
        Call iDaFile(LUnit,1,nLList(iPtr1,0),NodSiz,iDskPt)
      End Do
!
      End SubRoutine DmpLst
!----------------------------------------------------------------------*


      SubRoutine RclLst(iLList,LUnit,lDskPt,NoAllo)
      use LnkLst, only: nLList, SCF_V, lLList, NodSiz, MaxNodes
      use stdalloc, only: mma_allocate
      Implicit None
      Integer iLList,LUnit,lDskPt,NoAllo
!
#include "SysDef.fh"
      Integer iPtr1, iPtr2, iRoot, lislen, MaxMem, lVec, incore
!
! load listhead...
      lLList=lLList+1
      iLList=lLList

      Call iDaFile(LUnit,2,nLList(iLList,0),NodSiz,lDskPt)
      iroot=nLList(iLList,1)

      If (iroot.le.0) Then
        Write (6,*) 'RclLst: linked list has zero length, that''s strange!'
! linked list has zero length, that's strange
!       Call Quit(20)
        Return
      End If
!
! restore linked list...
      lislen=1
      lLList=lLList+1
      iPtr1=lLList

      nLList(iLList,1)=iPtr1
      Call iDaFile(LUnit,2,nLList(iPtr1,0),NodSiz,lDskPt)

      iroot=iPtr1
      iPtr2=iroot
      Do While (nLList(iPtr2,0).ne.0)
        lislen=lislen+1
         lLList=lLList+1
         iPtr1=lLList
        nLList(iPtr2,0)=iPtr1

        Call iDaFile(LUnit,2,nLList(iPtr1,0),NodSiz,lDskPt)
        iPtr2=iPtr1
      End Do
      If (nLList(iLList,2).ne.lislen) Then
        Write(6,*) 'RclLst:LList length mismatch:',nLList(iLList,2),lislen
        Call Abend()
      End If
      Write (6,*) 'Let''s restore...'
! now we have restored the list, let's fetch some vectors
      incore=nLList(iLList,3)
      Call mma_maxDBLE(MaxMem)
      lvec=nLList(iroot,3)
      iPtr2=iroot
      Do While ((incore.gt.0).AND.(MaxMem-NoAllo.ge.lvec).AND.(iPtr2.gt.0))
        lDskPt=nLList(iPtr2,1)

         If (iPtr2.gt.Maxnodes) Then
            Write (6,*) 'iPtr2.gt.Maxnodes, restoring'
            Call Abend()
        End If
        If (Allocated(SCF_V(iPtr2)%A)) Then
           Write (6,*) 'Node already allocated while restoring'
           Write (6,*) 'iPtr2=',iPtr2
           Call Abend()
        End If
        Call mma_Allocate(SCF_V(iPtr2)%A,lvec,Label='LVec')
        Call dDaFile(LUnit,2,SCF_V(iPtr2)%A,lvec,lDskPt)
        nLList(iPtr2,1)=iPtr2
        nLList(iPtr2,2)=0
        nLList(iPtr2,5)=1
        iPtr2=nLList(iPtr2,0)
        incore=incore-1

        Call mma_maxDBLE(MaxMem)
      End Do
      If (iPtr2.gt.0) nLList(iLList,3)=nLList(iLList,3)-incore
!
      End SubRoutine RclLst
