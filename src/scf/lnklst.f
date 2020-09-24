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
* Copyright (C) 1994-1996, Martin Schuetz                              *
************************************************************************
************************************************************************
* This Module comprises 8 subroutines and 5 functions, which are all   *
* used to store, recall and administrate vectors from e.g. subsequent  *
* iterations, which can be assigned to an ordered keylist, e.g. to a   *
* list of iteration numbers.                                           *
* The internal data structure used to control the vectors is a linked  *
* list of nodes, each containing a pointer to the memory or disk       *
* location of the vector, together with some other information.        *
* The head of the list is a special node called CNOD which contains    *
* the following info:                                                  *
*     CNOD:    addr(CNOD)   -> error code, if 0: op performed OK       *
*              addr(CNOD)+1 -> ptr to first NODE in list (iroot)       *
*              addr(CNOD)+2 -> actual length of list (lislen)          *
*              addr(CNOD)+3 -> # of vectors wished to keep in memory   *
*                              (incore)                                *
*              addr(CNOD)+4 -> unused yet                              *
*              addr(CNOD)+5 -> unused yet                              *
* One NODE of the linked list itself contains six integers and is      *
* organized as follows:                                                *
*     NODE:    addr(NODE)   -> ptr to next node                        *
*              addr(NODE)+1 -> ptr to stored vector                    *
*              addr(NODE)+2 -> next free position, if on disk          *
*                              unused otherwise                        *
*              addr(NODE)+3 -> length of vector                        *
*              addr(NODE)+4 -> iteration number                        *
*              addr(NODE)+5 -> 1, if in core; 0, if on disk            *
* The 11 subroutines and functions operating on the linked lists have  *
* the following purposes:                                              *
* IniLst(LList,incore)                                                 *
*  -> initialize list                                                  *
* PutVec(vec,lvec,LUnit,iterat,NoAllo,opcode,LList)                    *
*  -> store vector on list, ev. move the tailnode vector on disk       *
* GetVec(LUnit,iterat,LList,inode,vec,lvec)                            *
*  -> fetch vector from list, which corresponds to iterat              *
* GetNod(iterat,LList,inode)                                           *
*  -> does not read out vector, just returns address of node in inode  *
* InfNod(inode,iterat,ipnext,ipvec,lvec)                               *
*  -> returns info of node indicated by inode                          *
* Logical InCore(inode)                                                *
*  -> returns true, if corresponding vector is incore, false otherwise *
* Integer iVPtr(LUnit,ivptr1,inode)                                    *
*  -> uses InfNod or GetVec to read out vector of inode, depending     *
*     if InCore or not. If the result is ivptr1, the vector was read   *
*     from disk. ivptr1 has to be allocated before (same vector length *
*     as stored in inode).                                             *
* Integer LstPtr(LUnit,iterat,LList)                                   *
*  -> uses GetNod and InfNod to get the pointer to the vector, which   *
*     corresponds to iterat. The pointer is the return value of the    *
*     function. If the vector is not InCore, the function terminates   *
*     with an Error message. This function is used only in cases,      *
*     where the vector for sure is stored in core (e.g. last entries   *
*     in LList).                                                       *
* Logical LLErr(LList)                                                 *
*  -> checks, if ErrCode was set in previous LL Operation              *
* Integer LLLen(LList)                                                 *
*  -> returns the actual length of the LL                              *
* KilLst(LList)                                                        *
*  -> kill list and free all memory                                    *
* DmpLst(LList,LUnit,lDskPt)                                           *
*  -> flush list and dump it on corresponding DA listfile              *
*     lDskPt: DA pointer where listnodes begin (on output)             *
* RclLst(LList,LUnit,lDskPt,NoAllo)                                    *
*  -> recall and restore list from corresponding DA listfile           *
*     lDskPt: DA pointer where listnodes begin  (on input)             *
*     NoAllo: amount of memory not to allocate...                      *
*----------------------------------------------------------------------*
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1994/96                              *
************************************************************************


      SubRoutine IniLst(iLList,incore)
      Implicit Real*8 (a-h,o-z)
      Integer iLList,incore
#include "real.fh"

#include "mxdm.fh"
#include "lnklst.fh"
*
      Debug_LnkLst=.False.
*
*     allocate list header CNOD
      lLList=lLList+1
      iLList=lLList
      nLList(iLList,0)=0
      nLList(iLList,1)=0
      nLList(iLList,2)=0
      nLList(iLList,3)=incore

*
      If (Debug_LnkLst) Then
         Write (6,*) 'IniLst'
         Call StlLst(iLList)
      End If
*
      Return
      End
*----------------------------------------------------------------------*


      SubRoutine PutVec(vec,lvec,LUnit,iterat,NoAllo,opcode,iLList)
*     NoAllo is the amount of memory (in DWords) one wants to keep
*     for other purposes.
*     opcode is a 4 character string:
*     'APND' -append nevertheless to list, if same iterat is found
*             on head node
*     'NOOP' -no operation, if same iterat is found on head node
*     'OVWR' -overwrite vector, if same iterat is found on head node
*             if lvec is different from lvec stored on node, old
*             vector is not overwritten and iWork(LList) is set to
*             ErrCode 1
*
      Implicit Real*8 (a-h,o-z)
*
*     declaration subroutine parameters
      Integer lvec,LUnit,iterat,iLList,NoAllo,iroot,lislen,incore
      Real*8 vec(lvec)
      Character opcode*4
*
*     declaration local variables
      Integer iPtr1,iPtr2,MaxMem
C     Integer iDskPt,len
*
#include "real.fh"
#include "WrkSpc.fh"
#include "mxdm.fh"
#include "lnklst.fh"

#include "SysDef.fh"
*
      If (Debug_LnkLst) Then
         Write (6,*) 'PutVec'
         Call StlLst(iLList)
         Call qTrace
      End If
#ifdef _DEBUGPRINT_
      Call qEnter('PutVec')
#endif
*
*     clear ErrCode
      nLList(iLList,0)=0
*     read listhead
      iroot=nLList(iLList,1)
      lislen=nLList(iLList,2)
      incore=nLList(iLList,3)

      If ((iroot.gt.0).AND.(nLList(iroot,4).eq.iterat)) Then
        If (opcode.eq.'NOOP') Then
*         that's all, folks
#ifdef _DEBUGPRINT_
          Call qExit('PutVec')
#endif
          Return
        Else If (opcode.eq.'OVWR') Then
          If (nLList(iroot,3).ne.lvec) Then
* Set error code: inconsistency in vector lengths
            nLList(iLList,0)=1
          Else
            call dcopy_(lvec,vec,1,Work(nLList(iroot,1)),1)
          End If
#ifdef _DEBUGPRINT_
          Call qExit('PutVec')
#endif
          Return
        Else If (opcode.ne.'APND') Then
*         opcode unknown
          Write (6,*) 'PutVec: opcode unknown'
          Write (6,'(A,A)') 'opcode=',opcode
          Call QTrace
          Call Abend()
        End If
      End If
*     check if there is still enough memory to store vec
      Call GetMem('LVec ','Max','Real',iPtr1,MaxMem)
cvv Enough memory
*       let's allocate some memory
        Call GetMem('LVec ','Allo','Real',iPtr1,lvec)
*      End If
*     allocate new node
      lLList=lLList+1
      iPtr2=lLList
      nLList(iPtr2,0)=iroot
      nLList(iPtr2,1)=iPtr1
      nLList(iPtr2,2)=0
      nLList(iPtr2,3)=lvec
      nLList(iPtr2,4)=iterat
      nLList(iPtr2,5)=1


      call dcopy_(lvec,vec,1,Work(iPtr1),1)
      iroot=iPtr2
      lislen=lislen+1
      nLList(iLList,1)=iroot
      nLList(iLList,2)=lislen

*
#ifdef _DEBUGPRINT_
      Call qExit('PutVec')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(LUnit)
         Call Unused_integer(NoAllo)
      End If
      End
*----------------------------------------------------------------------*


      SubRoutine GetVec(LUnit,iterat,iLList,inode,vec,lvec)
*     searches linked list for node corresponding to iterat, starting
*     from iroot=iWork(LList+1).
*     inode points to the node found after searching.
*     the vector stored at inode then is copied to vec, but only
*     if lvec is equal to the vector length stored in inode. Otherwise
*     the entry is considered as inconsistent.
*     inode is set to zero, if no entry was found, or to its negative
*     address, if an inconsistent entry was found.
*     if LList<0, then -LList is interpreted as a direct node address,
*     and not the address of the listhead (faster access).
      Implicit Real*8 (a-h,o-z)
*
*     declaration subroutine parameters
      Integer lvec,LUnit,iterat,iLList,inode
      Real*8 vec(lvec)
*
*     declaration local variables
c      Integer iDskPt
*
#include "real.fh"
#include "WrkSpc.fh"
#include "mxdm.fh"
#include "lnklst.fh"

#include "SysDef.fh"
*
#ifdef _DEBUGPRINT_
      Call qEnter('GetVec')
#endif
*
        inode=nLList(iLList,1)

 100  If ((nLList(inode,4).ne.iterat).and.(nLList(inode,0).ne.0)) Then
        inode=nLList(inode,0)

        GoTo 100
      End If
      If (nLList(inode,4).eq.iterat) Then
*       we've found matching entry, so check if consistent
        If (nLList(inode,3).eq.lvec) Then
*         everything's allright, we made it, let's copy to vec
            call dcopy_(lvec,Work(nLList(inode,1)),1,vec,1)
        Else
* inconsistency
          write(6,*)' Found inconsistency.'
*          inode=-inode
          inode=0
        End If
      Else
        inode=0
      End If
*
#ifdef _DEBUGPRINT_
      Call qExit('GetVec')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LUnit)
      End
*----------------------------------------------------------------------*


      SubRoutine GetNod(iterat,iLList,inode)
*     searches linked list for node corresponding to iterat, starting
*     from iroot. inode points to the node found after searching.
*     inode is set to zero and iWork(LList)=0 is set to ErrCode 1,
*     if no correspondance was found.
      Implicit Real*8 (a-h,o-z)
*
*     declaration subroutine parameters
      Integer iterat,iLList,inode
*
#include "real.fh"
#include "mxdm.fh"
#include "lnklst.fh"

*
#ifdef _DEBUGPRINT_
      Call qEnter('GetNod')
#endif
*
      If (Debug_LnkLst) Then
         Write (6,*) 'GetNod'
         Call StlLst(iLList)
      End If
*
*     clear ErrCode
      nLList(iLList,0)=0
*     set inode to iroot
      inode=nLList(iLList,1)

 100  If ((nLList(inode,4).ne.iterat).and.(nLList(inode,0).ne.0)) Then
        inode=nLList(inode,0)
        GoTo 100
      End If
      If (nLList(inode,4).eq.iterat) Then
*       we've found matching entry
      Else
        Write (6,*) 'GetNod: Warning!'
        inode=0
        nLList(iLList,0)=1
      End If
*
#ifdef _DEBUGPRINT_
      Call qExit('GetNod')
#endif
      Return
      End
*----------------------------------------------------------------------*


      SubRoutine InfNod(inode,iterat,ipnext,ipvec,lvec)
*     returns info of node indicated by inode. iterat,ipnext,ipvec,lvec
*     are overwritten with the corresponding info on the node
      Implicit Real*8 (a-h,o-z)
*
*     declaration of procedure parameters
      Integer inode,iterat,ipnext,ipvec,lvec
*

#include "mxdm.fh"
#include "lnklst.fh"
*
      iterat=nLList(inode,4)
      ipnext=nLList(inode,0)
      ipvec= nLList(inode,1)
      lvec=  nLList(inode,3)
*
      Return
      End
*----------------------------------------------------------------------*


      Logical Function InCore(inode)
*     returns true, if corresponding vector is incore, false otherwise
      Implicit Real*8 (a-h,o-z)
      Integer inode
*

#include "mxdm.fh"
#include "lnklst.fh"
*
      If (nLList(inode,5).eq.1) Then
        InCore=.TRUE.
      Else
        InCore=.FALSE.
      End If
      Return
      End
*----------------------------------------------------------------------*


      Logical Function LLErr(iLList)
*     checks, if ErrCode was set in previous LL Operation
      Implicit Real*8 (a-h,o-z)
      Integer iLList
*

#include "mxdm.fh"
#include "lnklst.fh"
*

      If (nLList(iLList,0).eq.0) Then
        LLErr=.FALSE.
      Else
        LLErr=.TRUE.
      End If
      Return
      End
*----------------------------------------------------------------------*


      Integer Function LLLen(iLList)
*     returns the actual length of the LL
      Implicit Real*8 (a-h,o-z)
      Integer iLList
*

#include "mxdm.fh"
#include "lnklst.fh"
*
      LLLen=nLList(iLList,2)
      Return
      End
*----------------------------------------------------------------------*


      Subroutine iVPtr(LUnit,vptr1,nvptr1,inode)
*     uses InfNod or GetVec to read out vector of inode, depending
*     if InCore or not. If the result is ivptr1, the vector was read
*     from disk. ivptr1 has to be allocated before (same vector length
*     as stored in inode). The return value points to the memory
*     location of the read vector. If the function fails (in GetVec),
*     the inode value of GetVec is returned.
*
*     2017-03-15:Converted to return the array in vptr1.
      Implicit Real*8 (a-h,o-z)
      Integer LUnit,nvptr1,ivptr2,inode,idum
      Logical InCore
      Real*8  vptr1(nvptr1)
*
#include "real.fh"
#include "WrkSpc.fh"
#include "mxdm.fh"
#include "lnklst.fh"
*
      If (InCore(inode)) Then
        Call InfNod(inode,idum,idum,ivptr2,idum)
        Call DCopy_(nvptr1,Work(ivptr2),1,vptr1,1)
      Else
        Call GetVec(LUnit,nLList(inode,4),inode,inode,
     &              vptr1,nLList(inode,3))
      End If
*
      Return
      End
*----------------------------------------------------------------------*


      Integer Function LstPtr(LUnit,iterat,iLList)
*     uses GetNod and InfNod to obtain the pointer to the vector, which
*     corresponds to iterat. The pointer is the return value of the
*     function. If the vector is not InCore, the function terminates
*     with an Error message. This function is used only in cases,
*     where the vector for sure is stored in core (e.g. last entries
*     in LList).
      Implicit Real*8 (a-h,o-z)
*
*     declaration subroutine parameters
      Integer LUnit,iterat,iLList
*
*     declaration local variables
      Integer inode,idum,ivptr
*     and functions
      Logical InCore
*
#ifdef _DEBUGPRINT_
      Call qEnter('LstPtr')
#endif
      LstPtr=-999999
      Call GetNod(iterat,iLList,inode)
      If (inode.eq.0) Then
* Hmmm, no entry found in LList, that's strange
        Write (6,*) 'LstPtr: inode.le.0'
        Write (6,*) 'inode=',inode
        Call QTrace
        Call Abend()
      Else If (InCore(inode)) Then
        Call InfNod(inode,idum,idum,ivptr,idum)
        LstPtr=ivptr
      Else
* Hmmm, no incore hit for this entry, that's strange
        Write (6,*) 'LstPtr: no incore hit for this entry'
        Write (6,*) 'inode=',inode
        Call QTrace
        Call Abend()
      End If
#ifdef _DEBUGPRINT_
      Call qExit('LstPtr')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LUnit)
      End
*----------------------------------------------------------------------*


      SubRoutine KilLst(iLList)
      Implicit Real*8 (a-h,o-z)
*     Free all memory of linked list LList
#include "real.fh"
#include "WrkSpc.fh"
#include "mxdm.fh"
#include "lnklst.fh"
*     local vars
      Integer iLList,iroot,iPtr1
#ifdef _DEBUGPRINT_
      Call qEnter('KilLst')
#endif
*
*
      If (Debug_LnkLst) Then
         Write (6,*) 'KilLst'
         Call StlLst(iLList)
      End If
*
      iroot=nLList(iLList,1)
 100  Continue
      If (iroot.ne.0) Then
        iPtr1=nLList(iroot,1)
        iFlag=nLList(iroot,5)

        If (iFlag.eq.1) Then
          Call GetMem('LVec ','Free','Real',iPtr1,nLList(iroot,3))
        End If
        iPtr1=iroot
        iroot=nLList(iroot,0)
        GoTo 100
      End If
*
#ifdef _DEBUGPRINT_
      Call qExit('KilLst')
#endif
      Return
      End
*----------------------------------------------------------------------*


      SubRoutine DmpLst(iLList,LUnit,lDskPt)
      Implicit Real*8 (a-h,o-z)
      Integer iLList,LUnit,lDskPt
*
#include "WrkSpc.fh"
#include "mxdm.fh"
#include "lnklst.fh"

#include "SysDef.fh"
*
#ifdef _DEBUGPRINT_
      Call QEnter('DmpLst')
#endif
*     clear ErrCode
      nLList(iLList,0)=0
*     read listhead
      iroot=nLList(iLList,1)
      lislen=nLList(iLList,2)
      incore=nLList(iLList,3)

      If (iroot.le.0) Then
* linked list has zero length, that's strange
* save in any case listhead...
        lDskPt=0
        iDskPt=lDskPt
        Call iDaFile(LUnit,1,nLList(iLList,0),NodSiz,iDskPt)
*       Call GetMem('CNOD ','Free','Inte',LList,NodSiz)
#ifdef _DEBUGPRINT_
        Call QExit('DmpLst')
#endif
        Return
      End If
 10   Continue
      If (nLList(iroot,5).eq.1) Then
        iPtr1=iroot
        iPtr2=iPtr1
* go either to last element or list or to last element in core
* and flush vector
 100    If ((nLList(iPtr1,0).ne.0).and.(nLList(iPtr1,5).eq.1)) Then
          iPtr2=iPtr1
          iPtr1=nLList(iPtr1,0)
          GoTo 100
        End If
        If (nLList(iPtr1,5).eq.1) Then
* nothing written on disk yet
          iDskPt=0
          iPtr2=iPtr1
          iPtr1=nLList(iPtr1,1)
        Else
          iDskPt=nLList(iPtr1,2)
          iPtr1=nLList(iPtr2,1)
        End If
* now set up node again and write vec to disk
        nLList(iPtr2,1)=iDskPt
        nLList(iPtr2,5)=0
        len=nLList(iPtr2,3)

        Call dDaFile(LUnit,1,Work(iPtr1),len,iDskPt)
        Call GetMem('LVec ','Free','Real',iPtr1,len)
        Go To 10
      End If
      lDskPt=iDskPt
*     now all vectors are flushed... so dump linked list...
      Call iDaFile(LUnit,1,nLList(iLList,0),NodSiz,iDskPt)
*
 200  If (iroot.ne.0) Then
        iPtr1=iroot
        iroot=nLList(iroot,0)
        Call iDaFile(LUnit,1,nLList(iPtr1,0),NodSiz,iDskPt)
*       Call GetMem('LNode','Free','Inte',iPtr1,NodSiz)
        GoTo 200
      End If
*     Call GetMem('CNOD ','Free','Inte',LList,NodSiz)
*
#ifdef _DEBUGPRINT_
      Call QExit('DmpLst')
#endif
      Return
      End
*----------------------------------------------------------------------*


      SubRoutine RclLst(iLList,LUnit,lDskPt,NoAllo)
      Implicit Real*8 (a-h,o-z)
      Integer iLList,LUnit,lDskPt,NoAllo
*
#include "WrkSpc.fh"
#include "mxdm.fh"
#include "lnklst.fh"

#include "SysDef.fh"
*
#ifdef _DEBUGPRINT_
      Call QEnter('RclLst')
#endif
* load listhead...
      lLList=lLList+1
      iLList=lLList

      Call iDaFile(LUnit,2,nLList(iLList,0),NodSiz,lDskPt)
      iroot=nLList(iLList,1)

      If (iroot.le.0) Then
        Write (6,*) 'RclLst: linked list has zero length,'
     &            //' that''s strange!'
* linked list has zero length, that's strange
*       Call Quit(20)
#ifdef _DEBUGPRINT_
        Call QExit('RclLst')
#endif
        Return
      End If
*
* restore linked list...
      lislen=1
      lLList=lLList+1
      iPtr1=lLList

      nLList(iLList,1)=iPtr1
      Call iDaFile(LUnit,2,nLList(iPtr1,0),NodSiz,lDskPt)

      iroot=iPtr1
      iPtr2=iroot
 100  If (nLList(iPtr2,0).ne.0) Then
        lislen=lislen+1
         lLList=lLList+1
         iPtr1=lLList
        nLList(iPtr2,0)=iPtr1

        Call iDaFile(LUnit,2,nLList(iPtr1,0),NodSiz,lDskPt)
        iPtr2=iPtr1
        Go To 100
      End If
      If (nLList(iLList,2).ne.lislen) Then
        Write(6,*) 'RclLst:LList length mismatch:',
     &              nLList(iLList,2),lislen
        Call Abend
      End If
      Write (6,*) 'Let''s restore...'
* now we have restored the list, let's fetch some vectors
      incore=nLList(iLList,3)
      Call GetMem('LVec ','Max','Real',iPtr1,MaxMem)
      lvec=nLList(iroot,3)
      iPtr2=iroot
 200  If ((incore.gt.0).AND.(MaxMem-NoAllo.ge.lvec).AND.(iPtr2.gt.0))
     &  Then
        lDskPt=nLList(iPtr2,1)
        Call GetMem('LVec ','Allo','Real',iPtr1,lvec)
        Call dDaFile(LUnit,2,Work(iPtr1),lvec,lDskPt)
        nLList(iPtr2,1)=iPtr1
        nLList(iPtr2,2)=0
        nLList(iPtr2,5)=1
        iPtr2=nLList(iPtr2,0)
        incore=incore-1

        Call GetMem('LVec ','Max','Real',iPtr1,MaxMem)
        Go To 200
      End If
      If (iPtr2.gt.0) nLList(iLList,3)=nLList(iLList,3)-incore
*
#ifdef _DEBUGPRINT_
      Call QExit('RclLst')
#endif
      Return
      End
