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
! Logical InCore_f(inode)                                              *
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
!----------------------------------------------------------------------*
!     Define linked lists for storage of vectors of subsequent iters   *
!----------------------------------------------------------------------*
!     LLGrad - linked list of gradients                                *
!     LLlGrd - linked list of local gradients                          *
!     LLdGrd - linked list of gradient diffs                           *
!     LLDelt - linked list of Delta vectors                            *
!     LLy    - linked list of y-vectors                                *
!     LLx    - linked list of x orbital rotation parameter vectors     *
!              (only for QNR/DIIS combination)                         *
!----------------------------------------------------------------------*

!#define _DEBUGPRINT_
module LnkLst

use InfSCF, only: MxIter
use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_maxDBLE
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp), parameter :: MAXnodes = (MxIter+1)*6, NodSiz = 6

integer(kind=iwp) :: LLDelt, LLdGrd, LLGrad, LLlGrd, lLList = 0, LLx, LLy, nLList(MAXnodes,0:NodSiz-1)
logical(kind=iwp) :: Init_LLs = .false.
type(Alloc1DArray_Type) :: SCF_V(Maxnodes)

public :: DmpLst, GetNod, GetVec, IniLst, Init_LLs, iVPtr, KilLst, LLDelt, LLdGrd, LLGrad, LLLen, LLlGrd, lLList, LLx, LLy, &
          LstPtr, nLList, NodSiz, PutVec, RclLst, SCF_V, StlLst

contains

subroutine IniLst(iLList,incore)

  integer(kind=iwp), intent(out) :: iLList
  integer(kind=iwp), intent(in) :: incore

  ! allocate list header CNOD
  lLList = lLList+1
  iLList = lLList
  nLList(iLList,0) = 0
  nLList(iLList,1) = 0
  nLList(iLList,2) = 0
  nLList(iLList,3) = incore

# ifdef _DEBUGPRINT_
  write(u6,*) 'IniLst'
  call StlLst(iLList)
# endif

  return

end subroutine IniLst

subroutine PutVec(vec,lvec,iterat,opcode,iLList)
  ! NoAllo is the amount of memory (in DWords) one wants to keep
  ! for other purposes.
  ! opcode is a 4 character string:
  ! 'NOOP' -no operation, if same iterat is found on head node
  ! 'OVWR' -overwrite vector, if same iterat is found in any node
  !         If lvec is different from lvec stored on node, old
  !         vector is not overwritten and iWork(LList) is set to
  !         ErrCode 1

  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: lvec, iterat, iLList
  real(kind=wp), intent(in) :: vec(lvec)
  character(len=4), intent(in) :: opcode
  integer(kind=iwp) :: iPtr2, iroot, lislen, MaxMem

# ifdef _DEBUGPRINT_
  write(u6,*) 'PutVec',iterat
  call StlLst(iLList)
# endif

  ! clear ErrCode
  nLList(iLList,0) = 0
  ! read listhead
  iroot = nLList(iLList,1)
  lislen = nLList(iLList,2)

  select case (opcode)

    case ('NOOP')

      if (iroot > 0) then
        if (nLList(iroot,4) == iterat) return
      endif

    case ('OVWR')
      do while ((iroot > 0))
        if (nLList(iroot,3) /= lvec) then
          ! Set error code: inconsistency in vector lengths
          nLList(iLList,0) = 1
        else if (nLList(iroot,4) == iterat) then
          SCF_V(iroot)%A(1:lVec) = vec(:)
          return
        end if
        iroot = nLList(iroot,0)
      end do

    case default

      ! opcode unknown
      write(u6,*) 'PutVec: opcode unknown'
      write(u6,'(A,A)') 'opcode=',opcode
      call Abend()

  end select

  iroot = nLList(iLList,1)
  ! check if there is still enough memory to store vec
  call mma_maxDBLE(MaxMem)
  ! let's allocate some memory
  ! allocate new node
  lLList = lLList+1
  iPtr2 = lLList
  if (iPtr2 > Maxnodes) then
    write(u6,*) 'PutVec: iPtr2 > Maxnodes'
    call Abend()
  end if
  if (allocated(SCF_V(iPtr2)%A)) then
    write(u6,*) 'Node already allocated'
    write(u6,*) 'iPtr2=',iPtr2
    call Abend()
  end if
  call mma_allocate(SCF_V(iPtr2)%A,lVec,Label='LVec')
  nLList(iPtr2,0) = iroot
  nLList(iPtr2,1) = iPtr2
  nLList(iPtr2,2) = 0
  nLList(iPtr2,3) = lvec
  nLList(iPtr2,4) = iterat
  nLList(iPtr2,5) = 1

  SCF_V(iPtr2)%A(:) = Vec(:)
  iroot = iPtr2
  lislen = lislen+1
  nLList(iLList,1) = iroot
  nLList(iLList,2) = lislen

  return

end subroutine PutVec

subroutine GetVec(iterat,iLList,inode,vec,lvec)
  ! searches linked list for node corresponding to iterat, starting
  ! from iroot=iWork(LList+1).
  ! inode points to the node found after searching.
  ! the vector stored at inode then is copied to vec, but only
  ! if lvec is equal to the vector length stored in inode. Otherwise
  ! the entry is considered as inconsistent.
  ! inode is set to zero, if no entry was found, or to its negative
  ! address, if an inconsistent entry was found.
  ! if LList<0, then -LList is interpreted as a direct node address,
  ! and not the address of the listhead (faster access).

  integer(kind=iwp), intent(in) :: iterat, iLList, lvec
  integer(kind=iwp), intent(out) :: inode
  real(kind=wp), intent(out) :: vec(lvec)

  inode = nLList(iLList,1)
  if (inode <= 0) then
    write(u6,*) 'GetVec: iNode<=0'
    call Abend()
  end if

  do while ((nLList(inode,4) /= iterat) .and. (nLList(inode,0) /= 0))
    inode = nLList(inode,0)
  end do

  if (nLList(inode,4) == iterat) then
    ! we've found matching entry, so check if consistent
    if (nLList(inode,3) == lvec) then
      ! everything's alright, we made it, let's copy to vec
      vec(:) = SCF_V(iNode)%A(1:lVec)
    else
      ! inconsistency
      write(u6,*) ' Found inconsistency.'
      !inode = -inode
      inode = 0
    end if
  else
    inode = 0
  end if

end subroutine GetVec

subroutine GetNod(iterat,iLList,inode)
  ! searches linked list for node corresponding to iterat, starting
  ! from iroot. inode points to the node found after searching.
  ! inode is set to zero and iWork(LList)=0 is set to ErrCode 1,
  ! if no correspondance was found.

  integer(kind=iwp), intent(in) :: iterat, iLList
  integer(kind=iwp), intent(out) :: inode

# ifdef _DEBUGPRINT_
  write(u6,*) 'GetNod',iterat
  call StlLst(iLList)
# endif

  ! clear ErrCode
  nLList(iLList,0) = 0
  ! set inode to iroot
  inode = nLList(iLList,1)
  if (inode <= 0) then
    write(u6,*) 'GetNod: iNode<=0'
    write(u6,*) 'iLList=',iLList
    call Abend()
  end if

  do while ((nLList(inode,4) /= iterat) .and. (nLList(inode,0) /= 0))
    inode = nLList(inode,0)
  end do
  if (nLList(inode,4) == iterat) then
    ! we've found matching entry
  else
    write(u6,*) 'GetNod: Warning!'
    inode = 0
    nLList(iLList,0) = 1
  end if

end subroutine GetNod

subroutine InfNod(inode,iterat,ipnext,ipvec,lvec)
  ! returns info of node indicated by inode. iterat,ipnext,ipvec,lvec
  ! are overwritten with the corresponding info on the node

  integer(kind=iwp), intent(in) :: inode
  integer(kind=iwp), intent(out) :: iterat, ipnext, ipvec, lvec

  iterat = nLList(inode,4)
  ipnext = nLList(inode,0)
  ipvec = nLList(inode,1)
  lvec = nLList(inode,3)

  return

end subroutine InfNod

function InCore_f(inode)
  ! returns true, if corresponding vector is incore, false otherwise

  logical(kind=iwp) :: InCore_f
  integer(kind=iwp), intent(in) :: inode

  if (nLList(inode,5) == 1) then
    InCore_f = .true.
  else
    InCore_f = .false.
  end if

  return

end function InCore_f

function LLErr(iLList)
  ! checks, if ErrCode was set in previous LL Operation

  logical(kind=iwp) :: LLErr
  integer(kind=iwp), intent(in) :: iLList

  if (nLList(iLList,0) == 0) then
    LLErr = .false.
  else
    LLErr = .true.
  end if

  return

end function LLErr

function LLLen(iLList)
  ! returns the actual length of the LL

  integer(kind=iwp) :: LLLen
  integer(kind=iwp), intent(in) :: iLList

  LLLen = nLList(iLList,2)

  return

end function LLLen

subroutine iVPtr(vptr1,nvptr1,inode)
  ! uses InfNod or GetVec to read out vector of inode, depending
  ! if InCore or not. If the result is ivptr1, the vector was read
  ! from disk. ivptr1 has to be allocated before (same vector length
  ! as stored in inode). The return value points to the memory
  ! location of the read vector. If the function fails (in GetVec),
  ! the inode value of GetVec is returned.
  !
  ! 2017-03-15:Converted to return the array in vptr1.

  integer(kind=iwp), intent(in) :: nvptr1, inode
  real(kind=wp), intent(out) :: vptr1(nvptr1)
  integer(kind=iwp) :: idum1, idum2, idum3, ivptr2

  if (InCore_f(inode)) then
    call InfNod(inode,idum1,idum2,ivptr2,idum3)
    vPtr1(:) = SCF_V(inode)%A(1:nvptr1)
  else
    call GetVec(nLList(inode,4),inode,idum1,vptr1,nLList(inode,3))
  end if

  return

end subroutine iVPtr

function LstPtr(iterat,iLList)
  ! uses GetNod and InfNod to obtain the pointer to the vector, which
  ! corresponds to iterat. The pointer is the return value of the
  ! function. If the vector is not InCore, the function terminates
  ! with an Error message. This function is used only in cases,
  ! where the vector for sure is stored in core (e.g. last entries
  ! in LList).

  integer(kind=iwp) :: LstPtr
  integer(kind=iwp), intent(in) :: iterat, iLList
  integer(kind=iwp) :: idum1, idum2, idum3, inode, ivptr

  LstPtr = -999999
  call GetNod(iterat,iLList,inode)
  if (inode == 0) then
    ! Hmmm, no entry found in LList, that's strange
    write(u6,*) 'LstPtr: inode <= 0'
    write(u6,*) 'inode=',inode
    call Abend()
  else if (InCore_f(inode)) then
    call InfNod(inode,idum1,idum2,ivptr,idum3)
    LstPtr = ivptr
  else
    ! Hmmm, no incore hit for this entry, that's strange
    write(u6,*) 'LstPtr: no incore hit for this entry'
    write(u6,*) 'inode=',inode
    call Abend()
  end if

end function LstPtr

subroutine KilLst(iLList)
  ! Free all memory of linked list LList

  use stdalloc, only: mma_deallocate

  integer(kind=iwp), intent(in) :: iLList
  integer(kind=iwp) :: iFlag, iroot

# ifdef _DEBUGPRINT_
  write(u6,*) 'KilLst'
  call StlLst(iLList)
# endif

  iroot = nLList(iLList,1)
  do while (iroot /= 0)
    iFlag = nLList(iroot,5)

    if (iFlag == 1) call mma_deallocate(SCF_V(iroot)%A)
    iroot = nLList(iroot,0)
  end do

end subroutine KilLst

subroutine DmpLst(iLList,LUnit,lDskPt)

  use stdalloc, only: mma_deallocate

  integer(kind=iwp), intent(in) :: iLList, LUnit
  integer(kind=iwp), intent(out) :: lDskPt
  integer(kind=iwp) :: iDskPt, iPtr1, iPtr2, iRoot, Length

  ! clear ErrCode
  nLList(iLList,0) = 0
  ! read listhead
  iroot = nLList(iLList,1)

  if (iroot <= 0) then
    ! linked list has zero length, that's strange
    ! save in any case listhead...
    lDskPt = 0
    iDskPt = lDskPt
    call iDaFile(LUnit,1,nLList(iLList,0),NodSiz,iDskPt)
    return
  end if

  do while (nLList(iroot,5) == 1)
    iPtr1 = iroot
    iPtr2 = iPtr1
    ! go either to last element or list or to last element in core
    ! and flush vector
    do while ((nLList(iPtr1,0) /= 0) .and. (nLList(iPtr1,5) == 1))
      iPtr2 = iPtr1
      iPtr1 = nLList(iPtr1,0)
    end do
    if (nLList(iPtr1,5) == 1) then
      ! nothing written on disk yet
      iDskPt = 0
      iPtr2 = iPtr1
      iPtr1 = nLList(iPtr1,1)
    else
      iDskPt = nLList(iPtr1,2)
      iPtr1 = nLList(iPtr2,1)
    end if
    ! now set up node again and write vec to disk
    nLList(iPtr2,1) = iDskPt
    nLList(iPtr2,5) = 0
    Length = nLList(iPtr2,3)

    call dDaFile(LUnit,1,SCF_V(iPtr2)%A,Length,iDskPt)
    call mma_deallocate(SCF_V(iPtr2)%A)
  end do
  lDskPt = iDskPt
  ! now all vectors are flushed... so dump linked list...
  call iDaFile(LUnit,1,nLList(iLList,0),NodSiz,iDskPt)

  do while (iroot /= 0)
    iPtr1 = iroot
    iroot = nLList(iroot,0)
    call iDaFile(LUnit,1,nLList(iPtr1,0),NodSiz,iDskPt)
  end do

end subroutine DmpLst

subroutine RclLst(iLList,LUnit,lDskPt,NoAllo)

  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(out) :: iLList
  integer(kind=iwp), intent(in) :: LUnit, NoAllo
  integer(kind=iwp), intent(inout) :: lDskPt
  integer(kind=iwp) :: incore, iPtr1, iPtr2, iRoot, lislen, lVec, MaxMem

  ! load listhead...
  lLList = lLList+1
  iLList = lLList

  call iDaFile(LUnit,2,nLList(iLList,0),NodSiz,lDskPt)
  iroot = nLList(iLList,1)

  if (iroot <= 0) then
    write(u6,*) 'RclLst: linked list has zero length, that''s strange!'
    ! linked list has zero length, that's strange
    !call Quit(20)
    return
  end if

  ! restore linked list...
  lislen = 1
  lLList = lLList+1
  iPtr1 = lLList

  nLList(iLList,1) = iPtr1
  call iDaFile(LUnit,2,nLList(iPtr1,0),NodSiz,lDskPt)

  iroot = iPtr1
  iPtr2 = iroot
  do while (nLList(iPtr2,0) /= 0)
    lislen = lislen+1
    lLList = lLList+1
    iPtr1 = lLList
    nLList(iPtr2,0) = iPtr1

    call iDaFile(LUnit,2,nLList(iPtr1,0),NodSiz,lDskPt)
    iPtr2 = iPtr1
  end do
  if (nLList(iLList,2) /= lislen) then
    write(u6,*) 'RclLst:LList length mismatch:',nLList(iLList,2),lislen
    call Abend()
  end if
  write(u6,*) 'Let''s restore...'
  ! now we have restored the list, let's fetch some vectors
  incore = nLList(iLList,3)
  call mma_maxDBLE(MaxMem)
  lvec = nLList(iroot,3)
  iPtr2 = iroot
  do while ((incore > 0) .and. (MaxMem-NoAllo >= lvec) .and. (iPtr2 > 0))
    lDskPt = nLList(iPtr2,1)

    if (iPtr2 > Maxnodes) then
      write(u6,*) 'iPtr2 > Maxnodes, restoring'
      call Abend()
    end if
    if (allocated(SCF_V(iPtr2)%A)) then
      write(u6,*) 'Node already allocated while restoring'
      write(u6,*) 'iPtr2=',iPtr2
      call Abend()
    end if
    call mma_Allocate(SCF_V(iPtr2)%A,lvec,Label='LVec')
    call dDaFile(LUnit,2,SCF_V(iPtr2)%A,lvec,lDskPt)
    nLList(iPtr2,1) = iPtr2
    nLList(iPtr2,2) = 0
    nLList(iPtr2,5) = 1
    iPtr2 = nLList(iPtr2,0)
    incore = incore-1

    call mma_maxDBLE(MaxMem)
  end do
  if (iPtr2 > 0) nLList(iLList,3) = nLList(iLList,3)-incore

end subroutine RclLst

subroutine StlLst(LLink)

  integer(kind=iwp), intent(in) :: LLink
  integer(kind=iwp) :: iRoot

  write(u6,*)
  write(u6,*) '*********** Status of Linked List *************'
  write(u6,*)
  write(u6,*) ' LLink:',LLink
  write(u6,*)
  write(u6,*) ' CNOD data'
  write(u6,*) 'Error code:                       ',nLList(LLink,0)
  write(u6,*) 'Pointer to first NODE in the list:',nLList(LLink,1)
  write(u6,*) 'Actual length of list:            ',nLList(LLink,2)
  write(u6,*) '# of vectors in core:             ',nLList(LLink,3)
  write(u6,*)
  iRoot = nLList(LLink,1)
  do while (iRoot /= 0)
    write(u6,*) ' NODE data'
    write(u6,*) 'NODE @:                         ',iRoot
    write(u6,*) 'Pointer to next NODE:           ',nLList(iRoot,0)
    write(u6,*) 'Pointer to stored vector:       ',nLList(iRoot,1)
    if (nLList(iRoot,5) >= 1) then
      write(u6,*) 'Vector status:                  in Core'
    else
      write(u6,*) 'Vector status:                  on Disk'
    end if
    write(u6,*) 'Next free position:             ',nLList(iRoot,2)
    write(u6,*) 'Length of vector:               ',nLList(iRoot,3)
    write(u6,*) 'Iteration number:               ',nLList(iRoot,4)
    write(u6,*)
    iRoot = nLList(iRoot,0)
  end do
  write(u6,*) '************ End of Status Report *************'
  write(u6,*)

end subroutine StlLst

end module LnkLst
