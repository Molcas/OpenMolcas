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
* Copyright (C) 2012, Victor P. Vysotskiy                              *
************************************************************************

*  GetShMem
*
*> @brief
*>   A simple work space manager, used by all programs in MOLCAS
*> @author  Victor P. Vysotskiy
*>
*> @details
*> \p NameIn, \p KeyIn, and \p TypeIn are strings of any size. They are
*> not case sensitive, and only the four first letters matter.
*> If \p KeyIn is '``allo``' (or '``ALLO``' or ...) then ::GETMEM will return the
*> position of a previously unused piece of shared memory workspace,
*> capable of holding at least \p LENGTH items, and register that piece as being in use.
*> If \p TypeIn is '``Real``', the items will be accessible in
*> \c WORK(IPOS) ... ``WORK(IPOS-1+LENGTH)``.
*> If \p TypeIn is '``Inte``', the items will be accessible in
*> \c IWORK(IPOS) ... ``IWORK(IPOS-1+LENGTH)``.
*> If \p TypeIn is '``Sngl``', the items will be accessible in
*> \c SWORK(IPOS) ... ``SWORK(IPOS-1+LENGTH)``.
*> If \p KeyIn is '``Free``', the piece will be returned to the free pool.
*> \p NameIn has no function, except that the user provides a label to the
*> field, which is used in error prints or listings.
*> \p Path is functional only in conjuction with the 'Files In Memory' I/O
*> layer. Otherwise, just simple use an empty string.
*> \p ShmId on enter must be a non-zero integer value which is being replaced on output
*> by an unique integer ID. This ID must be used for any further manipulation
*> on shared memory segment with the \p iPos offset.
*>
*> @note
*> An include file, WrkSpc.fh, declares commons ``/WrkSpc/`` and
*> ``/cWrkSpc/``. The first common contains three arrays,
*> \c WORK, \c SWORK and \c IWORK, which  are equivalenced. The vector, \c CWORK,
*> belongs to the second common.
*> ::GETMEM uses calls to the Molcas's MA memory allocator routines.
*>
*> @param[in]     NameIn Arbitrary label
*> @param[in]     KeyIn Allo $|$ Free
*> @param[in]     TypeIn Real $|$ Inte $|$ Char $|$ Sngl
*> @param[in,out] iPos   Position
*> @param[in,out] Length Nr of items
*> @param[in]     Path   An arbitrary path or empty
*> @param[in,out] ShmId  An unique shared memory id
************************************************************************
      Subroutine GetShMem (NameIn,KeyIn,TypeIn,iPos,Length,Path,ShmId)
************************************************************************
*                                                                      *
* History: Victor P. Vysotskiy                                         *
*    2012: Native Molcas's Memory Allocator; Thread safety             *
*                                                                      *
************************************************************************
#include "SysCtl.fh"
#include "warnings.fh"
#include "WrkSpc.fh"
#include "mama.fh"
#ifdef _OPENMP
      Include 'omp_lib.h'
#endif
*
*
      Character*(*)    NameIn,KeyIn,TypeIn
      Character*(*)    Path
      Character*4096   Cpath
      Character*8      FldNam,eopr,elbl,etyp
      Character*4      Key,VarTyp
      Integer          iPos,Length,ShmId
      Integer          c_getshmem
      External         c_getshmem

*----------------------------------------------------------------------*
*     Initialize the Common / MemCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      If ( MemCtl(ipStat).ne.ON ) then
         Call IniMem()
      End if
      If ( MemCtl(ipQuery).eq.ON ) Call qEnter('GetShMem')
*----------------------------------------------------------------------*
*     read default parameters from Common / MemCtl /                   *
*----------------------------------------------------------------------*
      iW=MemCtl(ipSysOut)
      If ( MemCtl(ipTrace).eq.ON ) then
         Write(iW,*) ' <<< Entering GetMem 5.0 >>>'
         Write(iW,'(A,2X,A4)') ' Clear  =      ',MemCtl(ipClear)
         Write(iW,'(A,2X,A4)') ' Key    =    ',KeyIn
         Write(iW,'(A,2X,A4)') ' Name   =    ',NameIn
         Write(iW,'(A,2X,A4)') ' Type   =    ',TypeIn
         Write(iW,'(A,I12)') ' length =    ',Length
         Write(iW,'(A,I12)') ' iPos   =    ',iPos
      End If
*----------------------------------------------------------------------*
*     convert input strings to standard format                         *
*----------------------------------------------------------------------*
      Call StdFmt(NameIn,FldNam)
      Call StdFmt(KeyIn,Key)
      Call StdFmt(TypeIn,VarTyp)
*----------------------------------------------------------------------*
*     prepare passed values to the C char format                       *
*----------------------------------------------------------------------*
      elbl=FldNam
      elbl(8:8)=char(0)
      eopr=Key
      eopr(8:8)=char(0)
      etyp=VarTyp
      etyp(8:8)=char(0)
C      n=LEN(TRIM(path)); cpath(n+1:n+1)=char(0)
      n=LEN(path(1:index(path,' ')-1))
      cpath(n+1:n+1)=char(0)
*----------------------------------------------------------------------*
*     Allocate/Free new memory                                              *
*----------------------------------------------------------------------*

      If( Key.ne.'ALLO') iPos=iPos-kind2goff(VarTyp)
      iRc=c_getshmem(elbl,eopr,etyp,iPos,Length,cpath,ShmId)
      If(iRc.lt.0) Then
        If ( Key.eq.'ALLO' ) Then
          Write (6,'(A)')
     &         'MMA failed to allocate a memory block.'
        Else If ( Key.eq.'FREE' ) Then
          Write (6,'(A)')
     &         'MMA failed to release the memory block for further use.'
        End If
        Go To 777
      End If

      If ( Key.eq.'ALLO') iPos=iPos+kind2goff(VarTyp)

      If ( MemCtl(ipQuery).eq.ON ) Call qExit('GetShMem')
      Return
*
 777  Continue
      Call QTrace()
      Call Quit(_RC_MEMORY_ERROR_)

      End
