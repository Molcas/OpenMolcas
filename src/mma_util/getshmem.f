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
c
      Subroutine GetShMem (NameIn,KeyIn,TypeIn,iPos,Length,Path,ShmId)
************************************************************
*
*   <DOC>
*     <Name>GetMem</Name>
*     <Syntax>Call GetShMem(NameIn,KeyIn,TypeIn,iPos,Length)</Syntax>
*     <Arguments>
*       \Argument{NameIn}{Arbitrary label}{Char*(*)}{in}
*       \Argument{KeyIn}{Allo $|$ Free }{Char*(*)}{in}
*       \Argument{TypeIn}{Real $|$ Inte $|$ Char $|$ Sngl}{Char*(*)}{in}
*       \Argument{iPos}{Position}{Integer}{inout}
*       \Argument{Length}{Nr of items}{Integer}{inout}
*       \Argument{Path}{an arbitrary path or empty }{Char*(*)}{in}
*       \Argument{ShmId}{An unique shared memory id}{Integer}{inout}
*     </Arguments>
*     <Purpose>
* A simple work space manager, used by all programs in MOLCAS.
*     </Purpose>
*     <Dependencies>
* An include file, WrkSpc.fh, declares a commons /WrkSpc/,
* /cWrkSpc/. The first common contains three arrays,
*  WORK(), SWORK(*) and IWORK(*), which  are equivalenced. The vector, CWORK,
* belongs to the second commons.
* GETMEM uses calls to the Molcas's MA memory allocator routine.
*     </Dependencies>
*     <Author></Author> Victor P. Vysotskiy
*     <Modified_by></Modified_by>
*     <Description>
* NameIn, KeyIn, and TypeIn are strings of any size. They are
* not case sensitive, and only the four first letters matter.
* If KeyIn is 'allo' (or 'ALLO' or ...) then GETMEM will return the
* position of a previously unused piece of shared memory workspace,
* capable of holding  at least LENGTH items, and register that piece as being in use.
* If TypeIn is 'Real', the items will be accessible in
* WORK(IPOS)..WORK(IPOS-1+LENGTH).
* If TypeIn is 'Inte', the items will be accessible in
* IWORK(IPOS)..IWORK(IPOS-1+LENGTH).
* If TypeIn is 'Sngl', the items will be accessible in
* SWORK(IPOS)..SWORK(IPOS-1+LENGTH).
* If KeyIn is 'Free', the piece will be returned to the free pool.
* NameIn has no function, except that the user provides a label to the
* field, which is used in error prints or listings.
* Path is functional only in conjuction with the 'Files In Memory' I/O
* layer. Otherwise, just simple use an empty string.
* ShmId on enter must be an non-zero integer value which is being replaced on output
* by an unique integer ID. This ID must be used for any further manipulation
* on shared memory segment with the 'iPos' offset.
*     </Description>
*    </DOC>
*----------------------------------------------------------------------*
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
