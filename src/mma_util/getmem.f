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

*  GetMem
*
*> @brief
*>   A simple work space manager, used by all programs in MOLCAS
*> @author Victor P. Vysotskiy
*>
*> @details
*> \p NameIn, \p KeyIn, and \p TypeIn are strings of any size. They are
*> not case sensitive, and only the four first letters matter.
*> If KeyIn is '``allo``' (or '``ALLO``' or ...) then ::GETMEM will return the
*> position of a previously unused piece of workspace, capable of holding
*> at least \p LENGTH items, and register that piece as being in use.
*> If \p TypeIn is '``Real``', the items will be accessible in
*> \c WORK(IPOS) ... ``WORK(IPOS-1+LENGTH)``.
*> If \p TypeIn is '``Inte``', the items will be accessible in
*> \c IWORK(IPOS) ... ``IWORK(IPOS-1+LENGTH)``.
*> If \p KeyIn is '``Free``', the piece will be returned to the free pool.
*> If \p KeyIn is '``Chec``', the boundaries of all allocated pieces will be
*> checked to see if the contents of guardian words, surrounding each
*> piece, are intact or not.
*> If \p KeyIn is '``List``', the allocated fields will be tabulated.
*> If \p KeyIn is '``Max``', there will be no allocation made, but the length
*> of the largest allocatable field is returned.
*> \p NameIn has no function, except that the user provides a label to the
*> field, which is used in error prints or listings.
*>
*> @note
*> An include file, WrkSpc.fh, declares common ``/WrkSpc/``,
*> containinig two arrays,
*> \c WORK and \c IWORK, which  are equivalenced.
*> ::GETMEM uses calls to the Molcas's MA memory allocator routines.
*>
*> @param[in]     NameIn Arbitrary label
*> @param[in]     KeyIn  ``Allo`` / ``Free`` / ``Check`` / ``List`` / ``Max``
*> @param[in]     TypeIn ``Real`` / ``Inte`` / ``Char``
*> @param[in,out] iPos   Position
*> @param[in,out] Length Nr of items
************************************************************************
      Subroutine GetMem (NameIn,KeyIn,TypeIn,iPos,Length)
************************************************************************
*                                                                      *
* History: Victor P. Vysotskiy                                         *
*    2012: Native Molcas's Memory Allocator; Thread safety             *
*                                                                      *
************************************************************************

#include "SysCtl.fh"
#include "warnings.h"
#include "WrkSpc.fh"
#include "mama.fh"
*
*
      Character*(*) NameIn,KeyIn,TypeIn
      Character*8   FldNam,eopr,eoprcc,elbl,etyp
      Character*4   Key,VarTyp
#ifdef _GARBLE_
      Character*5   xKey
      Logical       SkipGarble
#endif
      Interface
        Function c_getmem(name_,Op,dtyp,offset,len_)
     &           bind(C,name='c_getmem_')
          Use, Intrinsic :: iso_c_binding, only: c_char
          Use Definitions, only: MOLCAS_C_INT
          Integer(kind=MOLCAS_C_INT) :: c_getmem
          Character(kind=c_char) :: name_(*), Op(*), dtyp(*)
          Integer(kind=MOLCAS_C_INT) :: offset, len_
        End Function c_getmem
      End Interface


*----------------------------------------------------------------------*
*     Initialize the Common / MemCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      If ( MemCtl(ipStat).ne.ON ) then
         Call IniMem()
      End if
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
      eoprcc='CHECK'
      eoprcc(8:8)=char(0)

*----------------------------------------------------------------------*
*     Trace memory                                                     *
*----------------------------------------------------------------------*
      If (MemCtl(ipCheck).eq.ON .or. MemCtl(ipTrace).eq.ON) Then
         iRc=c_getmem(elbl,eoprcc,etyp,ip_iDummy,ip_iDummy)
      End If
#ifdef _GARBLE_
*----------------------------------------------------------------------*
*     Skip garble                                                      *
*----------------------------------------------------------------------*
      Call StdFmt(KeyIn,xKey)
      If ((xKey.eq.'ALLON').or.(xKey.eq.'RGSTN')) Then
        SkipGarble = .True.
      Else
        SkipGarble = .False.
      End If
#endif
*----------------------------------------------------------------------*
*     Allocate new memory                                              *
*----------------------------------------------------------------------*

      If(Key.ne.'ALLO') iPos=iPos-kind2goff(VarTyp)
      iRc=c_getmem(elbl,eopr,etyp,iPos,Length)
      If(iRc.lt.0) Then
        If ( Key.eq.'ALLO' ) Then
          Write (6,'(A)')
     &         'MMA failed to allocate a memory block.'
        Else If ( Key.eq.'FREE' ) Then
          Write (6,'(A)')
     &         'MMA failed to release the memory block for further use.'
          iRc=c_getmem(elbl,eoprcc,etyp,ip_iDummy,ip_iDummy)
        Else
          Write (6,*)
        End If
        Go To 777
      End If

      If ( Key.eq.'ALLO' .or. Key.eq.'LENG'.or.
     > Key.eq.'FLUS'. or. Key.eq.'MAX' .or.
     > Key.eq.'CHEC'.or.Key.eq.'LIST' .or.
     > Key.eq.'RGST') Then
         iPos=iPos+kind2goff(VarTyp)
*----------------------------------------------------------------------*
*     Release a memory block or return length or decrease length       *
*----------------------------------------------------------------------*
      End If

#ifdef _GARBLE_
      If ( Key.eq.'ALLO' .or. Key.eq.'RGST') Then
        If (.not. SkipGarble) Call Garble(iPos,Length,VarTyp)
      End If
#endif

      Return
*
 777  Continue
      Call Quit(_RC_MEMORY_ERROR_)

      End




      function kind2goff(var)
#include "mama.fh"
      character*(4) var
      kind2goff=0
      if(var.eq.'INTE') kind2goff=iofint
      if(var.eq.'REAL') kind2goff=iofdbl
      if(var.eq.'CHAR') kind2goff=iofchr
      return
      end

#ifdef _GARBLE_
      subroutine garble(ipos,length,vartyp)
      use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
      implicit none
*include "SysDef.fh"
#include "WrkSpc.fh"
      integer :: ipos, length
      character(len=*) :: vartyp
      real*8, parameter ::    dgarbage(1) = [huge(1.0d0)]
      integer, parameter ::   igarbage(1) = [huge(1)]
      integer*1, parameter :: i1garbage = huge(i1garbage)

      select case(vartyp)
      case ('REAL')
        call dcopy_(length,dgarbage,0,work(ipos),1)
      case ('INTE')
        call icopy(length,igarbage,0,iwork(ipos),1)
      case ('CHAR')
        call garble_char(Work)
      end select

      ! This is to allow type punning without an explicit interface
      contains

      subroutine garble_char(buf)
      use Definitions, only: RtoB
      real*8, target :: buf(*)
      integer*1, pointer :: ibuf(:)
      integer :: ioff1, ioff2, foff1
      ioff1 = (ipos-1)/RtoB+1
      ioff2 = mod(ipos-1,RtoB)+1
      foff1 = (ipos+length-2)/RtoB+1
      call c_f_pointer(c_loc(buf(ioff1)),ibuf,[(foff1-ioff1+1)*RtoB])
      ibuf(ioff2:ioff2+length-1) = i1garbage
      nullify(ibuf)
      end subroutine garble_char

      end subroutine
#endif
