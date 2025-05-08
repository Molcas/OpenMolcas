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
! Copyright (C) 2012, Victor P. Vysotskiy                              *
!               2025, Ignacio Fdez. Galvan                             *
!***********************************************************************

!  GetMem
!
!> @brief
!>   A simple work space manager, used by all programs in MOLCAS
!> @author Victor P. Vysotskiy
!>
!> @details
!> \p NameIn, \p KeyIn, and \p TypeIn are strings of any size. They are
!> not case sensitive, and only the four first letters matter.
!> If \p KeyIn is '``allo``' (or '``ALLO``' or ...) then ::GETMEM will return the
!> position of a previously unused piece of workspace, capable of holding
!> at least \p LENGTH items, and register that piece as being in use.
!> If \p TypeIn is '``Real``', the items will be accessible in
!> \c WORK(IPOS) ... ``WORK(IPOS-1+LENGTH)``.
!> If \p TypeIn is '``Inte``', the items will be accessible in
!> \c IWORK(IPOS) ... ``IWORK(IPOS-1+LENGTH)``.
!> If \p KeyIn is '``Free``', the piece will be returned to the free pool.
!> If \p KeyIn is '``Chec``', the boundaries of all allocated pieces will be
!> checked to see if the contents of guardian words, surrounding each
!> piece, are intact or not.
!> If \p KeyIn is '``List``', the allocated fields will be tabulated.
!> If \p KeyIn is '``Max``', there will be no allocation made, but the length
!> of the largest allocatable field is returned.
!> \p NameIn has no function, except that the user provides a label to the
!> field, which is used in error prints or listings.
!>
!> @note
!> An include file, WrkSpc.fh, declares common ``/WrkSpc/``,
!> containing two arrays,
!> \c WORK and \c IWORK, which  are equivalenced.
!> ::GETMEM uses calls to the Molcas's MA memory allocator routines.
!>
!> @param[in]     NameIn Arbitrary label
!> @param[in]     KeyIn  ``Allo`` / ``Free`` / ``Check`` / ``List`` / ``Max``
!> @param[in]     TypeIn ``Real`` / ``Inte`` / ``Char``
!> @param[in,out] iPos   Position
!> @param[in,out] Length Nr of items
!***********************************************************************

subroutine GetMem(NameIn,KeyIn,TypeIn,iPos,Length)
!***********************************************************************
!                                                                      *
! History: Victor P. Vysotskiy                                         *
!    2012: Native Molcas's Memory Allocator; Thread safety             *
!          Ignacio Fdez. Galvan                                        *
!    2025: Garble using C pointers                                     *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: NameIn, KeyIn, TypeIn
integer(kind=iwp), intent(inout) :: iPos
integer(kind=iwp), intent(in) :: Length
#include "warnings.h"
#include "WrkSpc.fh"
integer(kind=iwp) :: irc
character(len=8) :: elbl, eopr, etyp, FldNam
character(len=4) :: Key, VarTyp
#ifdef _GARBLE_
logical(kind=iwp) :: SkipGarble
character(len=5) :: xKey
#endif
interface
  function c_getmem(name_,Op,dtyp,offset,len_) bind(C,name='c_getmem_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_getmem
    character(kind=c_char), intent(in) :: name_(*), Op(*), dtyp(*)
    integer(kind=MOLCAS_C_INT), intent(inout) :: offset
    integer(kind=MOLCAS_C_INT), intent(in) :: len_
  end function c_getmem
end interface

!----------------------------------------------------------------------*
!     Initialize the Common / MemCtl / the first time it is referenced *
!----------------------------------------------------------------------*
if (.not. MemStat) call IniMem()
!----------------------------------------------------------------------*
!     convert input strings to standard format                         *
!----------------------------------------------------------------------*
call StdFmt(NameIn,FldNam)
call StdFmt(KeyIn,Key)
call StdFmt(TypeIn,VarTyp)
!----------------------------------------------------------------------*
!     prepare passed values to the C char format                       *
!----------------------------------------------------------------------*
elbl = FldNam
elbl(8:8) = char(0)
eopr = Key
eopr(8:8) = char(0)
etyp = VarTyp
etyp(8:8) = char(0)

#ifdef _GARBLE_
!----------------------------------------------------------------------*
!     Skip garble                                                      *
!----------------------------------------------------------------------*
call StdFmt(KeyIn,xKey)
if ((xKey == 'ALLON') .or. (xKey == 'RGSTN')) then
  SkipGarble = .true.
else
  SkipGarble = .false.
end if
#endif

!----------------------------------------------------------------------*
!     Allocate new memory                                              *
!----------------------------------------------------------------------*
if (Key /= 'ALLO') iPos = iPos-1
iRc = c_getmem(elbl,eopr,etyp,iPos,Length)
if (iRc < 0) then
  if (Key == 'ALLO') then
    write(u6,'(A)') 'MMA failed to allocate a memory block.'
  else if (Key == 'FREE') then
    write(u6,'(A)') 'MMA failed to release the memory block for further use.'
    call abend()
  else
    write(u6,*)
  end if
  call Quit(_RC_MEMORY_ERROR_)
end if

if ((Key == 'ALLO') .or. (Key == 'LENG') .or. (Key == 'FLUS') .or. (Key == 'MAX') .or. (Key == 'CHEC') .or. (Key == 'LIST') .or. &
    (Key == 'RGST')) iPos = iPos+1

#ifdef _GARBLE_
if ((Key == 'ALLO') .or. (Key == 'RGST')) then
  if (.not. SkipGarble) call Garble(iPos,Length,VarTyp)
end if
#endif

end subroutine GetMem
