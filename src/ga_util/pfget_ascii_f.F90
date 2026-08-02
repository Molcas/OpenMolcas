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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

#include "compiler_features.h"
#ifndef _HAVE_EXTRA_

! Broadcast a file from the master to the slaves
subroutine PFGet_ASCII(FName)

#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_rootid, King
use Definitions, only: iwp, u6, ItoB
#endif

implicit none
character(len=*), intent(in) :: FName
#ifdef _MOLCAS_MPP_
integer(kind=iwp), parameter :: LBuf = 4096
integer(kind=iwp) :: Err, FLen, LU, Num, Pos
logical(kind=iwp) :: Failed, Found
character(len=LBuf) :: Buf
integer(kind=iwp), external :: IsFreeUnit
#include "mafdecls.fh"
interface
  subroutine GA_Brdcst(tp,buf,lenbuf,root)
    import :: iwp
    integer(kind=iwp) :: tp, lenbuf, root
    type(*) :: buf
  end subroutine GA_Brdcst
end interface

! Note that each process opens only one file, so there is a single
! unit number LU
LU = IsFreeUnit(10)
! Check file existence and read size on the master
if (King()) then
  call f_Inquire(FName,Found)
  if (Found) then
    call Molcas_Open_Ext2(LU,FName,"stream","unformatted",Err,.false.,0,"old",Failed)
    if (Failed .or. (Err /= 0)) then
      write(u6,*) "Failed to open file ",trim(FName)
      call AbEnd()
    end if
    inquire(LU,Size=FLen)
  else
    FLen = 0
  end if
end if
! Broadcast the file size
Err = 0
call GA_Brdcst(MT_INT,FLen,1*ItoB,mpp_rootid)
if (FLen <= 0) return
! Open file for writing in the slaves
if (.not. King()) then
  call Molcas_Open_Ext2(LU,FName,"stream","unformatted",Err,.false.,0,"replace",Failed)
  if (Failed .or. (Err /= 0)) then
    write(u6,*) "Failed to open file ",trim(FName)
    call AbEnd()
  end if
end if
! Pass the file content in chunks
Pos = 0
do while (Pos < FLen)
  ! Length of this chunk
  Num = min(LBuf,FLen-Pos)
  ! The master reads the file
  if (King()) then
    read(LU,IOStat=Err) Buf(1:Num)
    if (Err /= 0) then
      write(u6,*) "Error reading the file ",trim(FName)
      call AbEnd()
    end if
  end if
  call GA_Brdcst(MT_BYTE,Buf,Num,mpp_rootid)
  ! The slaves write the file
  if (.not. King()) then
    write(LU,IOStat=Err) Buf(1:Num)
    if (Err /= 0) then
      write(u6,*) "Error writing the file ",trim(FName)
      call AbEnd()
    end if
  end if
  Pos = Pos+LBuf
end do
close(LU)
#else
#include "macros.fh"
unused_var(_str(FName))
#endif

end subroutine PFGet_ASCII

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(PFGet_ASCII)

#endif
