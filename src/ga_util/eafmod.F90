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
! Copyright (C) Martin Schuetz                                         *
!               Roland Lindh                                           *
!               2015, Steven Vancoillie                                *
!               2017,2026, Ignacio Fdez. Galvan                        *
!***********************************************************************

!***********************************************************************
! Wrapper routines against the GA-ChemIO or MPI-IO package.            *
! (EAF = exclusive access file)                                        *
!                                                                      *
! If not an mpp installation the AIX-IO facility will be used.         *
!                                                                      *
! dEAFxxxx subroutines simply provides an interface for the EAF        *
! subroutines with real buffers instead of integers, note that nBuf is *
! still the length of the integer buffer.                              *
!***********************************************************************

module EAFmod

#ifdef ADD_
#define eaf_open eaf_open_
#define eaf_close eaf_close_
#define eaf_awrite eaf_awrite_
#define eaf_aread eaf_aread_
#define eaf_write eaf_write_
#define eaf_read eaf_read_
#define eaf_wait eaf_wait_
#endif

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use Definitions, only: u6, ItoB
#endif
use Definitions, only: iwp, wp

implicit none
private

#ifdef _MOLCAS_MPP_
#ifdef _HAVE_EXTRA_
integer(kind=iwp), external :: molcas_eaf_aread, molcas_eaf_awrite, molcas_eaf_close, molcas_eaf_open, molcas_eaf_read, &
                               molcas_eaf_wait molcas_eaf_write
#else
! Include file from GlobalArrays
#include "eaf.fh"
#endif
#endif

public :: dEAFARead, dEAFAWrite, dEAFRead, dEAFWrite, EAFARead, EAFAWrite, EAFClose, EAFOpen, EAFRead, EAFWait, EAFWrite

contains

subroutine EAFOpen(Lu,FName)

# ifdef _MOLCAS_MPP_
  use Para_Info, only: MyRank
# endif

  integer(kind=iwp) :: Lu
  character(len=*) :: FName
  integer(kind=iwp), external :: IsFreeUnit
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
  character(len=200) :: FN
# ifdef _HAVE_EXTRA_
  integer(kind=iwp) :: n
# endif
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    write(FN,'(A,"_",I0)') trim(FName),MyRank
#   ifdef _HAVE_EXTRA_
    n = len_trim(FName)
    FN(n+1:n+1) = char(0)
    iRC = molcas_eaf_open(FN,Lu)
#   else
    iRC = eaf_open(FN,eaf_rw,Lu)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFOpen: Abort!'
      write(u6,*) 'iRC=',iRC
      write(u6,*) 'Lu=',Lu
      write(u6,*) 'Fname=',FName
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  else
# endif
    Lu = IsFreeUnit(7)
    call DaName_MF(Lu,FName)
# ifdef _MOLCAS_MPP_
  end if
# endif
end subroutine EAFOpen

subroutine EAFClose(Lu)

  integer(kind=iwp) :: Lu
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
#   ifdef _HAVE_EXTRA_
    iRC = molcas_eaf_close(Lu)
#   else
    iRC = eaf_close(Lu)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFClose: Abort!'
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  else
# endif
    call DaClos(Lu)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine EAFClose

subroutine EAFAWrite(Lu,Buf,nBuf,Disk,id)

  integer(kind=iwp) :: Lu, nBuf, Buf(nBuf), id
  real(kind=wp) :: Disk
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
#   ifdef _HAVE_EXTRA_
    iRC = molcas_eaf_awrite(Lu,Disk,Buf,nBuf*ItoB,id)
#   else
    iRC = eaf_awrite(Lu,Disk,Buf,nBuf*ItoB,id)
    Disk = Disk+real(nBuf*ItoB,kind=wp)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFAWrite: Abort!'
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  else
# endif
    id = 0
    call EAFWrite(Lu,Buf,nBuf,Disk)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine EAFAWrite

subroutine EAFARead(Lu,Buf,nBuf,Disk,id)

  integer(kind=iwp) :: Lu, nBuf, Buf(nBuf), id
  real(kind=wp) :: Disk
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
#   ifdef _HAVE_EXTRA_
    iRC = molcas_eaf_aread(Lu,Disk,Buf,nBuf*ItoB,id)
#   else
    iRC = eaf_aread(Lu,Disk,Buf,nBuf*ItoB,id)
    Disk = Disk+real(nBuf*ItoB,kind=wp)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFARead: Abort!'
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  else
# endif
    id = 0
    call EAFRead(Lu,Buf,nBuf,Disk)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine EAFARead

subroutine EAFWrite(Lu,Buf,nBuf,Disk)

  integer(kind=iwp) :: Lu, nBuf, Buf(nBuf)
  real(kind=wp) :: Disk
  integer(kind=iwp) :: iDisk
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
#   ifdef _HAVE_EXTRA_
    iRC = molcas_eaf_write(Lu,Disk,Buf,nBuf*ItoB)
#   else
    iRC = eaf_write(Lu,Disk,Buf,nBuf*ItoB)
    Disk = Disk+real(nBuf*ItoB,kind=wp)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFWrite: Abort!'
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  else
# endif
    iDisk = int(Disk)
    call iDaFile(Lu,1,Buf,nBuf,iDisk)
    Disk = real(iDisk,kind=wp)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine EAFWrite

subroutine EAFRead(Lu,Buf,nBuf,Disk)

  integer(kind=iwp) :: Lu, nBuf, Buf(nBuf)
  real(kind=wp) :: Disk
  integer(kind=iwp) :: iDisk
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
#   ifdef _HAVE_EXTRA_
    iRC = molcas_eaf_read(Lu,Disk,Buf,nBuf*ItoB)
#   else
    iRC = eaf_read(Lu,Disk,Buf,nBuf*ItoB)
    Disk = Disk+real(nBuf*ItoB,kind=wp)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFRead: Abort!'
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  else
# endif
    iDisk = int(Disk)
    call iDaFile(Lu,2,Buf,nBuf,iDisk)
    Disk = real(iDisk,kind=wp)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine EAFRead

subroutine EAFWait(Lu,id)

  integer(kind=iwp) :: Lu, id
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: iRC
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
#   ifdef _HAVE_EXTRA_
    iRC = molcas_eaf_wait(Lu,id)
#   else
    iRC = eaf_wait(Lu,id)
#   endif
    if (iRC /= 0) then
      write(u6,*) 'EAFWait: Abort!'
      write(u6,*) 'EAF_Err_Code =',iRC
      call Abend()
    end if
  end if
# else
# include "macros.fh"
  unused_var(Lu)
  unused_var(id)
# endif

end subroutine EAFWait

subroutine dEAFARead(Lu,Buf,nBuf,Disk,id)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(out) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), intent(out) :: id
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFARead(Lu,iBuf,nBuf,Disk,id)
  nullify(iBuf)

end subroutine dEAFARead

subroutine dEAFAWrite(Lu,Buf,nBuf,Disk,id)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(in) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), intent(out) :: id
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFAWrite(Lu,iBuf,nBuf,Disk,id)
  nullify(iBuf)

end subroutine dEAFAWrite

subroutine dEAFRead(Lu,Buf,nBuf,Disk)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(out) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFRead(Lu,iBuf,nBuf,Disk)
  nullify(iBuf)

end subroutine dEAFRead

subroutine dEAFWrite(Lu,Buf,nBuf,Disk)

  integer(kind=iwp), intent(in) :: Lu, nBuf
  real(kind=wp), target, intent(in) :: Buf(nBuf)
  real(kind=wp), intent(inout) :: Disk
  integer(kind=iwp), pointer :: iBuf(:)

  call c_f_pointer(c_loc(Buf),iBuf,[nBuf])
  call EAFWrite(Lu,iBuf,nBuf,Disk)
  nullify(iBuf)

end subroutine dEAFWrite

end module EAFmod
