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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

module ipPage

use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
private

! n   : Length of CI-vector
! ida : disk address

integer(kind=iwp), parameter :: Max_CI_Vectors = 40, &
                                On_Disk = 0, In_Memory = 1, Null_Vector = 2, &
                                dWrite = 0, oWrite = 1, oRead = 2

integer(kind=iwp) :: ida(0:Max_CI_Vectors), iDisk_Addr_End = 0, Lu_ip, n(0:Max_CI_Vectors), n_CI_Vectors = 0, &
                     vStatus(0:Max_CI_Vectors)
logical(kind=iwp) :: DiskBased = .false.
type(Alloc1DArray_Type) :: W(0:Max_CI_Vectors)

public :: ipclose, ipget, ipin, ipin1, ipnout, ipopen, ipout, opout, W

contains

!=======================================================================
! Initiate the whole lot.
subroutine ipopen(page)

  logical(kind=iwp), intent(in) :: page
  integer(kind=iwp), external :: IsFreeUnit

  ! Ask how much memory is available

  if (Page) then

    ! Initiate for disk based storage.

    if (.not. DiskBased) then
      Lu_ip = IsFreeUnit(21)
      call Daname(Lu_ip,'TEMPCIV')
      DiskBased = .true.
    end if

    ! n  : Length of CI-vector
    ! ida: disk address

    n(0:Max_CI_Vectors) = 0
    ida(0:Max_CI_Vectors) = -1
    vStatus(0:Max_CI_Vectors) = Null_Vector

    ! iDisk_Addr_End: next free disk address
    ! n_CI_Vectors : number of CI-vectors

    iDisk_Addr_End = 0
    n_CI_Vectors = 0

  else

    call ipTerm()

  end if

end subroutine ipopen

!=======================================================================
! Object: release all vectors above and including the vector indexed ia.
subroutine ipclose(ia)

  integer(kind=iwp), intent(in) :: ia
  integer(kind=iwp) :: ii
  real(kind=wp) :: rdum(1)

  if (ia > Max_CI_Vectors) then
    write(u6,*) 'ipclose: ia > Max_CI_Vectors'
    write(u6,*) 'ia,Max_CI_Vectors=',ia,Max_CI_Vectors
    call Abend()
  end if

  ! Update iDisk_Addr_End

  iDisk_Addr_End = 0
  if (ia < 0) then

    n_CI_Vectors = 0
    call ipTerm()

  else

    n_CI_Vectors = ia-1
    if (DiskBased) then
      do ii=1,ia-1
        if (vStatus(ii) /= Null_Vector) call dDafile(Lu_ip,dWrite,rdum(1),n(ii),iDisk_Addr_End)
      end do
    end if

  end if

  ! Release memory and flag as a null vector

  do ii=max(ia,0),Max_CI_Vectors
    if (vStatus(ii) == In_Memory) then
      call mma_deallocate(W(ii)%A)
      ida(ii) = -1
      n(ii) = 0
      vStatus(ii) = Null_Vector
    end if
  end do

end subroutine ipclose

!=======================================================================
! Termination
subroutine ipterm()

  if (DiskBased) then
    call DaClos(Lu_ip)
    DiskBased = .false.
  end if

end subroutine ipterm

!=======================================================================
! Get the index of a vector with the length nn.
! Memory or disk space is allocated.
function ipget(nn)

  integer(kind=iwp) :: ipget
  integer(kind=iwp), intent(in) :: nn
  character(len=4) :: Label

  ! Take the next memory slot.

  n_CI_Vectors = n_CI_Vectors+1
  ipget = n_CI_Vectors

  if (n_CI_Vectors > Max_CI_Vectors) then
    write(u6,*) 'Number of CI vectors higher than Max_CI_Vectors'
    write(u6,*) 'Max_CI_Vectors=',Max_CI_Vectors
    call Abend()
  end if

  ida(ipget) = iDisk_Addr_End
  n(ipget) = nn

  ! Allocate memory for vector if of non-zero  length

  if (nn > 0) then
    write(Label,'(I3.3)') n_CI_Vectors
    call mma_allocate(W(ipget)%A,nn,Label='ipget'//Label)
    vStatus(ipget) = In_Memory
    W(ipget)%A(:) = Zero
  else
    !vStatus(ipget) = Null_Vector

    ! The calling code doesn't have the logic to handle the
    ! case that W(i)%A is not allocated. Hence, we have
    ! to make a dummy allocation to make sure that the compiler
    ! doesn't puke.
    n(ipget) = 1
    write(Label,'(I3.3)') n_CI_Vectors
    call mma_allocate(W(ipget)%A,1,Label='ipget'//Label)
    vStatus(ipget) = In_Memory
    W(ipget)%A(:) = Zero
  end if

  ! If DiskBased mode put vector on disc and release memory

  if (DiskBased) then
    if (vStatus(ipget) /= Null_Vector) then
      call dDafile(Lu_ip,oWrite,W(ipget)%A,nn,iDisk_Addr_End)
      vStatus(ipget) = On_Disk
      call mma_deallocate(W(ipget)%A)
    end if
  end if

end function ipget

! Object: retrieve vector ii with a length of nn and
!         make it available in memory as W(ii)%A
subroutine ipin1(ii,nn)

  integer(kind=iwp), intent(in) :: ii, nn
  integer(kind=iwp) :: idisk, nnn
  real(kind=wp), allocatable :: Tmp(:)

  if (ii > Max_CI_Vectors) then
    write(u6,*) 'ipin1: ii > Max_CI_Vectors'
    write(u6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
    call Abend()
  end if

  if (vStatus(ii) == In_Memory) then

    ! ii is in memory

    ! If the size of the vector is larger than was originally set
    ! resize the reservation and copy the content

    if (nn > n(ii)) then
      call mma_allocate(Tmp,nn,Label='Tmp')
      Tmp(:) = Zero
      Tmp(1:n(ii)) = W(ii)%A(:)
      call mma_deallocate(W(ii)%A)
      call mma_allocate(W(ii)%A,nn,Label='ipin1')
      W(ii)%A(:) = Tmp(:)
      call mma_deallocate(Tmp)
      n(ii) = nn
    end if

  else if (vStatus(ii) == On_Disk) then

    ! ii is on disk

    call mma_allocate(W(ii)%A,max(n(ii),nn),Label='ipin1')
    W(ii)%A(:) = Zero

    nnn = min(n(ii),nn)

    ! pick up from disk

    idisk = ida(ii)
    call dDafile(Lu_ip,oRead,W(ii)%A,nnn,idisk)
    vStatus(ii) = In_Memory

  else if (vStatus(ii) /= Null_Vector) then

    write(u6,*)
    write(u6,*) 'ipIn1: illegal vStatus(ii)'
    write(u6,*) 'ii=',ii
    write(u6,*)
    call Abend()

  end if

end subroutine ipin1

! Object: retrieve vector ii with a length of n(ii) and
!         make it available in memory as W(ii)%A
subroutine ipin(ii)

  integer(kind=iwp), intent(in) :: ii
  integer(kind=iwp) :: nn

  nn = n(ii)
  call ipin1(ii,nn)

end subroutine ipin

! ipout will page out vector ii to disk and free the memory area
subroutine ipout(ii)

  integer(kind=iwp), intent(in) :: ii
  integer(kind=iwp) :: idisk, nn

  if (DiskBased) then

    if ((vStatus(ii) == In_Memory) .and. (ii > 0)) then
      idisk = ida(ii)
      nn = n(ii)
      call dDafile(Lu_ip,oWrite,W(ii)%A,nn,idisk)
      vStatus(ii) = On_Disk
      call mma_deallocate(W(ii)%A)
    end if

  end if

end subroutine ipout

!=======================================================================
! Object: write all vectors in memory on disk but vector iii
subroutine ipnout(iii)

  integer(kind=iwp), intent(in) :: iii
  integer(kind=iwp) :: idisk, ii, nn

  if (iii > Max_CI_Vectors) then
    write(u6,*) 'ipout: iii > Max_CI_Vectors'
    write(u6,*) 'iii,Max_CI_Vectors=',iii,Max_CI_Vectors
    call Abend()
  end if

  if (DiskBased) then

    do ii=1,Max_CI_Vectors

      if ((vStatus(ii) == In_Memory) .and. (ii /= iii)) then
        idisk = ida(ii)
        nn = n(ii)
        call dDafile(Lu_ip,oWrite,W(ii)%A,nn,idisk)
        vStatus(ii) = On_Disk
        call mma_deallocate(W(ii)%A)
      end if

    end do

  end if

end subroutine ipnout

!=======================================================================
! opout will release the memory area of vector ii without updating the disk
subroutine opout(ii)

  integer(kind=iwp), intent(in) :: ii

  if (ii > Max_CI_Vectors) then
    write(u6,*) 'opout: ii > Max_CI_Vectors'
    write(u6,*) 'ii,Max_CI_Vectors=',ii,Max_CI_Vectors
    call Abend()
  end if

  if (DiskBased) then

    if ((vStatus(ii) == In_Memory) .and. (ii > 0)) then
      vStatus(ii) = On_Disk
      call mma_deallocate(W(ii)%A)
    end if

  end if

end subroutine opout

end module ipPage
