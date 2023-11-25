!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

#include "compiler_features.h"
#ifdef _GROMACS_

subroutine Fetch_QMMM(CastMM,nCastMM)

use, intrinsic :: iso_c_binding, only: c_char, c_int, c_loc, c_ptr
use espf_global, only: MMI, MMO, QM, TPRDefName
use Isotopes, only: PTab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten, Angstrom
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCastMM
integer(kind=iwp), intent(in) :: CastMM(nCastMM)
#include "LenIn.fh"
integer(kind=iwp) :: iAtGMX, iAtNmbGMX, iAtOut, iCastMM, iFirst, iGrpGMX, iLast, iOk, LuXYZ, nAtGMX, nAtIn, nAtOut
logical(kind=iwp) :: Exists
character(len=LenIn) :: Symbol
character(len=256) :: LogFileName, Message, TPRFileName
type(c_ptr) :: ipCR, ipGMS
integer(kind=iwp), allocatable :: AT(:)
real(kind=wp), allocatable :: CoordGMX(:,:), CoordMMO(:,:)
character(len=LenIn), allocatable :: LabMMO(:)
integer(kind=iwp), parameter :: QMGMX = 0, MMGMX = 1
real(kind=wp), parameter :: AuToNm = Angstrom/Ten, NmToAng = Ten
integer(kind=iwp), external :: isFreeUnit
interface
  subroutine mmslave_done(gms) bind(C,NAME='mmslave_done_')
    import :: c_ptr
    type(c_ptr), value :: gms
  end subroutine mmslave_done
  function mmslave_init(cr,log_) bind(C,NAME='mmslave_init_')
    import :: c_char, c_ptr
    type(c_ptr) :: mmslave_init
    type(c_ptr), value :: cr
    character(kind=c_char) :: log_(*)
  end function mmslave_init
  function mmslave_get_atomnumber(gms,id) bind(C,NAME='mmslave_get_atomnumber_')
    import :: c_int, c_ptr
    integer(kind=c_int) :: mmslave_get_atomnumber
    type(c_ptr), value :: gms
    integer(kind=c_int), value :: id
  end function mmslave_get_atomnumber
  function mmslave_get_group_id(gms,id) bind(C,NAME='mmslave_get_group_id_')
    import :: c_int, c_ptr
    integer(kind=c_int) :: mmslave_get_group_id
    type(c_ptr), value :: gms
    integer(kind=c_int), value :: id
  end function mmslave_get_group_id
  function mmslave_natoms(gms) bind(C,NAME='mmslave_natoms_')
    import :: c_int, c_ptr
    integer(kind=c_int) :: mmslave_natoms
    type(c_ptr), value :: gms
  end function mmslave_natoms
  function mmslave_read_tpr(tpr,gms) bind(C,NAME='mmslave_read_tpr_')
    import :: c_char, c_int, c_ptr
    integer(kind=c_int) :: mmslave_read_tpr
    character(kind=c_char) :: tpr(*)
    type(c_ptr), value :: gms
  end function mmslave_read_tpr
  function init_commrec() bind(C,NAME='init_commrec_')
    import :: c_ptr
    type(c_ptr) :: init_commrec
  end function init_commrec
end interface

! Initialize Gromacs mmslave
ipCR = init_commrec()
call prgmtranslate('GMX.LOG',LogFileName,iLast)
LogFileName(iLast+1:iLast+1) = char(0)
ipGMS = mmslave_init(ipCR,LogFileName)

! Tell Gromacs to read tpr file
TPRFileName = TPRDefName
iLast = len_trim(TPRFileName)
call f_inquire(TPRFileName,Exists)
if (.not. Exists) then
  Message = 'File '//TPRFileName(1:iLast)//' not found'
  call WarningMessage(2,Message)
  call Quit_OnUserError()
end if
TPRFileName(iLast+1:iLast+1) = char(0)
iOk = mmslave_read_tpr(TPRFileName,ipGMS)
if (iOk /= 1) then
  Message = 'Error reading tpr file'
  call WarningMessage(2,Message)
  call Quit_OnUserError()
end if

! Fetch coordinates from Gromacs
nAtGMX = mmslave_natoms(ipGMS)
call mma_allocate(CoordGMX,3,nAtGMX)
CoordGMX(:,:) = Zero
iOk = mmslave_copyX_wrapper(ipGMS,int(nAtGMX,kind=c_int),CoordGMX)
if (iOk /= 1) then
  Message = 'Fetch_QMMM: mmslave_copyx is not ok'
  call WarningMessage(2,Message)
  call Abend()
end if

! Find out and count atom types
nAtIn = 0
nAtOut = 0
call mma_allocate(AT,nAtGMX)
do iAtGMX=1,nAtGMX
  iGrpGMX = mmslave_get_group_id(ipGMS,int(iAtGMX-1,kind=c_int))
  if (iGrpGMX == QMGMX) then
    AT(iAtGMX) = QM
  else if (iGrpGMX == MMGMX) then
    AT(iAtGMX) = MMO
  else
    Message = 'Fetch_QMMM: unknown index group'
    call WarningMessage(2,Message)
    call Abend()
  end if
end do
do iCastMM=1,nCastMM
  if (AT(CastMM(iCastMM)) == QM) then
    Message = 'Attempting to cast atom of QM type'
    call WarningMessage(2,Message)
    call Quit_OnUserError()
  else if (AT(CastMM(iCastMM)) == MMI) then
    Message = 'Attempting to cast atom of MMI type'
    call WarningMessage(2,Message)
    call Quit_OnUserError()
  end if
  AT(CastMM(iCastMM)) = MMI
end do
do iAtGMX=1,nAtGMX
  if ((AT(iAtGMX) == QM) .or. (AT(iAtGMX) == MMI)) then
    nAtIn = nAtIn+1
  else
    nAtOut = nAtOut+1
  end if
end do

! Put QM and inner MM atoms in xyz file and outer MM atoms on runfile
LuXYZ = 1
LuXYZ = isFreeUnit(LuXYZ)
call molcas_open(LuXYZ,'GMX.XYZ')
write(LuXYZ,'(i5,/)') nAtIn
call mma_allocate(CoordMMO,3,nAtOut)
call mma_allocate(LabMMO,nAtOut)
iAtOut = 1
do iAtGMX=1,nAtGMX
  iAtNmbGMX = mmslave_get_atomnumber(ipGMS,int(iAtGMX-1,kind=c_int))
  iFirst = index(PTab(iAtNmbGMX),' ')+1
  if (AT(iAtGMX) == QM) then
    Symbol = PTab(iAtNmbGMX)(iFirst:2)
  else
    Symbol = PTab(iAtNmbGMX)(iFirst:2)//'_MM'
  end if
  if ((AT(iAtGMX) == QM) .or. (AT(iAtGMX) == MMI)) then
    write(LuXYZ,'(a8,3ES24.15E3)') Symbol(1:5),CoordGMX(:,iAtGMX)*NmToAng
  else
    CoordMMO(:,iAtOut) = CoordGMX(:,iAtGMX)
    LabMMO(iAtOut) = Symbol
    iAtOut = iAtOut+1
  end if
end do
close(LuXYZ)
call Put_iArray('Atom Types',AT,nAtGMX)
CoordMMO(:,:) = CoordMMO/AuToNm
call Put_dArray('MMO Coords',CoordMMO,3*nAtOut)
call Put_cArray('MMO Labels',LabMMO(1),LenIn*nAtOut)

! Clean up
call mma_deallocate(CoordGMX)
call mma_deallocate(AT)
call mma_deallocate(CoordMMO)
call mma_deallocate(LabMMO)
call mmslave_done(ipGMS)

return

contains

function mmslave_copyx_wrapper(gms,natoms,x)

  integer(kind=iwp) :: mmslave_copyx_wrapper
  type(c_ptr) :: gms
  integer(kind=c_int) :: natoms
  real(kind=wp), target :: x(*)
  interface
    function mmslave_copyx(gms,natoms,x) bind(C,NAME='mmslave_copyx_')
      import :: c_int, c_ptr
      integer(kind=c_int) :: mmslave_copyx
      type(c_ptr), value :: gms, x
      integer(kind=c_int), value :: natoms
    end function mmslave_copyx
  end interface

  mmslave_copyx_wrapper = mmslave_copyx(gms,natoms,c_loc(x(1)))

end function mmslave_copyx_wrapper

end subroutine Fetch_QMMM

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(Fetch_QMMM)

#endif
