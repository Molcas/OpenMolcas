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

subroutine RunGromacs(nAtIn,Coord,Mltp,DoFirst,MltOrd,Forces,Grad,Energy)

use, intrinsic :: iso_c_binding, only: c_int, c_loc, c_ptr
use espf_global, only: MMI, MMIterMax, MxExtPotComp, QM, TPRDefName
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten, Angstrom, auTokJmol, auTokJmolnm
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtIn, MltOrd
real(kind=wp), intent(in) :: Coord(3,nAtIn), Mltp(*)
logical(kind=iwp), intent(in) :: DoFirst, Forces
real(kind=wp), intent(out) :: Grad(3,nAtIn), Energy
integer(kind=iwp) :: iAtGMX, iAtIn, iAtOut, ic, iLast, iOk, j, LuExtPot, LuWr, nAtGMX, nAtOut
real(kind=wp) :: EnergyGMX, Energy2GMX, q
logical(kind=iwp) :: Found, isNotLast
character(len=256) :: LogFileName, Message, TPRFileName
character(len=12) :: ExtPotFormat
type(c_ptr) :: ipCR, ipGMS
integer(kind=iwp), allocatable :: AT(:)
real(kind=wp), allocatable :: CoordGMX(:,:), CoordMMO(:,:), Field2GMX(:,:), FieldGMX(:,:), Force2GMX(:,:), ForceGMX(:,:), &
                              GradMMO(:,:), Pot2GMX(:), PotGMX(:)
real(kind=wp), parameter :: AuToNm = Angstrom/Ten
integer(kind=iwp), external :: isFreeUnit
interface
  subroutine mmslave_done(gms) bind(C,NAME='mmslave_done_')
    use, intrinsic :: iso_c_binding, only: c_ptr
    type(c_ptr), value :: gms
  end subroutine mmslave_done
  function mmslave_init(cr,log_) bind(C,NAME='mmslave_init_')
    use, intrinsic :: iso_c_binding, only: c_char, c_ptr
    type(c_ptr) :: mmslave_init
    type(c_ptr), value :: cr
    character(kind=c_char) :: log_(*)
  end function mmslave_init
  function mmslave_read_tpr(tpr,gms) bind(C,NAME='mmslave_read_tpr_')
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr
    integer(kind=c_int) :: mmslave_read_tpr
    character(kind=c_char) :: tpr(*)
    type(c_ptr), value :: gms
  end function mmslave_read_tpr
  function mmslave_set_q(gms,id,q) bind(C,NAME='mmslave_set_q_')
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr
    integer(kind=c_int) :: mmslave_set_q
    type(c_ptr), value :: gms
    integer(kind=c_int), value :: id
    real(kind=c_double), value :: q
  end function mmslave_set_q
  function init_commrec() bind(C,NAME='init_commrec_')
    use, intrinsic :: iso_c_binding, only: c_ptr
    type(c_ptr) :: init_commrec
  end function init_commrec
end interface

LuExtPot = 1
LuWr = u6
Energy = Zero

write(ExtPotFormat,'(a4,i2,a6)') '(I4,',MxExtPotComp,'F13.8)'

! Get MMO coordinates and atom types from runfile
call Qpg_dArray('MMO Coords',Found,nAtOut)
nAtOut = nAtOut/3
call mma_allocate(CoordMMO,3,nAtOut)
if (Found) then
  call Get_dArray('MMO Coords',CoordMMO,3*nAtOut)
else
  Message = 'RunGromacs: MMO coordinates not found on runfile'
  call WarningMessage(2,Message)
  call Abend()
end if

call Qpg_iArray('Atom Types',Found,nAtGMX)
call mma_allocate(AT,nAtGMX)
if (Found) then
  call Get_iArray('Atom Types',AT,nAtGMX)
else
  Message = 'RunGromacs: Atom types not found on runfile'
  call WarningMessage(2,Message)
  call Abend()
end if

! Since this is the first time we use information from runfile, do a
! consistency check
if (nAtGMX /= nAtIn+nAtOut) then
  Message = 'RunGromacs: nAtGMX and nAtIn+nAtOut not equal'
  call WarningMessage(2,Message)
  call Abend()
end if

! Initialize Gromacs mmslave
ipCR = init_commrec()
call prgmtranslate('GMX.LOG',LogFileName,iLast)
LogFileName(iLast+1:iLast+1) = char(0)
ipGMS = mmslave_init(ipCR,LogFileName)

! Let Gromacs read tpr file
TPRFileName = TPRDefName
iLast = len_trim(TPRFileName)
TPRFileName(iLast+1:iLast+1) = char(0)
iOk = mmslave_read_tpr(TPRFileName,ipGMS)
if (iOk /= 1) then
  Message = 'RunGromacs: mmslave_read_tpr is not ok'
  call WarningMessage(2,Message)
  call Abend()
end if

! Perform microiterations, though only if (1) it is requested, (2) there
! are MMO atoms to optimize, (3) not last energy, (4) multipoles are
! available, and (5) no gradient calculation
if ((MMIterMax > 0) .and. (nAtOut > 0)) then
  isNotLast = SuperName(1:11) /= 'last_energy'
  if ((.not. Dofirst) .and. (.not. Forces) .and. isNotLast) call Opt_MMO(nAtIn,Coord,nAtOut,CoordMMO,nAtGMX,AT,ipGMS)
end if

! Trick: Set QM charges in Gromacs to zero (or, currently, to a very
! small number). This is needed to exclude their contribution to the
! external potential.
do iAtGMX=1,nAtGMX
  if (AT(iAtGMX) == QM) then
    iOk = mmslave_set_q(ipGMS,int(iAtGMX-1,kind=c_int),1.0e-10_wp)
    if (iOk /= 1) then
      Message = 'RunGromacs: mmslave_set_q is not ok'
      call WarningMessage(2,Message)
      call Abend()
    end if
  end if
end do

! Get energy, forces and external potential from Gromacs
call mma_allocate(CoordGMX,3,nAtGMX)
call mma_allocate(FieldGMX,3,nAtGMX)
call mma_allocate(ForceGMX,3,nAtGMX)
call mma_allocate(PotGMX,nAtGMX)
iAtIn = 1
iAtOut = 1
do iAtGMX=1,nAtGMX
  if ((AT(iAtGMX) == QM) .or. (AT(iAtGMX) == MMI)) then
    CoordGMX(:,iAtGMX) = Coord(:,iAtIn)
    iAtIn = iAtIn+1
  else
    CoordGMX(:,iAtGMX) = CoordMMO(:,iAtOut)
    iAtOut = iAtOut+1
  end if
end do
CoordGMX(:,:) = CoordGMX/AuToNm
FieldGMX(:,:) = Zero
ForceGMX(:,:) = Zero
PotGMX(:) = Zero
iOk = mmslave_calc_energy_wrapper(ipGMS,CoordGMX,ForceGMX,FieldGMX,PotGMX,EnergyGMX)
if (iOk /= 1) then
  Message = 'RunGromacs: mmslave_calc_energy is not ok'
  call WarningMessage(2,Message)
  call Abend()
end if

! Special case: forces on MM atoms
if (Forces) then
  if (DoFirst) then
    Message = 'RunGromacs: something is wrong in the code'
    call WarningMessage(2,Message)
    call Abend()
  end if
  call mma_allocate(Field2GMX,3,nAtGMX)
  call mma_allocate(Force2GMX,3,nAtGMX)
  call mma_allocate(Pot2GMX,nAtGMX)
  Field2GMX(:,:) = Zero
  Force2GMX(:,:) = Zero
  Pot2GMX(:) = Zero
  ic = 1
  do iAtGMX=1,nAtGMX
    if (AT(iAtGMX) == QM) then
      q = Mltp(ic)
      iOk = mmslave_set_q(ipGMS,int(iAtGMX-1,kind=c_int),q)
      if (iOk /= 1) then
        Message = 'RunGromacs: mmslave_set_q is not ok'
        call WarningMessage(2,Message)
        call Abend()
      end if
      ic = ic+MltOrd
    end if
  end do
  iOk = mmslave_calc_energy_wrapper(ipGMS,CoordGMX,Force2GMX,Field2GMX,Pot2GMX,Energy2GMX)
  if (iOk /= 1) then
    Message = 'RunGromacs: mmslave_calc_energy is not ok'
    call WarningMessage(2,Message)
    call Abend()
  end if
end if

! Store classical contributions to energy and gradient
Energy = EnergyGMX/auTokJmol
if (Forces) then
  call mma_allocate(GradMMO,3,nAtOut)
  iAtIn = 1
  iAtOut = 1
  do iAtGMX=1,nAtGMX
    if (AT(iAtGMX) == QM) then
      Grad(:,iAtIn) = ForceGMX(:,iAtGMX)
      iAtIn = iAtIn+1
    else if (AT(iAtGMX) == MMI) then
      Grad(:,iAtIn) = Force2GMX(:,iAtGMX)
      iAtIn = iAtIn+1
    else
      GradMMO(:,iAtOut) = Force2GMX(:,iAtGMX)
      iAtOut = iAtOut+1
    end if
  end do
  Grad(:,:) = Grad(:,:)/auTokJmolnm
  GradMMO(:,:) = GradMMO(:,:)/auTokJmolnm
  call Put_dArray('MMO Grad',GradMMO,3*nAtOut)
end if

! Write external potential to ESPF.EXTPOT file
LuExtPot = isFreeUnit(LuExtPot)
call molcas_open(LuExtPot,'ESPF.EXTPOT')
write(LuExtPot,'(I1)') 0
iAtIn = 1
do iAtGMX=1,nAtGMX
  if ((AT(iAtGMX) == QM) .or. (AT(iAtGMX) == MMI)) then
    write(LuExtPot,ExtPotFormat) iAtIn,PotGMX(iAtGMX)/auTokJmol,-FieldGMX(:,iAtGMX)/auTokJmolnm,(Zero,j=5,MxExtPotComp)
    iAtIn = iAtIn+1
  end if
end do
close(LuExtPot)

! Clean up
call mma_deallocate(CoordMMO)
call mma_deallocate(AT)
call mma_deallocate(CoordGMX)
call mma_deallocate(FieldGMX)
call mma_deallocate(ForceGMX)
call mma_deallocate(PotGMX)
if (Forces) then
  call mma_deallocate(Field2GMX)
  call mma_deallocate(Force2GMX)
  call mma_deallocate(Pot2GMX)
  call mma_deallocate(GradMMO)
end if
call mmslave_done(ipGMS)

return

contains

function mmslave_calc_energy_wrapper(gms,x,f,A,phi,energy)
  integer(kind=iwp) :: mmslave_calc_energy_wrapper
  type(c_ptr) :: gms
  real(kind=wp), target :: x(*), f(*), A(*), phi(*)
  real(kind=wp) :: energy
  interface
    function mmslave_calc_energy(gms,x,f,A,phi,energy) bind(C,NAME='mmslave_calc_energy_')
      use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr
      integer(kind=c_int) :: mmslave_calc_energy
      type(c_ptr), value :: gms, x, f, A, phi
      real(kind=c_double) :: energy
    end function mmslave_calc_energy
  end interface
  mmslave_calc_energy_wrapper = mmslave_calc_energy(gms,c_loc(x(1)),c_loc(f(1)),c_loc(A(1)),c_loc(phi(1)),energy)
end function mmslave_calc_energy_wrapper

end subroutine RunGromacs

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(RunGromacs)

#endif
