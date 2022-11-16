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

subroutine Opt_MMO(nAtIn,Coord,nAtOut,CoordMMO,nAtGMX,AT,ipGMS)

use, intrinsic :: iso_c_binding, only: c_loc, c_ptr
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtIn, nAtOut, nAtGMX, AT(nAtGMX)
real(kind=wp), intent(in) :: Coord(3,nAtIn)
real(kind=wp), intent(inout) :: CoordMMO(3,nAtOut)
type(c_ptr), intent(in) :: ipGMS
#include "espf.fh"
#include "opt_mmo.fh"
integer(kind=iwp) :: i, iAtIn, iAtOut, iOk, iPL, j, MMIter
real(kind=wp) :: EnergyGMX, MaxF, OldEn, PotGMX(1), Step
character(len=256) :: Message
real(kind=wp), allocatable :: CoordGMX(:,:), FieldGMX(:,:), ForceGMX(:,:), GradMMO(:,:), NewCoord(:,:), OldCoord(:,:)
real(kind=wp), parameter :: TinyStep = 1.0e-50_wp*AuToNm
integer(kind=iwp), external :: iPL_espf
real(kind=wp), external :: ddot_

iPL = iPL_espf()

if (nAtOut == 0) then
  Message = 'Opt_MMO: subroutine called with nAtOut=0'
  call WarningMessage(2,Message)
  call Abend()
end if

! Set up arrays for use with Gromacs
call mma_allocate(CoordGMX,3,nAtGMX)
call mma_allocate(ForceGMX,3,nAtGMX)
call mma_allocate(FieldGMX,3,nAtGMX)
iAtIn = 1
iAtOut = 1
do i=1,nAtGMX
  if ((AT(i) == QM) .or. (AT(i) == MMI)) then
    call dcopy_(3,Coord(1,iAtIn),1,CoordGMX(1,i),1)
    iAtIn = iAtIn+1
  else
    call dcopy_(3,CoordMMO(1,iAtOut),1,CoordGMX(1,i),1)
    iAtOut = iAtOut+1
  end if
end do
call dscal_(3*nAtGMX,AuToNm,CoordGMX,1)

! Set up arrays with MMO coordinates and gradient
call mma_allocate(NewCoord,3,nAtOut)
call mma_allocate(OldCoord,3,nAtOut)
call mma_allocate(GradMMO,3,nAtOut)
call dcopy_(3*nAtOut,CoordMMO,1,NewCoord,1)
call dscal_(3*nAtOut,AuToNm,NewCoord,1)

if (iPL >= 2) then
  call CollapseOutput(1,'Gromacs microiterations')
  write(u6,*)
  write(u6,*) 'Initial coordinates (angstrom)'
  write(u6,*) '------------------------------'
  do i=1,nAtGMX
    if (AT(i) == QM) then
      write(u6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'QM'
    else if (AT(i) == MMI) then
      write(u6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMI'
    else
      write(u6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMO'
    end if
  end do
  write(u6,*)
end if

! MMO optimization cycle, performed in Gromacs units
Step = 0.1_wp*AuToNm
OldEn = huge(OldEn)
MaxF = huge(MaxF)
MMIter = 0
do while ((MMIter < MMIterMax) .and. (MaxF > ConvF) .and. (Step > TinyStep))
  MMIter = MMIter+1
  ! Get gradient from Gromacs
  iOk = mmslave_calc_energy_wrapper(ipGMS,CoordGMX,ForceGMX,FieldGMX,PotGMX,EnergyGMX)
  if (iOk /= 1) then
    Message = 'Opt_MMO: mmslave_calc_energy is not ok'
    call WarningMessage(2,Message)
    call Abend()
  end if
  iAtOut = 1
  do i=1,nAtGMX
    if (AT(i) == MMO) then
      GradMMO(1,iAtOut) = -ForceGMX(1,i)
      GradMMO(2,iAtOut) = -ForceGMX(2,i)
      GradMMO(3,iAtOut) = -ForceGMX(3,i)
      iAtOut = iAtOut+1
    end if
  end do
# ifdef _DEBUGPRINT_
  if (MMIter == 1) then
    write(u6,*) 'Iter       E          |F|         Fmax       Step'
    write(u6,*) '                     (Atomic units)'
    write(u6,*) '-----------------------------------------------------'
  end if
  if ((mod(MMIter,10) == 0) .or. (MMIter == 1)) then
    write(u6,200) MMIter,EnergyGMX/AuToKjPerMol,sqrt(DDot_(3*nAtOut,GradMMO,1,GradMMO,1))/AuToKjPerMolNm,MaxF/AuToKjPerMolNm, &
                  Step/AuToNm
  end if
# endif
  ! Steepest descent with adaptive step a la Gromacs
  if ((EnergyGMX < OldEn) .or. (MMIter == 1)) then
    Step = 1.2_wp*Step
    MaxF = Zero
    do iAtOut=1,nAtOut
      do j=1,3
        MaxF = max(MaxF,abs(GradMMO(j,iAtOut)))
      end do
    end do
    OldEn = EnergyGMX
    call dcopy_(3*nAtOut,NewCoord,1,OldCoord,1)
    call daxpy_(3*nAtOut,-Step/MaxF,GradMMO,1,NewCoord,1)
  else
    Step = 0.2_wp*Step
    call dcopy_(3*nAtOut,OldCoord,1,NewCoord,1)
    call daxpy_(3*nAtOut,-Step/MaxF,GradMMO,1,NewCoord,1)
  end if
  ! Update coordinates for Gromacs
  iAtOut = 1
  do i=1,nAtGMX
    if (AT(i) == MMO) then
      CoordGMX(1,i) = NewCoord(1,iAtOut)
      CoordGMX(2,i) = NewCoord(2,iAtOut)
      CoordGMX(3,i) = NewCoord(3,iAtOut)
      iAtOut = iAtOut+1
    end if
  end do
end do

! Issue warning if calculation stopped for other reason than small F
if (MaxF > ConvF) then
  if (MMIter == MMIterMax) then
    Message = 'Maximum number of microiterations reached'
    call WarningMessage(1,Message)
  else if (Step <= TinyStep) then
    Message = 'Microiterations stopped due to zero step size'
    call WarningMessage(1,Message)
  end if
end if

! Undo the last move to recover the "converged" coordinates
if (MMIter > 0) then
  call dcopy_(3*nAtOut,OldCoord,1,NewCoord,1)
  iAtOut = 1
  do i=1,nAtGMX
    if (AT(i) == MMO) then
      CoordGMX(1,i) = NewCoord(1,iAtOut)
      CoordGMX(2,i) = NewCoord(2,iAtOut)
      CoordGMX(3,i) = NewCoord(3,iAtOut)
      iAtOut = iAtOut+1
    end if
  end do
end if

if (iPL >= 2) then
  write(u6,*)
  write(u6,300) MMIter
  write(u6,400) MaxF/AuToKjPerMolNm,ConvF/AuToKjPerMolNm
  call dscal_(3*nAtGMX,-One/AuToKjPerMolNm,ForceGMX,1)
  EnergyGMX = EnergyGMX/AuToKjPerMol
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Properties'
  write(u6,*) '----------'
  write(u6,*) 'Energy: ',EnergyGMX
  do i=1,nAtGMX
    if (AT(i) == MMO) then
      write(u6,500) (ForceGMX(j,i),j=1,3)
    else
      write(u6,500) Zero,Zero,Zero
    end if
  end do
# endif
  write(u6,*)
end if

! Put optimized MMO coordinates on runfile
call dcopy_(3*nAtOut,NewCoord,1,CoordMMO,1)
call dscal_(3*nAtOut,One/AuToNm,CoordMMO,1)
call Put_dArray('MMO Coords',CoordMMO,3*nAtOut)

if (iPL >= 2) then
  if (MMIter > 1) then
    write(u6,*) 'Final coordinates (angstrom)'
    write(u6,*) '----------------------------'
    do i=1,nAtGMX
      if (AT(i) == QM) then
        write(u6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'QM'
      else if (AT(i) == MMI) then
        write(u6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMI'
      else
        write(u6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMO'
      end if
    end do
  end if
  call CollapseOutput(0,'Gromacs microiterations')
  write(u6,*)
end if

! Clean up
call mma_deallocate(CoordGMX)
call mma_deallocate(ForceGMX)
call mma_deallocate(FieldGMX)
call mma_deallocate(NewCoord)
call mma_deallocate(OldCoord)
call mma_deallocate(GradMMO)

return

100 format(I6,3(F12.6,1X),2X,A)
300 format(' Performed: ',I6,' MM iterations')
400 format(' Max. force: ',F12.6,' (convergence at: ',F12.6,')')
#ifdef _DEBUGPRINT_
200 format(I5,6ES12.4)
500 format(3(F12.6,1X))
#endif

contains

function mmslave_calc_energy_wrapper(gms,x,f,A,phi,energy)
  integer :: mmslave_calc_energy_wrapper
  type(c_ptr) :: gms
  real*8, target :: x(*), f(*), A(*), phi(*)
  real*8 :: energy
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

end subroutine Opt_MMO

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(Opt_MMO)

#endif
