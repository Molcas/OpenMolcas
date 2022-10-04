!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Quan Phung                                       *
!***********************************************************************
! Main control file for DICE. Template from CheMPS2 interface.

Subroutine DiceCtl(W1,TUVX,IFINAL,IRST)

#ifdef _MOLCAS_MPP_
use MPI, only: MPI_COMM_WORLD
use Para_Info, only: Is_Real_Par, King
use Definitions, only: MPIInt
#endif
use rasscf_data, only: ENER, ITER, lroots, NAC, nref_dice, dice_eps1, dice_eps2, &
                       dice_sampleN, dice_iter, dice_stoc, dice_restart, diceocc, mxSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: W1(*), TUVX(*)
integer(kind=iwp), intent(in) :: IFINAL, IRST
#include "general.fh"
integer(kind=iwp) :: iChMolpro(8), LINSIZE, NUM_TEI, chemroot, nIrrep, iref_dice, &
                     iOper(0:7), iErr, iOrb, iSigma, iSym, jOrb, lSymMolpro, LUDICEIN, LUTOTE
integer(kind=iwp), allocatable :: OrbSym(:)
#ifdef _MOLCAS_MPP_
integer(kind=MPIInt) :: IERROR
#endif
real(kind=wp) :: pt2ener
logical(kind=iwp) :: Found
character(len=3) :: Label, dice_nprocs
character(len=10) :: rootindex
character(len=150) :: imp1, imp2
integer(kind=iwp), external :: isFreeUnit

! Quan: FIXME: Do we need this?
! Load symmetry info from RunFile
iOper = 0
call Get_iScalar('NSYM',nIrrep)
call Get_iArray('Symmetry operations',iOper,nIrrep)
call Get_iScalar('Rotational Symmetry Number',iSigma)

! Get character table to convert MOLPRO symmetry format
call MOLPRO_ChTab(nSym,Label,iChMolpro)

! Convert orbital symmetry into MOLPRO format
call mma_allocate(OrbSym,NAC,label='OrbSym')
iOrb = 1
do iSym=1,nSym
  do jOrb=1,NASH(iSym)
    OrbSym(iOrb) = iChMolpro(iSym)
    iOrb = iOrb+1
  end do
end do
lSymMolpro = iChMolpro(stSym)

!*********************
!  WRITEOUT FCIDUMP  *
!*********************

LINSIZE = (NAC*(NAC+1))/2
NUM_TEI = (LINSIZE*(LINSIZE+1))/2
call FCIDUMP_OUTPUT(NAC,NACTEL,ISPIN-1,lSymMolpro,OrbSym,Zero,W1,TUVX,LINSIZE,NUM_TEI)

call mma_deallocate(OrbSym)
! Dice only reads FCIDUMP file
call systemf('ln -sf FCIDUMP_CHEMPS2 FCIDUMP',iErr)

!************************
!  WRITEOUT INPUT FILE  *
!************************
#ifdef _MOLCAS_MPP_
if (KING() .or. (.not. Is_Real_Par())) then
#endif
  if (IRST == 0) then
      ! Cleanup dice.out.total
      imp1 = 'dice.out.total'
      call f_inquire(imp1,Found)
      if (Found) call aixrm(imp1)
  endif
#ifdef _MOLCAS_MPP_
endif
#endif

write(u6,*) 'DICE> INTERATION : ', ITER
LUDICEIN = isFreeUnit(30)
call molcas_open(LUDICEIN,'input.dat')

write(LUDICEIN,'(a4,i4)') 'nocc', NACTEL
  do iref_dice=1,nref_dice
  write(LUDICEIN,'(a)') trim(diceocc(iref_dice))
enddo
write(LUDICEIN,'(a3)') 'end'
write(LUDICEIN,'(a6,i3)') 'nroots', lroots
write(LUDICEIN,*)
write(LUDICEIN,'(a8)') 'schedule'
write(LUDICEIN,'(a1,e12.5)') '0', dice_eps1*1.0d1
write(LUDICEIN,'(a1,e12.5)') '3', dice_eps1*1.0d1
write(LUDICEIN,'(a1,e12.5)') '6', dice_eps1
write(LUDICEIN,'(a3)') 'end'
write(LUDICEIN,'(a7,i6)') 'maxiter', dice_iter
write(LUDICEIN,'(a5)') 'DoRDM'
write(LUDICEIN,'(a8)') 'dE 1.e-8'
write(LUDICEIN,*)
write(LUDICEIN,'(a7,i6)') 'SampleN', dice_sampleN
write(LUDICEIN,'(a8,e12.5)') 'epsilon2', dice_eps2
write(LUDICEIN,'(a18)') 'targetError 8.0e-5'
if (IRST>0 .or. dice_restart) then
  write(LUDICEIN,'(a11)') 'fullrestart'
endif

if (.not. dice_stoc) then
  write(LUDICEIN,'(a13)') 'deterministic'
else
  write(LUDICEIN,'(a13,e12.5)') 'epsilon2Large', dice_eps2*10.0d0
endif

close(LUDICEIN)

!*****************************
!  RUN DICE                  *
!*****************************

#ifdef _MOLCAS_MPP_
if (KING() .or. (.not. Is_Real_Par())) then
#endif
  call get_environment_variable("MOLCAS_DICE", dice_nprocs, status=ierr)
  if (ierr == 0) then
     imp2 = 'mpirun -np '//trim(adjustl(dice_nprocs))//' Dice >dice.out 2>dice.err'
  else
     imp2 = 'Dice >dice.out 2>dice.err'
  endif

  call systemf(imp2,iErr)
  if (iErr /= 0) then
     write(u6,*) 'DICE> DICE ends abnormally, check calculation'
  endif
  call systemf('cat dice.out >> dice.out.total',iErr)
#ifdef _MOLCAS_MPP_
end if

if (Is_Real_Par()) then
  call MPI_Barrier(MPI_COMM_WORLD,IERROR)
end if

if (Is_Real_Par() .and. (.not. KING())) then
  do chemroot=1,lroots
    write(rootindex,'(i2)') chemroot-1
    imp1 = 'ln -sf ../spatialRDM.'//trim(adjustl(rootindex))//'.'//trim(adjustl(rootindex))//'.txt .'
    call systemf(imp1,iErr)
  end do
  call systemf("ln -sf ../dice.out .", iErr)
end if
#endif

!*****************************
!  EXTRACT ENERGY            *
!*****************************

if (.not. dice_stoc) then
  imp1 = 'grep PTEnergy dice.out | cut -c 10- > dice.energy'
  call systemf(imp1,iErr)
else
  do chemroot=1,lroots
    write(rootindex,'(i2)') chemroot+2
    imp1 = 'grep -A '//trim(adjustl(rootindex))// &
           ' "VARIATIONAL CALCULATION RESULT" dice.out | tail -n 1 | cut -c 6-30 >> dice.energy'
    call systemf(imp1,iErr)
  enddo
  imp2 = 'grep +/- dice.out | cut -c 10- > PT2.energy'
  call systemf(imp2,iErr)
endif

LUTOTE = isFreeUnit(30)
call molcas_open(LUTOTE,'dice.energy')

do chemroot=1,lroots
  read(LUTOTE,*) ENER(chemroot,ITER)
  write(u6,*) 'DICE> Deterministic PT2 Energy: ', ENER(chemroot,ITER)
enddo
close(LUTOTE)

if (dice_stoc) then
  call molcas_open(LUTOTE,'PT2.energy')

  do chemroot=1,lroots
    read(LUTOTE,*) pt2ener
    write(u6,*) 'DICE> Stochastic PT2 Energy: ', pt2ener
  enddo
  close(LUTOTE)
endif

return

end subroutine DiceCtl
