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
! Copyright (C) 2016, Sebastian Wouters                                *
!               2016, Quan Phung                                       *
!***********************************************************************
! CheMPS2-Molcas main interface
! Based on Block interface, written by N. Nakatani
! Written by Quan Phung and Sebastian Wouters, Leuven, Aug 2016
! Adapted for Molcas 8.1 by Quan Phung, Leuven, Oct 2016

subroutine Chemps2Ctl(W1,TUVX,IFINAL,IRST)

#ifdef _MOLCAS_MPP_
use MPI, only: MPI_COMM_WORLD
use Para_Info, only: Is_Real_Par, King
use Definitions, only: MPIInt
#endif
use rasscf_data, only: CBLBM, chemps2_blb, chemps2_lrestart, chemps2_noise, chemps2_restart, davidson_tol, Do3RDM, ENER, iCIonly, &
                       iOrbTyp, ITER, lroots, max_canonical, max_sweep, mxSym, MxDMRG, NAC, THRE, hfocc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Five, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: W1(*), TUVX(*)
integer(kind=iwp), intent(in) :: IFINAL, IRST
#include "general.fh"
integer(kind=iwp) :: iChMolpro(8), LINSIZE, NUM_TEI, dtemp, nooctemp, labelpsi4, conversion(8), activesize(8), chemroot, &
                     chemps2_info, iOper(0:7), ihfocc, iErr, iOrb, iSigma, iSym, jOrb, lSymMolpro, LUCHEMIN, LUCONV, LUTOTE, &
                     nIrrep, NRDM_ORDER
integer(kind=iwp), allocatable :: OrbSym(:)
#ifdef _MOLCAS_MPP_
integer(kind=MPIInt) :: IERROR
#endif
real(kind=wp) :: chemps2_totale_4d, revdiff, chemps2_conv
logical(kind=iwp) :: fiedler, Found, mps0
character(len=3) :: Label
character(len=10) :: rootindex
character(len=100) :: imp1, imp2
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

NRDM_ORDER = 2
if (NACTEL == 1) NRDM_ORDER = 1

!*********************
!  WRITEOUT FCIDUMP  *
!*********************

LINSIZE = (NAC*(NAC+1))/2
NUM_TEI = (LINSIZE*(LINSIZE+1))/2
call FCIDUMP_OUTPUT(NAC,NACTEL,ISPIN-1,lSymMolpro,OrbSym,Zero,W1,TUVX,LINSIZE,NUM_TEI)

call mma_deallocate(OrbSym)

!*************************
!  WRITEOUT ACTIVE FOCK  *
!*************************

!write(u6,*) 'Currently the Fock matrix is printed in fckpt2.f'

!************************
!  WRITEOUT INPUT FILE  *
!************************

#ifdef _MOLCAS_MPP_
if (KING() .or. (.not. Is_Real_Par())) then
#endif
  if (IRST == 0) then
    ! Cleanup chemps2.log.total
    imp1 = 'chemps2.log.total'
    call f_inquire(imp1,Found)
    if (Found) call aixrm(imp1)
    ! Check if checkpoint files exist
    if (chemps2_restart) then
      call f_inquire('CHEMNATFIE',fiedler)
      call f_inquire('CHEMNATMPS0',mps0)
      if (fiedler .and. mps0) then
        write(u6,*) 'CHEMPS2> Found checkpoint files for DMRG-SCF'
        ! Copy CheMPS2_natorb_MPSxxx.h5 to CheMPS2_MPSxxx.h5
        call fcopy('CHEMNATFIE','CHEMFIE',iErr)
        do chemroot=1,lroots
          write(rootindex,'(i2)') chemroot-1
          imp1 = 'CheMPS2_natorb_MPS'//trim(adjustl(rootindex))//'.h5'
          imp2 = 'CheMPS2_MPS'//trim(adjustl(rootindex))//'.h5'
          call fcopy(imp1,imp2,iErr)
        end do
      else
        ! Reset chemps2_restart = .false. if not checkpoint files
        write(u6,*) 'CHEMPS2> No checkpoint files for DMRG-SCF'
        chemps2_restart = .false.
      end if
    end if
    ! Check if checkpoint files for 3RDM exist
    if (chemps2_lrestart == 1) then
      call f_inquire('CHEMCANFIE',fiedler)
      call f_inquire('CHEMCANMPS0',mps0)
      if (fiedler .and. mps0) then
        write(u6,*) 'CHEMPS2> Found checkpoint files for n-RDM'
      else
        write(u6,*) 'CHEMPS2> No checkpoint files for n-RDM'
        chemps2_lrestart = 0
      end if
    end if
  end if
#ifdef _MOLCAS_MPP_
end if
#endif

LUCHEMIN = isFreeUnit(29)
call molcas_open(LUCHEMIN,'chemps2.input')
write(LUCHEMIN,*) 'FCIDUMP = FCIDUMP_CHEMPS2'

call group_psi4number(Label,Labelpsi4)
write(LUCHEMIN,'(1x,a8,i1)') 'GROUP = ',Labelpsi4
write(LUCHEMIN,*)

write(LUCHEMIN,'(1x,a13,i2)') 'EXCITATION = ',lRoots-1
write(LUCHEMIN,*)

if (((abs(CBLBM) > chemps2_blb) .and. (IFINAL /= 2)) .or. &
    ((IRST == 0) .and. (.not. chemps2_restart)) .or. &
    ((IFINAL == 2) .and. Do3RDM .and. (chemps2_lrestart == 0)) .or. &
    ((IFINAL == 2) .and. (iOrbTyp == 2) .and. (chemps2_lrestart == 0))) then

  imp1 = 'molcas_fiedler.txt'
  call f_inquire(imp1,Found)
  if (Found) call aixrm(imp1)
  !Quan: FIXME: how to remove CheMPS2_MPS0.h5, etc with aixrm
  call systemf('rm -f CheMPS2_MPS*.h5',iErr)
  write(LUCHEMIN,*) 'MOLCAS_MPS     = TRUE'
  write(u6,*) 'CHEMPS2> Start DMRG from scratch'

  write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_STATES       = '
  dtemp = 500
  do
    write(LUCHEMIN,('(i7,a2)'),advance='NO') dtemp,','
    dtemp = dtemp+min(dtemp,1000)
    if (dtemp >= MxDMRG) then
      write(LUCHEMIN,'(i7,a2,i7)') MxDMRG,',',MxDMRG
      exit
    end if
  end do

  write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_ENERGY_CONV  = '
  dtemp = 500
  do
    if (dtemp == 500) then
      write(LUCHEMIN,'(es12.5,a3)',advance='NO') THRE*1.0e3_wp,','
    else
      write(LUCHEMIN,'(es12.5,a3)',advance='NO') THRE*1.0e2_wp,','
    end if

    dtemp = dtemp+min(dtemp,1000)
    if (dtemp >= MxDMRG) then
      write(LUCHEMIN,'(es12.5,a3,es12.5)') THRE*Five,',',THRE*Half
      exit
    end if
  end do

  write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_MAX_SWEEPS   = '
  dtemp = 500
  do
    if (dtemp == 500) then
      write(LUCHEMIN,'(i7,a3)',advance='NO') max_sweep,','
    else
      write(LUCHEMIN,'(i7,a3)',advance='NO') max_sweep/2,','
    end if

    dtemp = dtemp+min(dtemp,1000)
    if (dtemp >= MxDMRG) then
      if (IFINAL == 2) then
        if (Do3RDM .or. (iOrbTyp == 2)) then
          write(LUCHEMIN,'(i7,a3,i7)') max_sweep/2,',',max_canonical
        else
          write(LUCHEMIN,'(i7,a3,i7)') max_sweep/2,',',max_sweep*5
        end if
      else
        write(LUCHEMIN,'(i7,a3,i7)') max_sweep/2,',',max_sweep
      end if
      exit
    end if
  end do

  write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_NOISE_PREFAC = '
  dtemp = 500
  do
    write(LUCHEMIN,'(es12.5,a3)',advance='NO') chemps2_noise,','
    dtemp = dtemp+min(dtemp,1000)
    if (dtemp >= MxDMRG) then
      write(LUCHEMIN,'(es12.5,a10)') chemps2_noise,' ,   0.00'
      exit
    end if
  end do

  write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_DVDSON_RTOL  = '
  dtemp = 500
  do
    if (dtemp == 500) then
      write(LUCHEMIN,'(a9)',advance='NO') '1.0e-3 ,'
    else
      write(LUCHEMIN,'(a9)',advance='NO') '1.0e-4 ,'
    end if

    dtemp = dtemp+min(dtemp,1000)
    if (dtemp >= MxDMRG) then
      write(LUCHEMIN,'(a9,es12.5)') '1.0e-4 ,',davidson_tol
      exit
    end if
  end do
  write(LUCHEMIN,*)

else
  ! DMRG restart with fixed orbital order
  if (((abs(CBLBM) > chemps2_blb/Ten) .and. (IFINAL /= 2)) .or. ((IRST == 0) .and. chemps2_restart)) then
    write(u6,*) 'CHEMPS2> Partial restart DMRG from previous step'

    write(LUCHEMIN,*) 'MOLCAS_MPS     = TRUE'
    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_STATES       = '
    write(LUCHEMIN,'(i7,a2,i7)') MxDMRG,',',MxDMRG

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_ENERGY_CONV  = '
    write(LUCHEMIN,'(es12.5,a3,es12.5)') THRE*Five,',',THRE*Half

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_MAX_SWEEPS   = '
    if (IFINAL == 2) then
      if (Do3RDM .or. (iOrbTyp == 2)) then
        write(LUCHEMIN,'(i7,a3,i7)') max_sweep/2,',',max_canonical
      else
        write(LUCHEMIN,'(i7,a3,i7)') max_sweep/2,',',max_sweep*5
      end if
    else
      write(LUCHEMIN,'(i7,a3,i7)') max_sweep/2,',',max_sweep
    end if

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_NOISE_PREFAC = '
    write(LUCHEMIN,'(es12.5,a10)') chemps2_noise,' ,   0.00'

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_DVDSON_RTOL  = '
    write(LUCHEMIN,'(a9,es12.5)') '1.0e-4 ,',davidson_tol
    write(LUCHEMIN,*)
  else
    !write(u6,*) 'Full Restart'
    write(u6,*) 'CHEMPS2> Fully restart DMRG from previous step'
    write(LUCHEMIN,*) 'MOLCAS_MPS     = TRUE'
    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_STATES       = '
    write(LUCHEMIN,'(i7)') MxDMRG

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_ENERGY_CONV  = '
    write(LUCHEMIN,'(es12.5)') THRE*Half

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_MAX_SWEEPS   = '
    if ((IFINAL == 2) .or. ((IFINAL == 1) .and. (iCIonly == 1))) then
      if ((IFINAL == 2) .and. (Do3RDM .or. (iOrbTyp == 2))) then
        write(LUCHEMIN,'(i7)') max_canonical
      else
        write(LUCHEMIN,'(i7)') max_sweep*5
      end if
    else
      write(LUCHEMIN,'(i7)') max_sweep
    end if

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_NOISE_PREFAC = '
    write(LUCHEMIN,'(a16)') '   0.00'

    write(LUCHEMIN,'(1x,a21)',advance='NO') 'SWEEP_DVDSON_RTOL  = '
    write(LUCHEMIN,'(es12.5)') davidson_tol
    write(LUCHEMIN,*)
  end if
end if

write(LUCHEMIN,'(a8)',advance='nO') 'NOCC = '
do nooctemp=1,NSYM-1
  write(LUCHEMIN,'(a5)',advance='NO') '0 ,'
end do
write(LUCHEMIN,'(a3)') '0'

write(LUCHEMIN,'(a8)',advance='NO') 'NACT = '
call molpro2psi(Label,conversion)
do iSym=1,nSym
  activesize(conversion(iChMolpro(iSym))) = NASH(iSym)
end do
do nooctemp=1,NSYM-1
  write(LUCHEMIN,'(i3,a2)',advance='NO') activesize(nooctemp),' ,'
end do
write(LUCHEMIN,'(i3)') activesize(NSYM)

write(LUCHEMIN,'(a8)',advance='NO') 'NVIR = '
do nooctemp=1,NSYM-1
  write(LUCHEMIN,'(a5)',advance='NO') '0 ,'
end do
write(LUCHEMIN,'(a3)') '0'
write(LUCHEMIN,*)

write(LUCHEMIN,*) 'MOLCAS_FIEDLER = TRUE'
write(LUCHEMIN,*) 'MOLCAS_STATE_AVG = TRUE'
write(LUCHEMIN,*) 'MOLCAS_2RDM    = molcas_2rdm.h5'

if (sum(hfocc) == NACTEL) then
  write(u6,*) 'CHEMPS2> Using user-specified ROHF guess'
  write(LUCHEMIN,'(a13)',advance='NO') 'MOLCAS_OCC ='
  do ihfocc=1,NAC-1
    write(LUCHEMIN,'(i3,a2)',advance='NO') HFOCC(ihfocc),', '
  end do
  write(LUCHEMIN,'(i3)') HFOCC(NAC)
else
  write(u6,*) 'CHEMPS2> Using noise guess'
end if

if ((IFINAL == 2) .and. Do3RDM .and. (NACTEL > 2)) then
  write(u6,*) 'CHEMPS2> Running 3-RDM and F.4-RDM'
  write(LUCHEMIN,*) 'MOLCAS_3RDM    = molcas_3rdm.h5'
  write(LUCHEMIN,*) 'MOLCAS_F4RDM   = molcas_f4rdm.h5'
  write(LUCHEMIN,*) 'MOLCAS_FOCK    = FOCK_CHEMPS2'
end if

write(LUCHEMIN,*)

write(LUCHEMIN,*) 'PRINT_CORR = TRUE'
write(LUCHEMIN,*) 'TMP_FOLDER = ./'

close(LUCHEMIN)

#ifdef _MOLCAS_MPP_
write(u6,'(1x,a21,i3)') 'CHEMPS2> ITERATION : ',ITER
if (KING() .or. (.not. Is_Real_Par())) then
#endif

  ! Quan: overwrite CheMPS2_xxxorb_MPSX.h5 to CheMPS2_MPSX.h5
  if (((IFINAL == 2) .and. Do3RDM .and. (chemps2_lrestart > 0)) .or. &
      ((IFINAL == 2) .and. (iOrbTyp == 2) .and. (chemps2_lrestart > 0))) then
    if (chemps2_lrestart == 1) then
      write(u6,*) 'CHEMPS2> Using user-supplied checkpoint files'
      call fcopy('CHEMCANFIE','CHEMFIE',iErr)
      do chemroot=1,lroots
        write(rootindex,'(i2)') chemroot-1
        imp1 = 'CheMPS2_canorb_MPS'//trim(adjustl(rootindex))//'.h5'
        imp2 = 'CheMPS2_MPS'//trim(adjustl(rootindex))//'.h5'
        call fcopy(imp1,imp2,iErr)
      end do
    end if

    if (chemps2_lrestart == 2) then
      write(u6,*) 'CHEMPS2> Using checkpoint files from previous step (not recommended)'
      call fcopy('CHEMNATFIE','CHEMFIE',iErr)
      do chemroot=1,lroots
        write(rootindex,'(i2)') chemroot-1
        imp1 = 'CheMPS2_natorb_MPS'//trim(adjustl(rootindex))//'.h5'
        imp2 = 'CheMPS2_MPS'//trim(adjustl(rootindex))//'.h5'
        call fcopy(imp1,imp2,iErr)
      end do
    end if
  end if

  ! Quan: save CANORB before actually calculating
  if (((IFINAL == 2) .and. Do3RDM) .or. ((IFINAL == 2) .and. (iOrbTyp == 2))) then
    write(u6,*) 'CHEMPS2> Save CANORB'
    ! Quan: FIXME: Bug!
    !call OrbFiles(JOBIPH,IPRLEV)
  end if

  call systemf('chemps2 --file=chemps2.input > chemps2.log',iErr)
  call systemf('cat chemps2.log >> chemps2.log.total',iErr)

  ! Quan: save natorb checkpoint file in all iteration
  if (IFINAL < 2) then
    if (IFINAL == 1) then
      write(u6,*) 'CHEMPS2> Save natorb checkpoint files'
    end if
    call fcopy('CHEMFIE','CHEMNATFIE',iErr)
    do chemroot=1,lroots
      write(rootindex,'(i2)') chemroot-1
      imp1 = 'CheMPS2_natorb_MPS'//trim(adjustl(rootindex))//'.h5'
      imp2 = 'CheMPS2_MPS'//trim(adjustl(rootindex))//'.h5'
      call fcopy(imp2,imp1,iErr)
    end do
  end if

  ! Quan: save canorb checkpoint file if possible
  if (((IFINAL == 2) .and. Do3RDM) .or. ((IFINAL == 2) .and. (iOrbTyp == 2))) then

    write(u6,*) 'CHEMPS2> Save canorb checkpoint files'
    call fcopy('CHEMFIE','CHEMCANFIE',iErr)
    do chemroot=1,lroots
      write(rootindex,'(i2)') chemroot-1
      imp1 = 'CheMPS2_canorb_MPS'//trim(adjustl(rootindex))//'.h5'
      imp2 = 'CheMPS2_MPS'//trim(adjustl(rootindex))//'.h5'
      call fcopy(imp2,imp1,iErr)
    end do
  end if

  ! Quan: Cleanup checkpoint files
  if (IFINAL == 2) then
    imp1 = 'molcas_fiedler.txt'
    call f_inquire(imp1,Found)
    if (Found) call aixrm(imp1)
    !Quan: FIXME: how to remove CheMPS2_MPS0.h5, etc with aixrm
    call systemf('rm -f CheMPS2_MPS*.h5',iErr)
  end if

#ifdef _MOLCAS_MPP_
end if

if (Is_Real_Par()) then
  call MPI_Barrier(MPI_COMM_WORLD,IERROR)
end if

!Quan: FIXME: softlink all the n-RDM files
if (Is_Real_Par() .and. (.not. KING())) then
  do chemroot=1,lroots
    write(rootindex,'(i2)') chemroot-1
    imp1 = 'ln -sf ../molcas_2rdm.h5.r'//trim(adjustl(rootindex))//' .'
    call systemf(imp1,iErr)
    imp1 = 'ln -sf ../molcas_3rdm.h5.r'//trim(adjustl(rootindex))//' .'
    call systemf(imp1,iErr)
    imp1 = 'ln -sf ../molcas_f4rdm.h5.r'//trim(adjustl(rootindex))//' .'
    call systemf(imp1,iErr)
    imp1 = 'ln -sf ../CheMPS2_natorb_MPS0.h5 .'
    call systemf(imp1,iErr)
    imp1 = 'ln -sf ../CheMPS2_canorb_MPS0.h5 .'
    call systemf(imp1,iErr)
    imp1 = 'ln -sf ../molcas_natorb_fiedler.txt .'
    call systemf(imp1,iErr)
    imp1 = 'ln -sf ../molcas_canorb_fiedler.txt .'
    call systemf(imp1,iErr)
  end do
  call systemf('ln -sf ../chemps2.log .',iErr)
end if
#endif

!Quan: a very dirty way to extract the total energy
call systemf('grep "***  2-RDM" -B 5 chemps2.log | grep "all instructions" | cut -c 61- > chemps2_totale_4d',iErr)

!Quan: fix bug E(FCI) != E(CASSCF)
call systemf('grep "Econst" chemps2.log | cut -c 39- > chemps2_totale',iErr)

LUTOTE = isFreeUnit(30)
call molcas_open(LUTOTE,'chemps2_totale')

!Quan: write energy to ENER
do chemroot=1,lroots
  read(LUTOTE,*) ENER(chemroot,ITER)
end do
close(LUTOTE)

!Quan: check the difference between ener and ener_4d
LUTOTE = isFreeUnit(30)
call molcas_open(LUTOTE,'chemps2_totale_4d')
do chemroot=1,lroots
  read(LUTOTE,*) chemps2_totale_4d
  revdiff = abs(chemps2_totale_4d-ENER(chemroot,ITER))/chemps2_totale_4d
  if (revdiff > 1.0e-9_wp) then
    write(u6,*) 'CHEMPS2> large (E(4m) - E(m))/E(4m) = ',revdiff,'for root',chemroot,', consider increasing m!'
  end if
end do
close(LUTOTE)

!Quan: check chemps2 convergence
write(rootindex,'(i2)') lroots+9
imp1 = ''
imp1 = 'grep "***  2-RDM" -B '//trim(adjustl(rootindex))//' chemps2.log | grep "Energy difference" | cut -c 69- > chemps2_conv'
call systemf(imp1,iErr)

LUCONV = isFreeUnit(30)
call molcas_open(LUCONV,'chemps2_conv')
do chemroot=1,lroots
  read(LUCONV,*) chemps2_conv
  write(u6,'(1x,a14,i3,a30,es10.2)') 'CHEMPS2> Root ',chemroot,' :: DMRG energy convergence : ',chemps2_conv
  if (abs(chemps2_conv) > THRE*Half) then
    write(u6,*) 'CHEMPS2> DMRG not converged, consider increasing MXSWeep'
  end if
end do
close(LUCONV)

!Quan: check if CheMPS2 finished without error
imp1 = 'grep "Info on DMRG" chemps2.log | cut -c 43- > chemps2_info'
call systemf(imp1,iErr)

LUCONV = isFreeUnit(30)
call molcas_open(LUCONV,'chemps2_info')
read(LUCONV,*) chemps2_info
close(LUCONV)
if (chemps2_info /= 0) then
  write(u6,*) 'CHEMPS2> CheMPS2 ends abnormally, check calculation'
end if

return

end subroutine Chemps2Ctl
