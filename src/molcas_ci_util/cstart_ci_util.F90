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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine CStart_CI_Util(C,h0,TUVX,iSel,ExplE,ExplV,nMaxSel,iFinal)
!***********************************************************************
!                                                                      *
!     Find initial CI-vectors                                          *
!                                                                      *
!     calling arguments:                                               *
!     C       : array of real*8                                        *
!               CI vector                                              *
!     h0      : array of real*8                                        *
!               one-electron integrals                                 *
!     TUVX    : array of real*8                                        *
!               two-electron integrals                                 *
!     iSel    : array of integer                                       *
!               list of configuration included in the explicit Hamilt. *
!     ExplE   : array of real*8                                        *
!               eigenvalues of the explicit Hamiltonian                *
!     ExplV   : array of real*8                                        *
!               eigenvectors of the explicit Hamiltonian               *
!     iFinal  : integer                                                *
!               status of optimization process                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_fetch_dset,mh5_close_file
#endif
use csfbas, only: CONF
use glbbas, only: CFTP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: C(*), ExplE(*), ExplV(*)
real(kind=wp), intent(in) :: h0(*), TUVX(*)
integer(kind=iwp), intent(_OUT_) :: iSel(*)
integer(kind=iwp), intent(in) :: nMaxSel, iFinal
integer(kind=iwp) :: i, iDisk, iJOB, IPRLEV, j, k, l
#ifdef _HDF5_
integer(kind=iwp) :: mh5id
#endif
logical(kind=iwp) :: Exists
character(len=80) :: String
integer(kind=iwp), allocatable :: vkcnf(:)
real(kind=wp), allocatable :: Tmp1(:)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "output_ras.fh"

IPRLEV = IPRLOC(3)

! special case: nConf=1

if ((nConf == 1) .and. (nAc == 0)) then
  iDisk = IADR15(4)
  !ExplE(1) = C(1)  ! Commented out by Jesper
  !ExplV(1) = One   ! Commented out by Jesper
  C(1) = One
  call Save_tmp_CI_vec(1,nConf,C,LuDavid)
  return
end if

! compute the explicit Hamiltonian
! (needs to be out here as we want to redetermine dynamically the selection vector)

call ExplH2(C,h0,TUVX,iSel,ExplE,ExplV)

! special case: nSel=nConf

if (nSel == nMaxSel) then
  if (IPRLEV >= DEBUG) write(u6,*) ' Initial CI-vectors are obtained by diagonalizing the explicit Hamiltonian'
  iDisk = IADR15(4)
  do i=1,lRoots
    C(1:nConf) = Zero
    do j=1,nSel
      k = iSel(j)
      C(k) = ExplV(j+(i-1)*nSel)
    end do
    call Save_tmp_CI_vec(i,nConf,C,LuDavid)
    if (IPRLEV >= INSANE) then
      write(String,'(A,I2)') 'CI vector of root',i
      write(String,'(A,I4,A)') '(max. ',nSel,' elements)'
      l = min(nSel,nConf)
      call dVcPrt(String,' ',ExplV(1+(i-1)*nSel),l)
    end if
  end do
  return
end if

if (Start_Vectors) then
  Start_Vectors = .false.
  if (ICIRST /= 0) then
#   ifdef _HDF5_
    call f_Inquire(StartOrbFile,Exists)
    if (Exists) then
      if (mh5_is_hdf5(StartOrbFile)) then
        if (IPRLEV >= TERSE) then
          write(u6,'(1x,a)') 'reading initial CI vectors from '//StartOrbFile
        end if
        mh5id = mh5_open_file_r(StartOrbFile)

        call mma_allocate(Tmp1,nConf,label='Scr1')
        call mma_allocate(vkcnf,nactel,label='kcnf')
        do i=1,lRoots
          call mh5_fetch_dset(mh5id,'CI_VECTORS',Tmp1,[nconf,1],[0,i-1])
          call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp1,C,vkcnf)
          call Save_CI_vec(i,nConf,C,LuDavid)
        end do
        call mma_deallocate(Tmp1)
        call mma_deallocate(vkcnf)

        call mh5_close_file(mh5id)
      else
        Exists = .false.
      end if
    end if
    if (.not. Exists) then
#   endif
      ! get start vectors by reading JOBOLD
      iJOB = 0
      call f_Inquire('JOBOLD',Exists)
      if (Exists) iJOB = 1
      if (iJOB == 1) then
        if (IPRLEV >= TERSE) then
          write(u6,'(1x,a)') 'reading initial CI vectors from JOBOLD'
        end if
        if (JOBOLD <= 0) then
          JOBOLD = 20
          call DaName(JOBOLD,'JOBOLD')
        end if
      else
        if (IPRLEV >= TERSE) then
          write(u6,'(1x,a)') 'reading initial CI vectors from JOBIPH'
        end if
        JOBOLD = JOBIPH
      end if

      iDisk = 0
      call IDafile(JOBOLD,2,iToc,15,iDisk)
      iDisk = iToc(4)
      call mma_allocate(Tmp1,nConf,label='Scr1')
      call mma_allocate(vkcnf,nactel,label='kcnf')
      do i=1,lRoots
        call DDafile(JOBOLD,2,Tmp1,nConf,iDisk)
        call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp1,C,vkcnf)
        call Save_CI_vec(i,nConf,C,LuDavid)
        if (IPRLEV >= INSANE) then
          write(String,'(A,I2)') 'Start vector of root',i
          write(String,'(A,I4,A)') '(max. ',nSel,' elements)'
          l = min(nSel,nConf)
          call dVcPrt(String,' ',C,l)
        end if
      end do
      call mma_deallocate(Tmp1)
      call mma_deallocate(vkcnf)
      if (iJOB == 1) then
        if ((JOBOLD > 0) .and. (JOBOLD /= JOBIPH)) then
          call DaClos(JOBOLD)
          JOBOLD = -1
        else if (JOBOLD > 0) then
          JOBOLD = -1
        end if
      end if
#   ifdef _HDF5_
    end if
#   endif

  else
    ! no CI restart, get start vectors by diagonalizing the explicit Hamiltonian
    if (IPRLEV >= DEBUG) write(u6,*) ' Initial CI-vectors are obtained by diagonalizing the explicit Hamiltonian'
    do i=1,lRoots
      C(1:nConf) = Zero
      do j=1,nSel
        k = iSel(j)
        C(k) = ExplV(j+(i-1)*nSel)
      end do
      call Save_CI_vec(i,nConf,C,LuDavid)
      if (IPRLEV >= INSANE) then
        write(String,'(A,I2)') 'Start vector of root',i
        write(String,'(A,I4,A)') '(max. ',nSel,' elements)'
        l = min(nSel,nConf)
        call dVcPrt(String,' ',C,l)
      end if
    end do

  end if

else
  ! no external start vector needed, get start vectors by reading JOBIPH
  if (iFinal == 2) then
    if (IPRLEV >= DEBUG) write(u6,*) ' Initial CI-vectors are identical to the transformed CI-vectors of the previous RASSCF &
                                     &iteration'
  else
    if (IPRLEV >= DEBUG) write(u6,*) ' Initial CI-vectors are identical to the CI-vectors of the previous RASSCF iteration'
  end if
  iDisk = IADR15(4)
  do i=1,lRoots-hRoots
    call DDafile(JOBIPH,2,C,nConf,iDisk)
    call Save_CI_vec(i,nConf,C,LuDavid)
    if (IPRLEV > 10) then
      write(String,'(A,I2)') 'Start vector of root',i
      write(String,'(A,I4,A)') '(max. ',nSel,' elements)'
      l = min(nSel,nConf)
      call dVcPrt(String,' ',C,l)
    end if
  end do
  ! MGD simple guess for missing ones : explV
  ! dangerous if linear dependence with converged states
  do i=lRoots-hRoots+1,lRoots
    C(1:nConf) = Zero
    do j=1,nSel
      k = iSel(j)
      C(k) = ExplV(j+(i-1)*nSel)
    end do
    call Save_CI_vec(i,nConf,C,LuDavid)
  end do

end if

return

end subroutine CStart_CI_Util
