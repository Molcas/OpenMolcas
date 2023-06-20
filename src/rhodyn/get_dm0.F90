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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine get_dm0()
!***********************************************************************
!
! Purpose : prepare the initial density matrix in CSF basis
!           based on user's input
!
!***********************************************************************

use rhodyn_data, only: CSF2SO, DM0, DTOC, E, E_SF, E_SO, flag_so, ipglob, k_B, lroots, lrootstot, N, N_Populated, nconf, nconftot, &
                       ndet, ndet_tot, NState, p_style, prep_dm_i, prep_dm_r, runmode, sint, T, U_CI
use rhodyn_utils, only: transform, dashes
use mh5, only: mh5_put_dset
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, j, k, ii, jj, kk, lu !temporary io unit
complex(kind=wp) :: Z
complex(kind=wp), allocatable :: DET2CSF(:,:), DM0_bas(:,:), temp_dm(:)
integer(kind=iwp), external :: isFreeUnit

if ((p_style == 'SO_THERMAL') .and. (T == 0)) then
  p_style = 'SO'
end if
if ((p_style == 'SF_THERMAL') .and. (T == 0)) then
  p_style = 'SF'
end if

if (runmode /= 4) then
  call mma_allocate(DM0,nconftot,nconftot,label='DM0')
else ! in this case prepare DM in SF/SO basis
  call mma_allocate(DM0,lrootstot,lrootstot,label='DM0')
end if
DM0 = cZero

! If populate DET, read the transform matrix from DET to CSF to DET2CSF
if (p_style == 'DET') then
  call mma_allocate(DET2CSF,ndet_tot,nconftot)
  DET2CSF = cZero
  ii = 0
  jj = 0
  kk = 0
  do i=1,N
    if (i == 1) then
      ii = 0
    else
      ii = nconf(i-1)
    end if
    do j=1,ndet(i)
      jj = jj+1
      do k=1,nconf(i)
        kk = ii+kk
        DET2CSF(jj,kk) = DTOC(j,k,i)
      end do
    end do
  end do
end if

! Construct the initial density matrix for populating

if (flag_so) then
  select case (p_style)
    case ('CSF')
      if ((N_Populated > nconftot) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated states',N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      end if
      if (runmode /= 4) then
        DM0(N_Populated,N_Populated) = cOne
      else
        call mma_allocate(DM0_bas,nconftot,nconftot)
        DM0_bas = cZero
        DM0_bas(N_Populated,N_Populated) = cOne
        write(u6,*) 'DM element ',N_populated,' set to 1'
        call transform(DM0_bas,CSF2SO,DM0)
      end if

    case ('DET')
      call mma_allocate(DM0_bas,ndet_tot,ndet_tot)
      DM0_bas = cZero
      if ((N_Populated > ndet_tot) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated states',N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      end if
      DM0_bas(N_Populated,N_Populated) = cOne
      call transform(DM0_bas,DET2CSF,DM0)
      call mma_deallocate(DET2CSF)

    case ('SF')
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas = cZero
      if ((N_Populated > lrootstot) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated is wrong'
        call dashes()
        call abend()
      end if
      DM0_bas(N_Populated,N_Populated) = cOne
      call transform(DM0_bas,cmplx(U_CI,kind=wp),DM0,.false.)

    case ('SO')
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas = cZero
      if ((N_Populated > lrootstot) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated states',N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      end if
      DM0_bas(N_Populated,N_Populated) = cOne
      call transform(DM0_bas,CSF2SO,DM0,.false.)

    case ('SO_THERMAL')
      if (ipglob > 3) then
        call dashes()
        write(u6,*) 'Printout the SO eigenvalues'
        write(u6,sint) 'Number of the SO states:',lrootstot
        call dashes()
        do i=1,10
          write(u6,*) E_SO(i)
        end do
      end if
      Z = cZero
      do i=1,lrootstot
        Z = Z+exp(-(E_SO(i)-E_SO(1))/(k_B*T))
      end do
      call mma_allocate(DM0_bas,lrootstot,lrootstot,label='DM0_bas')
      DM0_bas = cZero
      if (ipglob > 3) then
        call dashes()
        write(u6,*) 'I            E_SO(I)-E_SO(1)               exp(-(E_SO(I)-E_SO(1))/(k_B*T)),   '// &
                    'exp(-(E_SO(I)-E_SO(1))/(k_B*T))/Z'
        call dashes()
        do i=1,20
          write(u6,*) i,E_SO(i)-E_SO(1),exp(-(E_SO(i)-E_SO(1))/(k_B*T)),exp(-(E_SO(i)-E_SO(1))/(k_B*T))/Z
        end do
      end if
      do i=1,lrootstot
        DM0_bas(i,i) = exp(-(E_SO(i)-E_SO(1))/(k_B*T))/Z
      end do
      if (runmode /= 4) then
        call transform(DM0_bas,CSF2SO,DM0,.false.)
      else
        DM0(:,:) = DM0_bas
      end if
    case default
      write(u6,*) 'Population style ',p_style,' is not recognized'
      call abend()
  end select

else
  ! SO is off
  ! has not checked yet
  select case (p_style)
    case ('CSF')
      if ((N_Populated > nconftot) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated states',N_populated,'read from the input file is wrong'
        call dashes()
        call abend()
      end if
      if (runmode /= 4) then
        DM0(N_Populated,N_Populated) = cOne
      else
        call mma_allocate(DM0_bas,nconftot,nconftot)
        DM0_bas = cZero
        DM0_bas(N_Populated,N_Populated) = cOne
        write(u6,*) 'DM element ',N_populated,' set to 1'
        ! transform right to SF basis for propagation
        call transform(DM0_bas,cmplx(U_CI,kind=wp),DM0)
      end if

    case ('DET')
      call mma_allocate(DM0_bas,ndet_tot,ndet_tot)
      DM0_bas = cZero
      if ((N_Populated > NDET_TOT) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated states',N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      end if
      DM0_bas(N_Populated,N_Populated) = cOne
      call transform(DM0_bas,DET2CSF,DM0)
      call mma_deallocate(DET2CSF)

    case ('SF')
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas = cZero
      if ((N_Populated > lrootstot) .or. (N_populated <= 0)) then
        call dashes()
        write(u6,*) 'WARNING!!! Nr of the populated states',N_populated,' read from the input file is wrong'
        call dashes()
        call abend()
      end if
      DM0_bas(N_Populated,N_Populated) = cOne
      ! transform DM to CSF basis by default
      if (runmode /= 4) call transform(DM0_bas,cmplx(U_CI,kind=wp),DM0,.false.)

    case ('SO')
      call dashes()
      write(u6,*) 'WARNING! IF IFSO=OFF, CANNOT POPULATE SO STATE'
      call dashes()
      call abend()

    case ('SF_THERMAL')
      ! Construct the total eigenvalue of SF states to E_SF
      ii = 0
      do i=1,N
        do j=1,lroots(i)
          ii = ii+1
          E_SF(ii) = E(i,j)
        end do
      end do
      Z = cZero
      do i=1,nconftot
        Z = Z+exp(-(E_SF(i)-E_SF(1))/(k_B*T))
      end do
      call mma_allocate(DM0_bas,lrootstot,lrootstot)
      DM0_bas = cZero
      call dashes()
      write(u6,*) ' I               (E_SF(I)-E_SF(1))/(k_B*T)              exp(-(E_SF(I)-E_SF(1))/(k_B*T))  '// &
                  'exp(-(E_SF(I)-E_SF(1))/(k_B*T))/Z'
      call dashes()
      do i=1,nconftot
        write(u6,*) i,(E_SF(i)-E_SF(1))/(k_B*T),exp(-(E_SF(i)-E_SF(1))/(k_B*T)),exp(-(E_SF(i)-E_SF(1))/(k_B*T))/Z
        DM0_bas(i,i) = exp(-(E_SF(i)-E_SF(1))/(k_B*T))/Z
      end do
      ! transform DM to CSF basis by default
      if (runmode /= 4) call transform(DM0_bas,cmplx(U_CI,kind=wp),DM0,.false.)
    case default
      write(u6,*) 'Population style ',p_style,' is not recognized'
      call abend()
  end select
end if !ifso

! it will be transformed as usual further
if (p_style == 'FROMFILE') then
  call mma_allocate(temp_dm,NState)
  lu = isFreeUnit(11)
  call molcas_open(lu,'INDENS')
  do i=1,NState
    do j=1,Nstate
      read(lu,'(E16.8)',advance='no') temp_dm(j)
    end do
    DM0(i,:) = temp_dm
    read(lu,*)
  end do
  close(lu) ! close INDENS file
  call mma_deallocate(temp_dm)
end if

if (allocated(DM0_bas)) call mma_deallocate(DM0_bas)

! Printout the initial density matrix to SDPREP in CSF basis
! not CM case
if (runmode /= 4) then
  if (ipglob > 3) then
    call dashes()
    write(u6,*) 'The initial density matrix is saved in RDPREP'
    call dashes()
  end if
  call mh5_put_dset(prep_dm_r,real(DM0))
  call mh5_put_dset(prep_dm_i,aimag(DM0))
end if

end subroutine get_dm0
