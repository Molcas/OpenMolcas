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

subroutine pop(time,popcount,dgl_csf,density_csf)
!***********************************************************************
! prints diagonal of the density matrix densityt in required basis
! at the current time
!***********************************************************************

use rhodyn_data, only: a_einstein, basis, CSF2SO, d, density0, densityt, dgl, dipole_basis, DM_basis, emiss, flag_dipole, &
                       flag_emiss, lu_csf, lu_dip, lu_sf, lu_so, n_freq, nconftot, Nstate, out_dm_csf, out_dm_sf, out_dm_so, &
                       out_emiss, out_fmt, out_fmt_csf, out_tout, pulse_vec, SO_CI, tmp, U_CI_compl
use rhodyn_utils, only: transform
use linalg_mod, only: mult
use mh5, only: mh5_put_dset
use Constants, only: Zero, auToFs
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: time
integer(kind=iwp), intent(in) :: popcount
real(kind=wp), intent(out) :: dgl_csf(nconftot)
complex(kind=wp), intent(out) :: density_csf(nconftot,nconftot)
integer(kind=iwp) :: i, j, l
real(kind=wp) :: norm
character(len=64) :: sline

! here notation d is dimension of all the basis matrices
!!! density0 (can't be) used as a temporary storage for dm in required basis

write(sline,'(f10.3)') time*auToFs
call StatusLine('RhoDyn: current time ',trim(sline))

call mh5_put_dset(out_tout,[time*auToFs],[1],[popcount-1])

if (basis == 'CSF') then

  if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
    ! the density in CSF basis
    dgl_csf = [(real(densityt(i,i)),i=1,d)]
    norm = sum(dgl_csf)
    write(lu_csf,out_fmt) time*auToFs,(dgl_csf(i),i=1,d),norm
    call mh5_put_dset(out_dm_csf,dgl_csf,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from CSF to SF states basis
    call transform(densityt,U_CI_compl,density0)
    dgl(:) = [(real(density0(i,i)),i=1,d)]
    norm = sum(dgl)
    write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_sf,dgl,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from CSF to SO states basis
    call transform(densityt,CSF2SO,density0)
    dgl(:) = [(real(density0(i,i)),i=1,d)]
    norm = sum(dgl)
    write(lu_so,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_so,dgl,[1,d],[popcount-1,0])
  end if

else if (basis == 'SF') then

  if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SF to CSF states basis
    call transform(densityt,U_CI_compl,density_csf,.false.)
    dgl_csf = [(real(density_csf(i,i)),i=1,nconftot)]
    norm = sum(dgl_csf)
    write(lu_csf,out_fmt_csf) time*auToFs,(dgl_csf(i),i=1,nconftot),norm
    call mh5_put_dset(out_dm_csf,dgl_csf,[1,nconftot],[popcount-1,0])
  end if

  if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! the density in SF basis
    dgl(:) = [(real(densityt(i,i)),i=1,Nstate)]
    norm = sum(dgl)
    write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,Nstate),norm
    call mh5_put_dset(out_dm_sf,dgl,[1,Nstate],[popcount-1,0])
  end if

  if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SF to SO states basis
    call transform(densityt,SO_CI,density0)
    dgl(:) = [(real(density0(i,i)),i=1,d)]
    norm = sum(dgl)
    write(lu_so,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_so,dgl,[1,d],[popcount-1,0])
  end if

else if (basis == 'SO') then

  if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SO to CSF states basis
    call transform(densityt,CSF2SO,density_csf,.false.)
    dgl_csf = [(real(density_csf(i,i)),i=1,nconftot)]
    norm = sum(dgl_csf)
    write(lu_csf,out_fmt_csf) time*auToFs,(dgl_csf(i),i=1,nconftot),norm
    call mh5_put_dset(out_dm_csf,dgl_csf,[1,nconftot],[popcount-1,0])
  end if

  if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! transform the density from SO to SF states basis
    call transform(densityt,SO_CI,density0,.false.)
    dgl(:) = [(real(density0(i,i)),i=1,d)]
    norm = sum(dgl)
    write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_sf,dgl,[1,d],[popcount-1,0])
  end if

  if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
    ! the density in SO basis
    dgl(:) = [(real(densityt(i,i)),i=1,d)]
    norm = sum(dgl)
    write(lu_so,out_fmt) time*auToFs,(dgl(i),i=1,d),norm
    call mh5_put_dset(out_dm_so,dgl,[1,d],[popcount-1,0])
  end if

else if (basis == 'SPH') then
  dgl(:) = [(real(densityt(i,i)),i=1,Nstate)]
  norm = sum(dgl)
  write(lu_sf,out_fmt) time*auToFs,(dgl(i),i=1,Nstate),norm
end if

! time-dependent dipole moment
! pulse_vec is used as storage for tr(rho d)
if (flag_dipole) then
  do i=1,3
    call mult(densityt,dipole_basis(:,:,i),tmp)
    dgl(:) = [(real(tmp(i,i)),i=1,d)]
    pulse_vec(i) = sum(dgl)
  end do
  write(lu_dip,'(7(g25.15e3,2x))') time*auToFs,(real(pulse_vec(i)),aimag(pulse_vec(i)),i=1,3)
end if

! emission spectra
if (flag_emiss) then
  ! in dgl the SOC populations are left
  l = 1
  emiss = Zero
  do j=1,(Nstate-1)
    do i=(j+1),Nstate
      emiss(l) = a_einstein(i,j)*dgl(i)
      l = l+1
    end do
  end do
  call mh5_put_dset(out_emiss,emiss,[1,n_freq],[popcount-1,0])
end if

end subroutine pop
