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

subroutine cre_out()
!***********************************************************************
! creates output h5 file *.rhodyn.h5 as well as ASCII files for
! densities, pulse, dipole moments if needed
!***********************************************************************

use rhodyn_data, only: basis, d, DM_basis, flag_decay, flag_dipole, flag_emiss, flag_fdm, flag_pulse, len_sph, lu_csf, lu_dip, &
                       lu_pls, lu_sf, lu_so, n_freq, nconftot, Npop, Nstate, Nstep, Ntime_tmp_dm, out2_fmt, out3_fmt, out_decay_i, &
                       out_decay_r, out_dm_csf, out_dm_sf, out_dm_so, out_emiss, out_fdmi, out_fdmr, out_fmt, out_fmt_csf, &
                       out_freq, out_ham_i, out_ham_r, out_id, out_pulse, out_t, out_tfdm, out_tout
use mh5, only: mh5_create_dset_real, mh5_create_file, mh5_init_attr
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i
integer(kind=iwp), external :: isfreeunit
character(len=64) :: out1_fmt, out1_fmt_csf

! write formats for output files SFDENS, SODENS, CSFDENS
! header formats
write(out1_fmt,'(a,i5,a)') '(a22,',Nstate,'(i8,14x),a)'
write(out1_fmt_csf,'(a,i5,a)') '(a22,',nconftot,'(i8,14x),a)'
! line formats
write(out_fmt,'(a,i5,a)') '(1x,',Nstate+2,'(f22.16))'
write(out_fmt_csf,'(a,i5,a)') '(1x,',nconftot+2,'(f22.16))'
out2_fmt = '(2x,a,28x,a,28x,a,28x,a)'
out3_fmt = '(2x,f22.16,1x,f22.16,1x,f22.16,1x,i1,a1,i2.2,a1,i2.2)'

! creating main output file of rhodyn
out_id = mh5_create_file('RDOUT')

call mh5_init_attr(out_id,'MOLCAS_MODULE','RHODYN')

! pulse
if (flag_pulse) then
  out_pulse = mh5_create_dset_real(out_id,'PULSE',2,[Nstep,6])
  call mh5_init_attr(out_pulse,'description','Pulse')
  ! ascii file
  lu_pls = isfreeunit(lu_pls)
  call molcas_open(lu_pls,'PULSE')
end if

! time
out_t = mh5_create_dset_real(out_id,'TIME',1,[Nstep])
call mh5_init_attr(out_t,'description','Complete time grid')

! this set/file is for the diagonal density matrix in CSF basis
if ((DM_basis == 'CSF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'ALL')) then
  out_dm_csf = mh5_create_dset_real(out_id,'DM_CSF',2,[Npop,nconftot])
  call mh5_init_attr(out_dm_csf,'description','Density matrix in CSF basis')
  ! ascii file
  lu_csf = isfreeunit(lu_csf)
  call molcas_open(lu_csf,'CSFDEN')
  write(lu_csf,out1_fmt_csf) '#time(fs)',(i,i=1,nconftot),'Norm'
end if

! this set/file is for the diagonal of density matrix in SO basis
if ((DM_basis == 'SO') .or. (DM_basis == 'CSF_SO') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
  out_dm_so = mh5_create_dset_real(out_id,'DM_SO',2,[Npop,Nstate])
  call mh5_init_attr(out_dm_so,'description','Density matrix in SO basis')
  ! ascii file
  lu_so = isfreeunit(lu_so)
  call molcas_open(lu_so,'SODENS')
  write(lu_so,out1_fmt) '#time(fs)',(i,i=1,Nstate),'Norm'
end if

! this set/file is for the diagonal density matrix in SF basis
if ((DM_basis == 'SF') .or. (DM_basis == 'CSF_SF') .or. (DM_basis == 'SF_SO') .or. (DM_basis == 'ALL')) then
  out_dm_sf = mh5_create_dset_real(out_id,'DM_SF',2,[Npop,Nstate])
  call mh5_init_attr(out_dm_sf,'description','Density matrix in SF basis')
  ! ascii file
  lu_sf = isfreeunit(lu_sf)
  call molcas_open(lu_sf,'SFDENS')
  write(lu_sf,out1_fmt) '#time(fs)',(i,i=1,Nstate),'Norm'
end if

! time points for stored populations
out_tout = mh5_create_dset_real(out_id,'TOUT',1,[Npop])
call mh5_init_attr(out_tout,'description','TOUT step time grid')

! Hamiltonian used for propagation
out_ham_r = mh5_create_dset_real(out_id,'HAM_R',2,[d,d])
call mh5_init_attr(out_ham_r,'description','Hamiltonian used for propagation, real part')
out_ham_i = mh5_create_dset_real(out_id,'HAM_I',2,[d,d])
call mh5_init_attr(out_ham_i,'description','Hamiltonian used for propagation, imaginary part')

! decay matrix
if (flag_decay) then
  out_decay_r = mh5_create_dset_real(out_id,'DECAY_R',2,[d,d])
  call mh5_init_attr(out_decay_r,'description','Decay matrix, real part')
  out_decay_i = mh5_create_dset_real(out_id,'DECAY_I',2,[d,d])
  call mh5_init_attr(out_decay_i,'description','Decay matrix, imaginary part')
end if

if (flag_emiss) then
  ! frequencies
  out_freq = mh5_create_dset_real(out_id,'FREQ',1,[n_freq])
  call mh5_init_attr(out_freq,'description','frequencies')
  ! emission spectra
  out_emiss = mh5_create_dset_real(out_id,'EMISSION',2,[Npop,n_freq])
  call mh5_init_attr(out_emiss,'description','emission spectrum')
end if

if (flag_fdm) then
  ! time points for stored full density matrix
  out_tfdm = mh5_create_dset_real(out_id,'TFDM',1,[Ntime_tmp_dm])
  call mh5_init_attr(out_tfdm,'description','TFDM grid')
  ! full density matrix (only abs values for state multipoles in spherical basis)
  if (basis == 'SPH') then
    out_fdmr = mh5_create_dset_real(out_id,'SPH_MULTIPOLES',4,[Ntime_tmp_dm,len_sph,d,d])
    call mh5_init_attr(out_fdmr,'description','Absolute values of state multipoles stored as [n_time,sph_length,d,d]')
  else
    out_fdmr = mh5_create_dset_real(out_id,'FDM_R',3,[Ntime_tmp_dm,d,d])
    out_fdmi = mh5_create_dset_real(out_id,'FDM_I',3,[Ntime_tmp_dm,d,d])
    call mh5_init_attr(out_fdmr,'description','Real part of the full density matrix')
    call mh5_init_attr(out_fdmi,'description','Imaginary part of the full density matrix')
  end if
end if

! this file is for the TD-dipole moment data
if (flag_dipole) then
  lu_dip = isfreeunit(lu_dip)
  call molcas_open(lu_dip,'DIPOLE')
end if

end subroutine cre_out
