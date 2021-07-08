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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine cre_out

  use rhodyn_data
  use mh5
  implicit none

! creating main output file of rhodyn
  out_id = mh5_create_file('RDOUT')  

  call mh5_init_attr (out_id,'MOLCAS_MODULE','RHODYN')  
! PULSE
  out_pulse = mh5_create_dset_real (out_id,'PULSE',2,[Nstep,6])
  call mh5_init_attr(out_pulse,'description','Pulse')
! TIME
  out_t = mh5_create_dset_real (out_id,'TIME',1,[Nstep])
  call mh5_init_attr(out_t,'description','Complete time grid')
! density in csf basis
  out_dm_csf = mh5_create_dset_real(out_id, &
        'DM_CSF', 2, [Npop,nconftot])
  call mh5_init_attr(out_dm_csf, 'description', &
        'Density matrix in CSF basis')
! density in so basis
  out_dm_so = mh5_create_dset_real (out_id, &
        'DM_SO', 2, [Npop,Nstate])
  call mh5_init_attr(out_dm_so, 'description', &
        'Density matrix in SO basis')
! density in spin free sf basis
  out_dm_sf = mh5_create_dset_real(out_id,'DM_SF',2,[Npop,Nstate])
  call mh5_init_attr(out_dm_sf, 'description', &
        'Density matrix in SF basis')
! TIME FOR DENSITY OUT
  out_tout = mh5_create_dset_real (out_id,'TOUT',1,[Npop])
  call mh5_init_attr(out_tout,'description','TOUT step time grid')
! Hamiltonian used for propagation
  out_ham_r = mh5_create_dset_real (out_id, 'HAM_R', 2, [d,d])
  call mh5_init_attr(out_ham_r, 'description', &
         'Hamiltonian used for propagation, real part')
  out_ham_i = mh5_create_dset_real (out_id, 'HAM_I', 2, [d,d])
  call mh5_init_attr(out_ham_i, 'description', &
         'Hamiltonian used for propagation, imaginary part')
! Decay matrix
  out_decay_r = mh5_create_dset_real (out_id, 'DECAY_R', 2, [d,d])
  call mh5_init_attr(out_decay_r, 'description', &
         'Decay matrix, real part')
  out_decay_i = mh5_create_dset_real (out_id, 'DECAY_I', 2, [d,d])
  call mh5_init_attr(out_decay_i, 'description', &
         'Decay matrix, imaginary part')
! frequencies
  out_freq = mh5_create_dset_real(out_id,'FREQ',1,[n_freq])
  call mh5_init_attr(out_freq,'description','frequencies')
! emission spectra
  out_emiss=mh5_create_dset_real(out_id,'EMIS',2,[Npop,n_freq])
  call mh5_init_attr(out_emiss,'description','emission spec')
! TIME steps FOR FULL DENSITY OUT
  out_tfdm = mh5_create_dset_real (out_id,'TFDM',1,[Ntime_tmp_dm])
  call mh5_init_attr(out_tfdm,'description','TFDM grid')
! Full density matrix out
  out_fdm = mh5_create_dset_real (out_id,'FDM',3,[Ntime_tmp_dm,d,d])
  call mh5_init_attr(out_fdm,'description','Full density matrix ABS')
end
