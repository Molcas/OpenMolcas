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
! Copyright (C) 2017, Stefan Knecht                                    *
!***********************************************************************
  subroutine set_dmrg_settings()

! module dependencies
  use active_space_solver_cfg

  implicit none

! ----------------------------------------------------------------------
  integer            :: luspool, irc, nr_lines
  character(len=180) :: line, blank

  external           :: get_ln, isfreeunit
  character(len=180) :: get_ln
  integer            :: isfreeunit
! ----------------------------------------------------------------------


  nr_lines = 0
  luspool = 18
  luspool = isfreeunit(luspool)
  call spoolinp(luspool)

  rewind(luspool)
  call rdnlst(luspool,'DMRGSCF')
  blank=' '

  call setpos(luspool,'DMRG',line,irc)

  if(irc /= 0)then
    call warningmessage(2,'Error in input processing.')
    write(6,*)' SET_DMRG_SETTUNGS: active space solver settings are missing but'
    write(6,*)' are required for DMRG calculations in OpenMOLCAS to set compulsory'
    write(6,*)' DMRG internal parameters, for example:'
    write(6,*)' max_bond_dimension (aka m), nsweeps, ...'
    write(6,*)' Please consult the manual for the active space solver ',as_solver(1:8),' for further details.'
    iRc=112
    call abend()
  end if

  if(as_solver(1:8) == 'qcmaquis')then
#ifdef _DMRG_
    call qcmaquis_rdinp(luspool,1,nr_lines)
    call setpos(luspool,'DMRG',line,irc)
    call qcmaquis_rdinp(luspool,2,nr_lines)
#endif
    as_solver_inp_proc = .true.
  end if

  call close_luspool(luspool)

  end subroutine set_dmrg_settings
! ----------------------------------------------------------------------
