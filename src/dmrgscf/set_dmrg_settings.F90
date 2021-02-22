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
use active_space_solver_cfg, only: as_solver, as_solver_inp_proc
use Definitions, only: iwp, u6

implicit none
! ----------------------------------------------------------------------
integer(kind=iwp)  :: luspool, irc, nr_lines
character(len=180) :: line

character(len=180), external :: get_ln
integer(kind=iwp), external  :: isfreeunit
! ----------------------------------------------------------------------

nr_lines = 0
luspool = isfreeunit(18)
call spoolinp(luspool)

rewind(luspool)
call rdnlst(luspool,'DMRGSCF')

call setpos(luspool,'DMRG',line,irc)

if (irc /= 0) then
  call warningmessage(2,'Error in input processing.')
  write(u6,*) ' SET_DMRG_SETTUNGS: active space solver settings are missing but'
  write(u6,*) ' are required for DMRG calculations in OpenMOLCAS to set compulsory'
  write(u6,*) ' DMRG internal parameters, for example:'
  write(u6,*) ' max_bond_dimension (aka m), nsweeps, ...'
  write(u6,*) ' Please consult the manual for the active space solver ',as_solver(1:8),' for further details.'
  iRc = 112
  call abend()
end if

if (as_solver(1:8) == 'qcmaquis') then
# ifdef _DMRG_
  call qcmaquis_rdinp(luspool,1,nr_lines)
  call setpos(luspool,'DMRG',line,irc)
  call qcmaquis_rdinp(luspool,2,nr_lines)
# endif
  as_solver_inp_proc = .true.
end if

call close_luspool(luspool)

end subroutine set_dmrg_settings
! ----------------------------------------------------------------------
