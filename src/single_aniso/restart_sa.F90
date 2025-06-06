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

subroutine restart_sa(input_to_read,input_file_name,nss,nstate)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: input_to_read
character(len=180), intent(in) :: input_file_name
integer(kind=iwp), intent(inout) :: nss, nstate
integer(kind=iwp) :: iDisk, idum(1), luaniso
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: dbg = _DBG_
integer(kind=iwp), external :: IsFreeUnit

if (input_to_read == 1) then
  ! read the binary file "$Project.aniso":
  luaniso = 8
  call daname(luaniso,'POLYFILE')
  iDisk = 0
  call idafile(luaniso,2,idum,1,iDisk)
  nstate = idum(1)
  call idafile(luaniso,2,idum,1,iDisk)
  nss = idum(1)
  call daclos(luaniso)
  ! put them on RunFile:
  call Put_iScalar('NSTATE_SINGLE   ',nstate)
  call Put_iScalar('NSS_SINGLE      ',nss)
  call Put_iScalar('MXJOB_SINGLE    ',1)
  call Put_iScalar('NJOB_SINGLE     ',1)

else if ((input_to_read == 2) .or. (input_to_read == 4)) then
  ! read the ascii formatted "aniso.input" file:
  luaniso = IsFreeUnit(18)
  call molcas_open(luaniso,input_file_name)
  read(luaniso,*) nstate,nss
  close(luaniso)
  ! put them on RunFile:
  call Put_iScalar('NSTATE_SINGLE   ',nstate)
  call Put_iScalar('NSS_SINGLE      ',nss)
  call Put_iScalar('MXJOB_SINGLE    ',1)
  call Put_iScalar('NJOB_SINGLE     ',1)

else if (input_to_read == 3) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'restart_sa: file h5=',trim(input_file_name)
# endif
# ifdef _HDF5_
  ! NSS and NSTATE are also placed on RunFile
  call read_hdf5_init(input_file_name,nstate,nss)
# ifdef _DEBUGPRINT_
  write(u6,*) 'restart_sa:    nss=',nss
  write(u6,*) 'restart_sa: nstate=',nstate
# endif
  call Put_iScalar('NSTATE_SINGLE   ',nstate)
  call Put_iScalar('NSS_SINGLE      ',nss)
  call Put_iScalar('MXJOB_SINGLE    ',1)
  call Put_iScalar('NJOB_SINGLE     ',1)
# else
  write(u6,'(A)') 'Warning:: restart option was set to 3: i.e. from an HDF5 file'
  write(u6,'(A,A)') 'file id =',trim(input_file_name)
  call WarningMessage(2,'MOLCAS was compiled without _HDF5_ option.')
  call Quit_OnUserError()
# endif

else if (input_to_read == 6) then
  ! read the ascii formatted "aniso.input" file (NEW):
  luaniso = IsFreeUnit(18)
  call molcas_open(luaniso,input_file_name)

  call read_nss(luaniso,nss,dbg)
  call read_nstate(luaniso,nstate,dbg)

  ! put them on RunFile:
  call Put_iScalar('NSTATE_SINGLE   ',nstate)
  call Put_iScalar('NSS_SINGLE      ',nss)
  call Put_iScalar('MXJOB_SINGLE    ',1)
  call Put_iScalar('NJOB_SINGLE     ',1)

  close(luaniso)

else
  call WarningMessage(2,'SINGLE_ANISO:: RESTART  option is not known.')
  write(u6,'(A,I6)') 'restart_option =',input_to_read
  write(u6,'(A,I6)') 'restart_option can only take integer values:'
  write(u6,'(A,I6)') '1 - from binary $Project.aniso'
  write(u6,'(A,I6)') '2 - from formatted file "aniso.input" (filename can be given in the input)'
  write(u6,'(A,I6)') '3 - from an HDF5 type file generated by RASSI code (filename can be given in the input)'
  write(u6,'(A,I6)') '4 - from formatted file "aniso.input" (filename can be given in the input) in molcas-8.0 format'
  call Quit_OnUserError()

end if

return

end subroutine restart_sa
