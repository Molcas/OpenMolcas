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
! Copyright (C) 1998, Markus P. Fuelscher                              *
!***********************************************************************
!  readvc_m
!
!> @brief get start MO coeffients
!>
!> @author M. P. Fuelscher
!> @author Matthew R. Hennefarth
!>
!> @param[out] cmo mo coefficients
!***********************************************************************

subroutine ReadVC_m(CMO)

use PrintLevel, only: DEBUG, TERSE, VERBOSE
use mcpdft_output, only: iPrGlb, iPrLoc
use mcpdft_input, only: mcpdft_options
use general_data, only: invec, jobiph, jobold, ntot2
#ifdef _HDF5_
use mh5, only: mh5_close_file, mh5_fetch_dset, mh5_open_file_r
#endif
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: CMO(ntot2)
integer(kind=iwp) :: iad19, IADR19(30), ijob, iprlev
#ifdef _HDF5_
integer(kind=iwp) ::mh5id
#endif
logical(kind=iwp) :: Found

IPRLEV = IPRLOC(1)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering READVC'

! invec can either be 3 or 4 based off of proc_inpx
! can also be 5 (if FileOrb points to a non hdf5 reference file)

! read from unit formatted ascii file with starting orbitals

! Note: Inside RDVEC, the file wfn_file is opened, but uses blindly
! the unit number provided here. So that should better be a usable
! number, or else!
! read from unit JOBOLD (binary file)
if (InVec == 3) then
  IAD19 = 0
  iJOB = 0
  call f_Inquire('JOBOLD',Found)
  if (Found) iJOB = 1
  if (iJOB == 1) then
    if (JOBOLD <= 0) then
      JOBOLD = 20
      call DaName(JOBOLD,'JOBOLD')
    end if
  else
    if (IPRLEV >= TERSE) write(u6,*) '  File JOBOLD not found -- use JOBIPH.'
    if (JOBIPH > 0) then
      JOBOLD = JOBIPH
    else
      call DaName(JOBOLD,mcpdft_options%wfn_file)
    end if
  end if
  call IDaFile(JOBOLD,2,IADR19,15,IAD19)
  if (IADR19(15) == -1) then
    IAD19 = 0
    call IDAFILE(JOBOLD,2,IADR19,30,IAD19)
  else
    IADR19(16:30) = 0
    if (IPRGLB >= VERBOSE) call WarningMessage(1,'Old JOBIP file layout.')
  end if
  if (IPRLEV >= TERSE) then
    if (iJOB == 1) then
      write(u6,'(6X,A)') 'The MO-coefficients are taken from the file:'
      write(u6,'(6X,A)') 'JOBOLD'
    else
      write(u6,'(6X,A)') 'The MO-coefficients are taken from the file:'
      write(u6,'(6X,A)') trim(mcpdft_options%wfn_file)
    end if
  end if

  iAd19 = iAdr19(2)
  call DDaFile(JobOld,2,CMO,NTOT2,iAd19)

  if ((JOBOLD > 0) .and. (JOBOLD /= JOBIPH)) then
    call DaClos(JOBOLD)
    JOBOLD = -1
  else if (JOBOLD > 0) then
    JOBOLD = -1
  end if

else if (InVec == 4) then
  ! read from a HDF5 wavefunction file
# ifdef _HDF5_
  if (IPRLEV >= TERSE) then
    write(u6,'(6X,A)') 'The MO-coefficients are taken from the file:'
    write(u6,'(6X,A)') trim(mcpdft_options%wfn_file)
  end if

  mh5id = mh5_open_file_r(mcpdft_options%wfn_file)
  call mh5_fetch_dset(mh5id,'MO_VECTORS',CMO)
  call mh5_close_file(mh5id)
# else
  write(u6,*) 'Orbitals requested from HDF5, but this'
  write(u6,*) 'installation does not support that, abort!'
  call abend()
# endif
else if (invec == 5) then
  write(u6,*) 'FileOrb specified, but does not point to hdf5 file'
  write(u6,*) 'This has not been implemented, aborting'
  call abend()
end if

end subroutine ReadVC_m
