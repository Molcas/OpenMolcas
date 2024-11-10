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

!> @brief get start MO coeffients
!>
!> @author M. P. Fuelscher
!> @author Matthew R. Hennefarth
!>
!> @param[out] cmo mo coefficients
Subroutine ReadVC_m(CMO)
  use definitions,only:wp,iwp,u6
  use printlevel,only:terse,verbose,debug
  use mcpdft_output,only:iPrGlb,iPrLoc
  use mcpdft_input,only:mcpdft_options
  use general_data,only:invec,jobiph,jobold,ntot2
#ifdef _HDF5_
  use mh5,only:mh5_open_file_r,mh5_fetch_dset,mh5_close_file
#endif
  implicit none

  real(kind=wp),intent(out) :: CMO(*)
  logical(kind=iwp) :: Found
  integer(kind=iwp) :: IADR19(30),i,iad19,ijob,iprlev

#ifdef _HDF5_
  integer(kind=iwp) mh5id
#endif

  IPRLEV = IPRLOC(1)
  IF(IPRLEV >= DEBUG) THEN
    WRITE(u6,*) ' Entering READVC'
  ENDIF

  ! invec can either be 3 or 4 based off of proc_inpx
  ! can also be 5 (if FileOrb points to a non hdf5 reference file)

  ! read from unit formatted ascii file with starting orbitals

  ! Note: Inside RDVEC, the file wfn_file is opened, but uses blindly
  ! the unit number provided here. So that should better be a usable
  ! number, or else!
  ! read from unit JOBOLD (binary file)
  If(InVec == 3) then
    IAD19 = 0
    iJOB = 0
    Call f_Inquire('JOBOLD',Found)
    If(Found) iJOB = 1
    If(iJOB == 1) Then
      if(JOBOLD <= 0) Then
        JOBOLD = 20
        Call DaName(JOBOLD,'JOBOLD')
      endif
    Else
      If(IPRLEV >= TERSE) then
        Write(u6,*) '  File JOBOLD not found -- use JOBIPH.'
      EndIf
      If(JOBIPH > 0) Then
        JOBOLD = JOBIPH
      Else
        Call DaName(JOBOLD,mcpdft_options%wfn_file)
      EndIf
    EndIf
    Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
    IF(IADR19(15) == -1) THEN
      IAD19 = 0
      CALL IDAFILE(JOBOLD,2,IADR19,30,IAD19)
    ELSE
      DO I = 16,30
        IADR19(I) = 0
      ENDDO
      IF(IPRGLB >= VERBOSE) then

        Call WarningMessage(1,'Old JOBIP file layout.')
      endif
    ENDIF
    IF(IPRLEV >= TERSE) THEN
      If(iJOB == 1) Then
        Write(u6,'(6X,A)') 'The MO-coefficients are taken from the file:'
        Write(u6,'(6X,A)') 'JOBOLD'
      Else
        Write(u6,'(6X,A)') 'The MO-coefficients are taken from the file:'
        Write(u6,'(6X,A)') trim(mcpdft_options%wfn_file)
      EndIf
    ENDIF

    iAd19 = iAdr19(2)
    Call DDaFile(JobOld,2,CMO,NTOT2,iAd19)

    If(JOBOLD > 0 .and. JOBOLD /= JOBIPH) Then
      Call DaClos(JOBOLD)
      JOBOLD = -1
    Else If(JOBOLD > 0) Then
      JOBOLD = -1
    EndIf

    ! read from a HDF5 wavefunction file
  Else If(InVec == 4) then
#ifdef _HDF5_
    IF(IPRLEV >= TERSE) THEN
      Write(u6,'(6X,A)') 'The MO-coefficients are taken from the file:'
      Write(u6,'(6X,A)') trim(mcpdft_options%wfn_file)
    ENDIF

    mh5id = mh5_open_file_r(mcpdft_options%wfn_file)
    call mh5_fetch_dset(mh5id,'MO_VECTORS',CMO)
    call mh5_close_file(mh5id)
#else
    write(u6,*) 'Orbitals requested from HDF5, but this'
    write(u6,*) 'installation does not support that, abort!'
    call abend
#endif
  else if(invec == 5) then
    write(u6,*) "FileOrb specified, but does not point to hdf5 file"
    write(u6,*) "This has not been implemented, aborting"
    call abend
  EndIf

ENDSubroutine ReadVC_m
