************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Markus P. Fuelscher                              *
************************************************************************
      Subroutine ReadVC_m(CMO,OCC,D,DS,P,PA)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     get start MO-coefficients                                        *
*     (if the CI is restarted from a previous calculation get          *
*      also the density matrices)                                      *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : array of real*8                                        *
*               MO-coefficients                                        *
*     OCC     : array of real*8                                        *
*               occupation numbers                                     *
*     INVEC   : integer                                                *
*               flag indicating orbital type                           *
*     D       : array of real*8                                        *
*               averaged one-body density matrix                       *
*     DS      : array of real*8                                        *
*               averaged one-body spin density matrix                  *
*     P       : array of real*8                                        *
*               averaged two body density matrix                       *
*     PA      : array of real*8                                        *
*               averaged antisymmetric twobody density matrix          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1998                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     none                                                             *
*                                                                      *
************************************************************************

#ifdef _HDF5_
      use mh5, only: mh5_open_file_r, mh5_fetch_dset, mh5_close_file
#endif
      use definitions, only: wp, iwp
      use printlevel, only: terse, verbose, debug
      use mcpdft_output, only: lf, iPrGlb, iPrLoc
      use mcpdft_input, only: mcpdft_options
      use stdalloc, only: mma_allocate, mma_deallocate

      implicit none

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "SysDef.fh"
#include "warnings.h"

      real(kind=wp), Dimension(*) :: CMO,OCC,D,DS,P,PA

!     local data declarations
      integer(kind=iwp), dimension(30) :: IADR19
      logical :: Found
      integer(kind=iwp) :: i, iad19, idisk, ijob, iprlev
      integer(kind=iwp) :: jroot, kroot

      real(kind=wp) :: scal
      real(kind=wp), Allocatable :: Scr(:), Ene(:)

#ifdef _HDF5_
      integer mh5id
#endif

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
! Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV >= DEBUG) THEN
        WRITE(LF,*)' Entering READVC'
      END IF

! invec can either be 3 or 4 based off of proc_inpx
! can also be 5 (if FileOrb points to a non hdf5 reference file)

! read from unit formatted ascii file with starting orbitals

! Note: Inside RDVEC, the file wfn_file is opened, but uses blindly
! the unit number provided here. So that should better be a usable
! number, or else!
!     read from unit JOBOLD (binary file)
      If ( InVec == 3 ) then
        IAD19=0
        iJOB=0
        Call f_Inquire('JOBOLD',Found)
        If (Found) iJOB=1
        If (iJOB.eq.1) Then
           if(JOBOLD <= 0) Then
             JOBOLD=20
             Call DaName(JOBOLD,'JOBOLD')
           end if
        Else
           If (IPRLEV >= TERSE) then
              Write(LF,*) '  File JOBOLD not found -- use JOBIPH.'
           End If
           If (JOBIPH > 0) Then
              JOBOLD=JOBIPH
           Else
              Call DaName(JOBOLD, mcpdft_options%wfn_file)
           End If
        End If
        Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
        IF(IADR19(15) == -1) THEN
          IAD19=0
          CALL IDAFILE(JOBOLD,2,IADR19,30,IAD19)
        ELSE
          DO I=16,30
            IADR19(I)=0
          END DO
          IF(IPRGLB >= VERBOSE)
     &               Call WarningMessage(1,'Old JOBIP file layout.')
        END IF
        IF(IPRLEV >= TERSE) THEN
         If (iJOB == 1) Then
            Write(LF,'(6X,A)')
     &      'The MO-coefficients are taken from the file:'
            Write(LF,'(6X,A)') 'JOBOLD'
         Else
            Write(LF,'(6X,A)')
     &      'The MO-coefficients are taken from the file:'
            Write(LF,'(6X,A)') trim(mcpdft_options%wfn_file)
         End If
        END IF

        iAd19=iAdr19(2)
        Call DDaFile(JobOld,2,CMO,NTOT2,iAd19)
        Call DDaFile(JobOld,2,OCC,nTot,iAd19)
        If ( IPRLEV >= VERBOSE) then
          If (iJOB == 1) Then
             Write(LF,'(6X,A)')
     &       'The active density matrices (D,DS,P,PA) are read from'//
     &       ' file JOBOLD and weighted together.'
          Else
             Write(LF,'(6X,A)')
     &       'The active density matrices (D,DS,P,PA) are read from'//
     &       ' file '//trim(mcpdft_options%wfn_file)//
     &       ' and weighted together.'
          End If
        End If
        Call mma_allocate(scr,NACPR2,Label='Scr')
        iDisk = IADR19(3)
        Do jRoot = 1,lRoots
          scal = 0.0d0
          Do kRoot = 1,nRoots
            If ( iRoot(kRoot).eq.jRoot ) then
              scal = Weight(kRoot)
            End If
          End Do
          Call DDaFile(JOBOLD,2,scr,NACPAR,iDisk)
          call daxpy_(NACPAR,Scal,scr,1,D,1)
          Call DDaFile(JOBOLD,2,scr,NACPAR,iDisk)
          call daxpy_(NACPAR,Scal,scr,1,DS,1)
          Call DDaFile(JOBOLD,2,scr,NACPR2,iDisk)
          call daxpy_(NACPR2,Scal,scr,1,P,1)
          Call DDaFile(JOBOLD,2,scr,NACPR2,iDisk)
          call daxpy_(NACPR2,Scal,scr,1,PA,1)
        End Do
        Call mma_deallocate(Scr)

        If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
          Call DaClos(JOBOLD)
          JOBOLD=-1
        Else If(JOBOLD.gt.0) Then
          JOBOLD=-1
        End If


!     read from a HDF5 wavefunction file
      Else If (InVec.eq.4) then
#ifdef _HDF5_
        IF(IPRLEV.ge.TERSE) THEN
          Write(LF,'(6X,A)')
     &            'The MO-coefficients are taken from the file:'
          Write(LF,'(6X,A)') trim(mcpdft_options%wfn_file)
        END IF

        mh5id = mh5_open_file_r(mcpdft_options%wfn_file)
        call mh5_fetch_dset(mh5id, 'MO_VECTORS', CMO)
        call mh5_close_file(mh5id)
#else
        write (6,*) 'Orbitals requested from HDF5, but this'
        write (6,*) 'installation does not support that, abort!'
        call abend
#endif
      else if (invec .eq. 5) then
        write(lf,*) "FileOrb specified, but does not point to hdf5 file"
        write(lf,*) "This has not been implemented, aborting"
        call abend
      End If

! Should we ever need to do this????
!     print start orbitals
      IF(IPRLEV >= DEBUG) THEN
        ! This lene can actually be removed since it is not needed..
        ! Also, it will override the orbital energies in the runfile.
        CALL mma_allocate(ene,nTot,Label='Ene')
        CALL DCOPY_(nTot,[0.0D0],0,ene,1)
        CALL PRIMO_RASSCF_m('Input orbitals',ene,OCC,CMO)
        CALL mma_deallocate(ene)
      END IF

      END Subroutine ReadVC_m
