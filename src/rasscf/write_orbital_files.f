************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      module write_orbital_files

        use stdalloc, only : mma_allocate, mma_deallocate
        use general_data, only : nSym
        use gas_data, only : iDoGas, nGAS

        implicit none
        private
        public :: OrbFiles, get_typeidx, putOrbFile,
     &      write_orb_per_iter
        save

        interface get_typeidx
          module procedure :: RAS_get_typeidx, GAS_get_typeidx
        end interface

        logical :: write_orb_per_iter = .false.

        interface
          integer function isfreeunit(iseed)
            integer, intent(in) :: iseed
          end function
        end interface
      contains

      Subroutine OrbFiles(JOBIPH, IPRLEV)
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      use fortran_strings, only : str

      use rasscf_global, only : iToc, BName, header, title, lRoots,
     & nRoots,
     &  iRoot, iPt2, Weight, iOrbTyp,
     &  FDiag, E2Act, maxorbout
      use general_data, only : nActel, iSpin, stSym,
     &  nFro, nIsh, nAsh, nDel, nBas, nRs1, nRs2, nRs3, nHole1, nElec3,
     &  nTot, nTot2, nConf
      use gas_data, only : nGssh
      use stdalloc, only: mma_allocate, mma_deallocate
      use printlevel, only: USUAL
      use output_ras, only: LF
      Implicit None

#include "rasdim.fh"
      integer, intent(in) :: JobIph, iPrlev

      integer :: iDisk, iRt, iNDType(7, 8), lUVVVec
      real*8 :: Energy, PotNucDummy

      character(len=80) :: VecTyp
      character(len=128) :: Filename
      real*8, allocatable:: CMO(:), Occ(:), Ene(:), EDum(:)
      interface
        integer function isfreeunit(iseed)
          integer, intent(in) :: iseed
        end function
      end interface

* This routine is used at normal end of a RASSCF optimization, or
* when using the OrbOnly keyword to create orbital files.
*-------------------------------------------------------------------
* NOTE:
* PAM 2008: Before this subroutine replaced RASREAD, the orbital energies
* were sent as subroutine argument when rasread was called from rasscf.
* But the real argument, in rasscf, was FDIAG, which turns out to be a
* fixed array FDIAG(mxorb) in rasscf_global.F90.
* Rather than checking why that array is there, and how the values are
* put there, I simply use the array FDIAG in rasscf_global.F90, when needing
* orbital energies to WrVec calls. This may need checking later...
*----------------------------------------------------------------------*
*     Read the JobIph file to get the required information             *
*----------------------------------------------------------------------*
* PAM Jan 2014 -- do not take POTNUC from JOBIPH; take it directly
* from runfile, where it was stored by seward.
      iDisk = 0
      Call iDaFile(JobIph,2,iToc,15,iDisk)
      iDisk = iToc(1)
      Call WR_RASSCF_Info(JobIph,2,iDisk,
     &                    nActEl,iSpin,nSym,stSym,
     &                    nFro,nIsh,nAsh,nDel,
     &                    nBas,mxSym,BName,LENIN8*mxOrb,nConf,
     &                    Header,144,Title,4*18*mxTit,PotNucDummy,
     &                    lRoots,nRoots,iRoot,mxRoot,
     &                    nRs1,nRs2,nRs3,
     &                    nHole1,nElec3,iPt2,Weight)
      ntot = sum(nBas(:nSym))
      ntot2 = sum(nBas(:nSym)**2)
*----------------------------------------------------------------------*
*     Allocate CMO array                                               *
*----------------------------------------------------------------------*
      call mma_allocate(CMO,ntot2,Label='CMO')
      call mma_allocate(Occ,ntot,Label='Occ')
*----------------------------------------------------------------------*
*     Make typeindex information                                       *
*----------------------------------------------------------------------*
      if (.not. iDoGas) then
        IndType(:, :) =
     &      get_typeidx(nFro, nIsh, nRs1, nRs2, nRs3, nBas, nDel)
      else
        IndType(:, :) = get_typeidx(nFro, nIsh, nGSSH, nBas, nDel)
      endif
*----------------------------------------------------------------------*
*     First, write orbitals to RasOrb:                                 *
* IORBTYP=1 for 'Average' orbitals... Default!
* IORBTYP=2 for 'Canonical' orbitals.
* IORBTYP=3 for 'Natural' orbitals... in this case the number of roots need to be specified
* IORBTYP=4 for 'Spin' orbitals... in this case the number of roots need to be specified
*----------------------------------------------------------------------*
      filename = 'RASORB'
      If (iOrbTyp .ne. 2) then
        iDisk=iToc(2)
        Call dDaFile(JobIph,2,CMO,ntot2,iDisk)
        IF(IPRLEV.GE.USUAL) then
          Write(LF,'(6X,3A)') 'Average orbitals are written to the ',
     &                          trim(filename),' file'
        end if
        VecTyp = '* RASSCF average (pseudo-natural) orbitals'
        Call dDaFile(JobIph,2,Occ,ntot,iDisk)
      Else
        iDisk=iToc(9)
        Call dDaFile(JobIph,2,CMO,ntot2,iDisk)
        IF(IPRLEV.GE.USUAL) then
          Write(LF,'(6X,3A)') 'Canonical orbitals are written to the ',
     &          trim(filename),' file'
        end if
        VecTyp = '* RASSCF canonical orbitals for CASPT2'
        call dcopy_(ntot,[1.0D0],0,Occ,1)
      End If
*----------------------------------------------------------------------*
*     Write  orbitals                                                  *
*----------------------------------------------------------------------*
      LuvvVec=50
      LuvvVec=isfreeunit(LuvvVec)
c      Call WrVec(filename,LuvvVec,'COE',nSym,nBas,nBas,
c     &           CMO, Occ, FDIAG, iDummy,VecTyp)
      Call WrVec_(filename,LuvvVec,'COET',0,nSym,nBas,nBas,
     &            CMO,CMO,Occ,Occ,FDIAG,[E2act],indType,VecTyp,0)
c      Call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,
c     &           CMO, Occ, FDIAG, IndType,VecTyp)
      Call WrVec_(filename,LuvvVec,'AIT',0,nSym,nBas,nBas,
     &            CMO,CMO,Occ,Occ,FDIAG,[E2act],
     &            indType,VecTyp,0)
*----------------------------------------------------------------------*
*     Second, write natural orbitals
*----------------------------------------------------------------------*
      Call mma_allocate(Ene,mxRoot*mxIter,Label='Ene')
      Call Get_dArray('Last energies',Ene,lRoots)

      iDisk=iToc(12)
      DO IRT=1, MIN(MAXORBOUT, LROOTS, 999)
        energy=Ene(IRT)
        if (irt < 999) then
          filename = 'RASORB.'//str(IRT)
        else
          filename = 'RASORB.x'
        end if
        Call dDaFile(JobIph,2,CMO,ntot2,iDisk)
        Call dDaFile(JobIph,2,Occ,ntot,iDisk)
        IF(IPRLEV.GE.USUAL) then
          Write(LF,'(6X,A,I3,3A)') 'Natural orbitals for root ', IRT,
     &    ' are written to the ', trim(filename), ' file'
        end if
        Write(VecTyp,'(A41,I3,A3,f22.12)')
     &   '* RASSCF natural orbitals for root number',IRT,
     &   ' E=',Energy
*----------------------------------------------------------------------*
*     Write  orbitals                                                  *
*----------------------------------------------------------------------*
        LuvvVec=50
        LuvvVec=isfreeunit(LuvvVec)
        Call WrVec(filename,LuvvVec,'COE',nSym,nBas,nBas,
     &             CMO, Occ, FDIAG, IndType,VecTyp)
        Call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,
     &             CMO, Occ, FDIAG, IndType,VecTyp)
      END DO
      Call mma_deallocate(Ene)
*----------------------------------------------------------------------*
*     Third, write spin density orbitals
*----------------------------------------------------------------------*
      iDisk=iToc(14)
      Call mma_allocate(EDum,NTot,Label='EDum')
      EDum(:)=0.0D0
      DO IRT=1,MIN(MAXORBOUT,LROOTS,999)
        if (irt < 999) then
          filename = 'SPDORB.'//str(IRT)
        else
          filename = 'SPDORB.x'
        end if
        Call dDaFile(JobIph,2,CMO,ntot2,iDisk)
        Call dDaFile(JobIph,2,Occ,ntot,iDisk)
        IF (IPRLEV.GE.USUAL) then
          Write(LF,'(6X,A,I3,3A)')'Spin density orbitals for root ',
     &      irt, ' are written to the ',trim(filename),' file'
        end if
        Write(VecTyp,'(A,I3)')
     &   '* RASSCF spin density orbitals for root number',IRT
*----------------------------------------------------------------------*
*     Write  orbitals                                                  *
*----------------------------------------------------------------------*
        LuvvVec = 50
        LuvvVec = isfreeunit(LuvvVec)
        Call WrVec(filename,LuvvVec,'CEO',nSym,nBas,nBas,
     &             CMO, Occ, EDum, IndType,VecTyp)
        Call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,
     &             CMO, Occ, EDum, IndType,VecTyp)
      END DO
      Call mma_deallocate(EDum)
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
      call mma_deallocate(CMO)
      call mma_deallocate(Occ)

      Return
      End subroutine


      function RAS_get_typeidx(
     &    nFro, nIsh, nRs1, nRs2, nRs3, nBas, nDel) result(typeidx)
        integer, intent(in) ::  nFro(:), nIsh(:), nRs1(:), nRs2(:),
     &    nRs3(:), nBas(:), nDel(:)
        integer :: typeidx(7, 8)


        typeidx(1, :nSym) = nFro(:nSym)
        typeidx(2, :nSym) = nIsh(:nSym)
        typeidx(3, :nSym) = nRS1(:nSym)
        typeidx(4, :nSym) = nRS2(:nSym)
        typeidx(5, :nSym) = nRS3(:nSym)
        typeidx(7, :nSym) = nDel(:nSym)

        typeidx(6, :nSym) = 0
        typeidx(6, :nSym) = nBas(:nSym) - sum(typeidx(:, :nSym), dim=1)
      end function RAS_get_typeidx

      function GAS_get_typeidx(
     &      nFro, nIsh, nGSSH, nBas, nDel) result(typeidx)
        integer, intent(in) ::
     &    nFro(:), nIsh(:), nBas(:), nGSSH(:, :), nDel(:)
        integer :: typeidx(7, 8)

        typeidx(1, :nSym) = nFro(:nSym)
        typeidx(2, :nSym) = nIsh(:nSym)
        typeidx(3, :nSym) = 0
        typeidx(4, :nSym) = sum(nGssh(1:nGAS, :nSym), dim=1)
        typeidx(5, :nSym) = 0
        typeidx(7, :nSym) = nDel(:nSym)

        typeidx(6, :nSym) = 0
        typeidx(6, :nSym) = nBas(:nSym) - sum(typeidx(:, :nSym), dim=1)
      end function GAS_get_typeidx

      subroutine putOrbFile(CMO, orbital_E, iDoGAS)
        use general_data, only : ntot,
     &    nFro, nIsh, nRs1, nRs2, nRs3, nDel, nBas
        use gas_data, only : nGSSH
        real*8, intent(in) :: CMO(:), orbital_E(:)
        logical, intent(in) :: iDoGAS

        character(len=*), parameter :: filename = 'ORTHORB'
        real*8, allocatable :: occ_number(:)
        integer, parameter :: arbitrary_magic_number = 50
        integer :: file_id, typeidx(7, 8)
        character(len=80) ::
     &    orbfile_title = 'Orbitals after Orthonormalization.'

        file_id = arbitrary_magic_number
        file_id = isfreeunit(file_id)
        if (.not. iDoGas) then
          typeidx = get_typeidx(nFro, nIsh, nRs1, nRs2, nRs3, nBas,nDel)

        else
          typeidx = get_typeidx(nFro, nIsh, nGSSH, nBas, nDel)
        endif

! TODO(Oskar): Implement proper occupation number reading.
        call mma_allocate(occ_number, nTot)
        occ_number(:) = 1.d0
        call WrVec(filename, file_id, 'COIE', nSym, nBas, nBas,
     &             CMO, occ_number, orbital_E, typeidx, orbfile_title)
        call mma_deallocate(occ_number)
      end subroutine putOrbFile

      end module write_orbital_files
