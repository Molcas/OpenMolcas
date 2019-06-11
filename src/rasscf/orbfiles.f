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
      Subroutine OrbFiles(JOBIPH, IPRLEV)
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      use fortran_strings, only : str

      use rasscf_data, only : iToc, name, header, title, lRoots, nRoots,
     &  iRoot, LENIN8, mXORB, mxTit, mXroot, iPt2, Weight, iOrbTyp,
     &  FDiag, E2Act, mxiter, maxorbout
      use general_data, only : nActel, iSpin, nSym, lSym, mXSym,
     &  nFro, nIsh, nAsh, nDel, nBas, nRs1, nRs2, nRs3, nHole1, nElec3,
     &  nTot, nTot2, nConf
      use gugx_data, only : ifCas
      use gas_data, only : nGssh, iDoGas, nGas

      implicit none
#include "output_ras.fh"
      Parameter (ROUTINE='ORBFILES')
#include "WrkSpc.fh"
      integer, intent(in) :: JobIph, iPrlev

      integer :: iDisk, iRt, iSym, lCMO, iShift, iNDT, iNDType(56),
     &    lUVVVec, isfreeunit, lEdum, nGast, ipEne,ipOcc
      real*8 :: Energy, PotNucDummy

      character(len=80) :: VecTyp
      character(len=128) :: Filename
#ifndef _DMRG_
      logical :: doDMRG = .false.
#endif

      call qEnter(routine)
* This routine is used at normal end of a RASSCF optimization, or
* when using the OrbOnly keyword to create orbital files.
*-------------------------------------------------------------------
* NOTE:
* PAM 2008: Before this subroutine replaced RASREAD, the orbital energies
* were sent as subroutine argument when rasread was called from rasscf.
* But the real argument, in rasscf, was FDIAG, which turns out to be a
* fixed array FDIAG(mxorb) in rasscf.fh.
* Rather than checking why that array is there, and how the values are
* put there, I simply use the array FDIAG in rasscf.fh, when needing
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
     &                    nActEl,iSpin,nSym,lSym,
     &                    nFro,nIsh,nAsh,nDel,
     &                    nBas,mxSym,Name,LENIN8*mxOrb,nConf,
     &                    Header,144,Title,4*18*mxTit,PotNucDummy,
     &                    lRoots,nRoots,iRoot,mxRoot,
     &                    nRs1,nRs2,nRs3,
     &                    nHole1,nElec3,iPt2,Weight)
      ntot = sum(nBas(:nSym))
      ntot2 = sum(nBas(:nSym)**2)
*----------------------------------------------------------------------*
*     Allocate CMO array                                               *
*----------------------------------------------------------------------*
      call getmem('CMO','allo','real',LCMO,ntot2)
      call getmem('Occ','allo','real',ipOcc,ntot)
*----------------------------------------------------------------------*
*     Make typeindex information                                       *
*----------------------------------------------------------------------*
      iShift=0
      DO ISYM=1,NSYM
        IndT=0
        IndType(1+iShift)= NFRO(ISYM)
        IndT=IndT+NFRO(ISYM)
        IndType(2+iShift)= NISH(ISYM)
        IndT=IndT+NISH(ISYM)
        If (.not. iDoGas) Then
          IndType(3+iShift)= NRS1(ISYM)
          IndT=IndT+NRS1(ISYM)
          IndType(4+iShift)= NRS2(ISYM)
          IndT=IndT+NRS2(ISYM)
          IndType(5+iShift)= NRS3(ISYM)
          IndT=IndT+NRS3(ISYM)
        Else
          IndType(3+iShift)=0
          NGAST=SUM(NGSSH(1:NGAS,ISYM))
          IndType(4+iShift)=ngast
          IndT=IndT+ngast
          IndType(5+iShift)=0
        End If
        IndType(7+iShift)= NDEL(ISYM)
        IndT=IndT+NDEL(ISYM)
        IndType(6+iShift)= NBAS(ISYM)-IndT
        iShift=iShift+7
      EndDo
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
        Call dDaFile(JobIph,2,Work(lCMO),ntot2,iDisk)
        IF(IPRLEV.GE.USUAL) then
          Write(LF,'(6X,3A)') 'Average orbitals are written to the ',
     &                          trim(filename),' file'
        end if
        VecTyp = '* RASSCF average (pseudo-natural) orbitals'
        Call dDaFile(JobIph,2,Work(ipOcc),ntot,iDisk)
      Else
        iDisk=iToc(9)
        Call dDaFile(JobIph,2,Work(lCMO),ntot2,iDisk)
        IF(IPRLEV.GE.USUAL) then
          Write(LF,'(6X,3A)') 'Canonical orbitals are written to the ',
     &          trim(filename),' file'
        end if
        VecTyp = '* RASSCF canonical orbitals for CASPT2'
        call dcopy_(ntot,[1.0D0],0,Work(ipOcc),1)
      End If
*----------------------------------------------------------------------*
*     Write  orbitals                                                  *
*----------------------------------------------------------------------*
      LuvvVec = 50
      LuvvVec = isfreeunit(LuvvVec)
      call my_WrVec_('COET')
      call my_WrVec_('AIT')

*----------------------------------------------------------------------*
*     Second, write natural orbitals
*----------------------------------------------------------------------*
      Call GetMem('Ene','Allo','Real',ipEne,mxRoot*mxIter)
      Call Get_dArray('Last energies',Work(ipEne),lRoots)

      iDisk=iToc(12)
      DO IRT=1, MIN(MAXORBOUT, LROOTS, 999)
        energy=Work(ipEne+IRT-1)

        if(doDMRG)then
#ifdef _DMRG_
          energy = dmrg_energy%dmrg_state_specific(irt)
#endif
        end if
        filename = 'RASORB.'//merge(str(IRT), 'x', irt < 999)
        Call dDaFile(JobIph,2,Work(lCMO),ntot2,iDisk)
        Call dDaFile(JobIph,2,Work(ipOcc),ntot,iDisk)
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
     &    Work(lCMO), Work(ipOcc), FDIAG, IndType,VecTyp)
        Call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,
     &   Work(lCMO), Work(ipOcc), FDIAG, IndType,VecTyp)
      END DO
      Call GetMem('Ene','Free','Real',ipEne,mxRoot*mxIter)
*----------------------------------------------------------------------*
*     Third, write spin density orbitals
*----------------------------------------------------------------------*
      iDisk=iToc(14)
      Call GetMem('EDummy','Allo','Real',LEDum,NTot)
      call dcopy_(NTot,[0.0D0],0,Work(LEDum),1)
      DO IRT=1,MIN(MAXORBOUT,LROOTS,999)
        filename = 'SPDORB.'//merge(str(IRT), 'x', irt < 999)
        Call dDaFile(JobIph,2,Work(lCMO),ntot2,iDisk)
        Call dDaFile(JobIph,2,Work(ipOcc),ntot,iDisk)
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
     &    Work(lCMO), Work(ipOcc), Work(LEDum), IndType,VecTyp)
        Call WrVec(filename,LuvvVec,'AI',NSYM,NBAS,NBAS,
     &   Work(lCMO), Work(ipOcc), Work(LEDum), IndType,VecTyp)
      END DO
      Call GetMem('EDummy','Free','Real',LEDUM,NTOT)
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
      call getmem('CMO','free','real',LCMO,ntot2)
      call getmem('Occ','free','real',ipOcc,ntot)

      Call qExit(routine)
      Return
      contains
        subroutine partial_WrVec_(label)
          character(*), intent(in) :: label
          Call WrVec_(filename,LuvvVec,label,0,nSym,nBas,nBas,
     &            Work(lCMO),Work(lCMO),
     &            Work(ipOcc),Work(ipOcc),
     &            FDIAG,E2act,
     &            indType,VecTyp,0)
        end subroutine
      End subroutine
