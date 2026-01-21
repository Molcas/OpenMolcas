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
*               2018, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine ReadVC(CMO,OCC,D,DS,P,PA,scheme)
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
*     FI      : array of real*8                                        *
*               inactive Fock matrix                                   *
*     D1I     : array of real*8                                        *
*               inactive one body density matrix                       *
*     D1A     : array of real*8                                        *
*               active one body density matrix                         *
*     TUVX    : array of real*8                                        *
*               two-electron integrals (tu!vx)                         *
*     IFINAL  : integer                                                *
*               termination flag                                       *
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
      use stdalloc, only : mma_allocate, mma_deallocate

      use rasscf_global, only : lRoots, nRoots,
     &  iRoot, Weight,
     &  nAcPar, iXsym, iAlphaBeta,
     &  iOverwr, iSUPSM, iCIrst, iPhName, nAcpr2, nOrbT,
     &  purify, iAdr15
      use general_data, only : nSym,
     &  nDel, nBas, nOrb,
     &  nTot, nTot2, Invec, LuStartOrb, StartOrbFile, JobOld,
     &  JobIph, nSSH
      use casvb_global, only: ifvb

      use orthonormalization, only : t_ON_scheme, ON_scheme_values,
     &  orthonormalize

#ifdef _HDF5_
      use mh5, only: mh5_open_file_r, mh5_exists_dset, mh5_fetch_dset,
     &               mh5_close_file
#endif
!     See comment below why this is commented out.
!     use sxci, only: IDXCI, IDXSX
      use general_data, only: CleanMask
      use PrintLevel, only: DEBUG,TERSE,VERBOSE
      use output_ras, only: LF,IPRGLB,IPRLOC
      use Definitions, only: RtoI
      use rasdim, only: MxOrb, LenIn8, MxTit, MaxBfn, MxRoot, MxSym

      implicit none

*     global data declarations
      Character(LEN=16), Parameter :: ROUTINE='READVC  '
#include "warnings.h"

      real*8 :: CMO(*),OCC(*),D(*),DS(*),P(*),PA(*)
      type(t_ON_scheme), intent(in) :: scheme

      logical :: found, changed
      integer :: iPrlev, nData, i, j, NNwOrd, iSym, iErr, IAD19, iJOB,
     &           lll, iDisk, jRoot, kRoot,
     &           iDummy(1), IADR19(30), iAD15, nTmp(8)
      real*8 :: Dummy(1), Scal
      real*8, allocatable :: CMO_copy(:)
#ifdef _HDF5_
      integer mh5id
      character(Len=maxbfn) typestring
#endif
      character(len=LENIN8*mxOrb) :: lJobH1
      character(len=2*72) :: lJobH2
      character(len=72) :: JobTit(mxTit)
      character(len=80) :: VecTit
      character(len=4) :: Label
      Integer, Allocatable:: TIND(:), NewOrd(:), TmpXSym(:), JobH(:)
      Real*8, Allocatable:: Scr(:), Ene(:), JobR(:)

      interface
        integer function isfreeunit(seed)
          integer, intent(in) :: seed
        end function
      end interface

*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
C Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*----------------------------------------------------------------------*
* Do we use default orbitals?                                          *
*----------------------------------------------------------------------*
      If(InVec.eq.0) Then
         Call qpg_darray('RASSCF orbitals',Found,nData)
         If(Found) Then
            InVec=6
            IF(IPRLEV.ge.TERSE) THEN
            Write(6,'(6x,a)') 'Orbitals from runfile: rasscf orbitals'
            END IF
         End If
      End If
      Call Check_InVec(InVec)
      If(InVec.eq.0) Then
         Call qpg_darray('SCF orbitals',Found,nData)
         If(Found) Then
            InVec=7
            IF(IPRLEV.ge.TERSE) THEN
            Write(6,'(6x,a)') 'Orbitals from runfile: scf orbitals'
            END IF
         End If
      End If
      Call Check_InVec(InVec)
      If(InVec.eq.0) Then
         Call qpg_darray('Guessorb',Found,nData)
         If(Found) Then
            InVec=5
            IF(IPRLEV.ge.TERSE) THEN
            Write(6,'(6x,a)') 'Orbitals from runfile: guessorb orbitals'
            END IF
         End If
      End If
      Call Check_InVec(InVec)
      If(Invec.eq.0) Then
         InVec=1
      End If
*----------------------------------------------------------------------*
* read from unit formatted ascii file with starting orbitals

* Note: Inside RDVEC, the file StartOrbFile is opened, but uses blindly
* the unit number provided here. So that should better be a usable
* number, or else!
      LUStartOrb=19
      LUStartOrb=IsFreeUnit(LUStartOrb)
      if(ifvb.eq.2)invec=3
      If ( InVec.eq.2 ) then
       Label='CO  '
       If (iAlphaBeta.eq.1) Label(3:3)='A'
       If (iAlphaBeta.eq.-1) Label(3:3)='B'
       if(iOverwr.eq.1) then
        CALL RDVEC(StartOrbFile,LUStartOrb,Label,NSYM,NBAS,NBAS,
     &             CMO, OCC, Dummy, iDummy, VECTIT, 0, iErr)
       else
        Label(4:4)="I"
        Call mma_allocate(TIND,maxbfn,Label='TIND')
        CALL RDVEC(StartOrbFile,LUStartOrb,Label,NSYM,NBAS,NBAS,
     &          CMO, OCC, Dummy,TIND, VECTIT, 0, iErr)
* If the typeindex array is used to resort orbitals, then if also
* a supersymmetry array is used, it has to be changed.
* The supersymmtry array is IXSYM().
* VECSORT is a utility that does not know about supersymmetry.
* So changing any orbital indices in IXSYM (or potentially any
* other orbital indices -- what about ALTER??) must be done HERE
* immediately togather with the VecSort, but not *inside* VecSort.

* But VecSort does not return any indexing information -- how are
* we to know how to change IXSYM?
* VecSort changed to include a reindexing array!
        NNwOrd=0
        Do ISym=1,NSym
         NNwOrd=NNwOrd+NBas(ISym)
        End Do
        Call mma_allocate(NEWORD,NNwOrd,Label='NewOrd')
*        Call VecSort(NSYM,NBAS,NBAS,CMO,OCC,TIND,iErr)
        Call VecSort(NSYM,NBAS,NBAS,CMO,OCC,TIND,NNwOrd,NewOrd,iErr)
* If there is a supersymmetry array, use the orbital mapping:
      If (iSUPSM.ne.0) Then
       Call mma_allocate(TMPXSYM,NNwOrd,Label='TmpXSym')
       Do I=1,NNwOrd
        J=NewOrd(I)
        TmpXSym(I)=IXSYM(J)
       End Do
       Call ICopy(NNwOrd,TmpXSym,1,IXSYM,1)
       Call mma_deallocate(TMPXSYM)
      End If

        Call mma_deallocate(NewOrd)
        Call mma_deallocate(TIND)
       endif
       Close(LUStartOrb)
       if(iErr.eq.1) then
       Write(LF,*) 'RASSCF tried to read input orbitals from a'
       Write(LF,*) 'file, but encountered an error in the'
       Write(LF,*) 'TypeIndex data.'
       call Abend()
       return
       endif

       IF(IPRLEV.ge.TERSE) THEN
        Write(LF,'(6X,A)')
     &         'The MO-coefficients are taken from the file:'
        Write(LF,'(6X,A)') trim(StartOrbFile)
        Write(LF,'(6X,A,A)') 'Title:', VecTit(2:80)
       END IF

*     read from unit JOBOLD (binary file)

      Else If ( InVec.eq.3 ) then
        IAD19=0
        iJOB=0
        Call f_Inquire('JOBOLD',Found)
        If (Found) iJOB=1
        If (iJOB.eq.1) Then
           if(JOBOLD.le.0) Then
             JOBOLD=20
             Call DaName(JOBOLD,'JOBOLD')
           end if
        Else
           If (IPRLEV.ge.TERSE) then
              Write(LF,'(6X,A)') 'File JOBOLD not found -- use JOBIPH.'
           End If
           If (JOBIPH.gt.0) Then
              JOBOLD=JOBIPH
           Else
              Call DaName(JOBOLD,IPHNAME)
           End If
        End If
        Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
        IF(IADR19(15).EQ.-1) THEN
          IAD19=0
          CALL IDAFILE(JOBOLD,2,IADR19,30,IAD19)
        ELSE
          DO I=16,30
            IADR19(I)=0
          END DO
          IF(IPRGLB.GE.VERBOSE)
     &               Call WarningMessage(1,'Old JOBIP file layout.')
        END IF
        lll = 10+RtoI
        lll = MAX(lll,mxSym)
        lll = MAX(lll,mxOrb)
        lll = MAX(lll,mxRoot)
        CALL mma_allocate(JobH,lll,Label='JobH')
        CALL mma_allocate(JobR,MxRoot,Label='JobH')
        iAd19=iAdr19(1)
        CALL WR_RASSCF_Info(JobOld,2,iAd19,
     &                      JobH(1),JobH(2),JobH(3),
     &                      JobH(4),JobH,JobH,
     &                      JobH,JobH,JobH,
     &                      mxSym,
     &                      lJobH1,LENIN8*mxOrb,JobH(5),
     &                      lJobH2,2*72,JobTit,72*mxTit,
     &                      JobR(1),JobH(6),
     &                      JobH(7),JobH,mxRoot,
     &                      JobH,JobH,JobH,
     &                      JobH(8),JobH(9),
     &                      JobH(10),
     &                      JobR(:))
        IF(IPRLEV.ge.TERSE) THEN
         If (iJOB.eq.1) Then
            Write(LF,'(6X,A)')
     &      'The MO-coefficients are taken from the file:'
            Write(LF,'(6X,A)') 'JOBOLD'
         Else
            Write(LF,'(6X,A)')
     &      'The MO-coefficients are taken from the file:'
            Write(LF,'(6X,A)') trim(iPhName)
         End If
         Write(VecTit(1:72),'(A72)') JobTit(1)
         Write(LF,'(6X,2A)') 'Title:',VecTit(1:72)
        END IF
        CALL mma_deallocate(JobH)
        CALL mma_deallocate(JobR)
        iAd19=iAdr19(2)
        Call DDaFile(JobOld,2,CMO,NTOT2,iAd19)
        Call DDaFile(JobOld,2,OCC,nTot,iAd19)
        If ( ICIRST.eq.1) then
         If ( IPRLEV.ge.VERBOSE) then
           If (iJOB.eq.1) Then
              Write(LF,'(6X,A)')
     &        'The active density matrices (D,DS,P,PA) are read from'//
     &        ' file JOBOLD and weighted together.'
           Else
              Write(LF,'(6X,A)')
     &        'The active density matrices (D,DS,P,PA) are read from'//
     &        ' file '//trim(iPhName)//' and weighted together.'
           End If
         End If
         Call mma_allocate(Scr,NACPR2,Label='Scr')
         iDisk = IADR19(3)
         Do jRoot = 1,lRoots
           Scal = 0.0d0
           Do kRoot = 1,nRoots
             If ( iRoot(kRoot).eq.jRoot ) then
               Scal = Weight(kRoot)
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
        End If
CSVC: read the L2ACT and LEVEL arrays from the jobiph file
!IFG: disabled, since it breaks when changing active space specification
        !IAD19=IADR19(18)
        !IF (IAD19.NE.0) THEN
        !  CALL IDAFILE(JOBOLD,2,IDXSX,mxAct,IAD19)
        !  CALL IDAFILE(JOBOLD,2,IDXCI,mxAct,IAD19)
        !END IF
        If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
          Call DaClos(JOBOLD)
          JOBOLD=-1
        Else If(JOBOLD.gt.0) Then
          JOBOLD=-1
        End If

*     read from a HDF5 wavefunction file
      Else If (InVec.eq.4) then
#ifdef _HDF5_
        IF(IPRLEV.ge.TERSE) THEN
          Write(LF,'(6X,A)')
     &            'The MO-coefficients are taken from the file:'
          Write(LF,'(6X,A)') trim(StartOrbFile)
        END IF

        mh5id = mh5_open_file_r(StartOrbFile)
        typestring=''
        Select Case (iAlphaBeta)
          Case (1)
            Label='CA  '
            VecTit='MO_ALPHA_TYPEINDICES'
          Case (-1)
            Label='CB  '
            VecTit='MO_BETA_TYPEINDICES'
          Case default
            Label='C   '
            VecTit='MO_TYPEINDICES'
        End Select
        Call RdVec_HDF5(mh5id,Label,NSYM,NBAS,CMO,Dummy,Dummy,iDummy)
        If (mh5_exists_dset(mh5id,Trim(VecTit)))
     &    Call mh5_fetch_dset(mh5id,Trim(VecTit),typestring)
        Call mh5_close_file(mh5id)
* Reorder orbitals based on typeindex
        If (typestring.ne.'') Then
          NNwOrd=0
          Do iSym=1,nSym
            NNwOrd=NNwOrd+NBas(iSym)
          End Do
          Call mma_allocate(TInd,NNWOrd,Label='TIND')
          Call mma_allocate(NewOrd,NNwOrd,Label='NewOrd')
          Call tpstr2tpidx(typestring,TInd,NNWOrd)
          Call VecSort(NSYM,NBAS,NBAS,CMO,OCC,TInd,
     &                                NNwOrd,NewOrd,iErr)
* If there is a supersymmetry array, use the orbital mapping:
          If (iSUPSM.ne.0) Then
            Call mma_allocate(TmpXSym,NNwOrd,Label='TmpXSym')
            Do i=1,NNwOrd
              j=NewOrd(i)
              TmpXSym(i)=iXSym(j)
            End Do
            Call iCopy(NNwOrd,TmpXSym,1,iXSym,1)
            Call mma_deallocate(TmpXSym)
          End If
          Call mma_deallocate(NewOrd)
          Call mma_deallocate(TInd)
        End If
#else
        write (6,*) 'Orbitals requested from HDF5, but this'
        write (6,*) 'installation does not support that, abort!'
        call abend()
#endif

*     guess MO-coefficients

      Else If (InVec.eq.5) then
         IF(IPRLEV.ge.VERBOSE) Write(LF,'(6x,a)')
     &                               'Detected guessorb orbitals'
         Call Qpg_dArray('Guessorb',Found,nData)
         Call Get_dArray('Guessorb',CMO,nData)
         Call Qpg_iArray('nDel_go',Found,nData)
         If(Found) Then
            Call Get_iArray('nDel_go',nTmp,nData)
            Changed=.false.
            Do iSym=1,nSym
               If(nTmp(iSym).gt.nDel(iSym)) Changed=.true.
            End Do
            If(Changed) Then
               Write(6,'(5x,a,8i5)')'Number of deleted orbitals '//
     &                            'changed from',(nDel(i),i=1,nSym)
               Write(6,'(5x,a,8i5)')'                           '//
     &                            'changed to  ',(nTmp(i),i=1,nSym)
            End If
            Do iSym=1,nSym
               If(nTmp(iSym).gt.nDel(iSym)) Then
                  nSsh(iSym)=nSsh(iSym)-nTmp(iSym)+nDel(iSym)
                  nOrb(iSym)=nOrb(iSym)-nTmp(iSym)+nDel(iSym)
                  nOrbT=nOrbT-nTmp(iSym)+nDel(iSym)
                  nDel(iSym)=nTmp(iSym)
               End If
            End Do
            IF(IPRLEV.ge.TERSE) THEN
             Write(LF,'(6X,A)')
     &       'The MO-coefficients are taken from guessorb on runfile'
            END IF
         End If
      Else If (InVec.eq.6) then
         IF(IPRLEV.ge.VERBOSE) Write(LF,'(6x,a)')
     &                               'Detected old RASSCF orbitals'
         Call qpg_darray('RASSCF orbitals',Found,nData)
         Call get_darray('RASSCF orbitals',CMO,nData)
         Call Qpg_iArray('nDel',Found,nData)
         If(Found) Then
            Call Get_iArray('nDel',nTmp,nData)
            Changed=.false.
            Do iSym=1,nSym
               If(nTmp(iSym).gt.nDel(iSym)) Changed=.true.
            End Do
            If(Changed) Then
               Write(6,'(5x,a,8i5)')'Number of deleted orbitals '//
     &                            'changed from',(nDel(i),i=1,nSym)
               Write(6,'(5x,a,8i5)')'                           '//
     &                            'changed to  ',(nTmp(i),i=1,nSym)
            End If
            Do iSym=1,nSym
               If(nTmp(iSym).gt.nDel(iSym)) Then
                  nSsh(iSym)=nSsh(iSym)-nTmp(iSym)+nDel(iSym)
                  nOrb(iSym)=nOrb(iSym)-nTmp(iSym)+nDel(iSym)
                  nOrbT=nOrbT-nTmp(iSym)+nDel(iSym)
                  nDel(iSym)=nTmp(iSym)
               End If
            End Do
            IF(IPRLEV.ge.TERSE) THEN
             Write(LF,'(6X,A,A)') 'The MO-coefficients are taken from',
     &                                    ' rasscf orbitals on runfile'
            END IF
         End If
      Else If (InVec.eq.7) then
         IF(IPRLEV.ge.VERBOSE) Write(LF,'(6x,a)')
     &                               'Detected SCF orbitals'
         Call qpg_darray('SCF orbitals',Found,nData)
         Call get_darray('SCF orbitals',CMO,nData)
         Call Qpg_iArray('nDel',Found,nData)
         If(Found) Then
            Call Get_iArray('nDel',nTmp,nData)
            Changed=.false.
            Do iSym=1,nSym
               If(nTmp(iSym).gt.nDel(iSym)) Changed=.true.
            End Do
            If(Changed) Then
               Write(6,'(5x,a,8i5)')'Number of deleted orbitals '//
     &                            'changed from',(nDel(i),i=1,nSym)
               Write(6,'(5x,a,8i5)')'                           '//
     &                            'changed to  ',(nTmp(i),i=1,nSym)
            End If
            Do iSym=1,nSym
               If(nTmp(iSym).gt.nDel(iSym)) Then
                  nSsh(iSym)=nSsh(iSym)-nTmp(iSym)+nDel(iSym)
                  nOrb(iSym)=nOrb(iSym)-nTmp(iSym)+nDel(iSym)
                  nOrbT=nOrbT-nTmp(iSym)+nDel(iSym)
                  nDel(iSym)=nTmp(iSym)
               End If
            End Do
         End If
         IF(IPRLEV.ge.TERSE) THEN
           Write(LF,'(6X,A,A)') 'The MO-coefficients are taken from',
     &                                 ' scf orbitals on runfile'
         END IF
      Else If (InVec.eq.1) then
        IF(IPRLEV.ge.VERBOSE) Write(LF,'(6X,2A)')
     &  'The MO-coefficients are obtained by diagonalizing ',
     &  'the core Hamiltonian'
        Call Guess(CMO)
      Else
       Write(LF,*) 'Severe internal bug prevents further calculation.'
       Write(LF,*) 'Invalid value for INVEC in READVC. Program stops.'
       Write(LF,*) 'Please issue bug report. INVEC=',INVEC
       CALL QUIT(_RC_GENERAL_ERROR_)
      End If
*     print start orbitals
      IF(IPRLEV.GE.DEBUG) THEN
        CALL mma_allocate(ENE,nTot,Label='ENE')
        CALL DCOPY_(nTot,[0.0D0],0,ENE,1)
        CALL PRIMO_RASSCF('Input orbitals',ENE,OCC,CMO)
        CALL mma_deallocate(ENE)
      END IF

*     cleaning orbitals for high symmetry cases

      If(Allocated(CleanMask)) Call ClnMO(CMO)
      If(PURIFY(1:6).eq.'LINEAR') CALL LINPUR(CMO)
      If(PURIFY(1:4).eq.'ATOM') CALL SPHPUR(CMO)

      if (scheme%val /= ON_scheme_values%no_ON) then
        call mma_allocate(CMO_copy, nTot2)
        CMO_copy(:nTot2) = CMO(:nTot2)
        call orthonormalize(CMO_copy, scheme, CMO(:nTot2))
        call mma_deallocate(CMO_copy)
      end if

*     save start orbitals

      IAD15=IADR15(2)
      CALL DDAFILE(JOBIPH,1,CMO,NTOT2,IAD15)
      CALL DDAFILE(JOBIPH,1,OCC,NTOT,IAD15)

*     exit

      RETURN
      END
