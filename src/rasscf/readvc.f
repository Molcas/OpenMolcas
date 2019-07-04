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
      Subroutine ReadVC(CMO,OCC,D,DS,P,PA)
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

      use rasscf_data, only : lRoots, nRoots,
     &  iRoot, LENIN8, mxTit, Weight, mXOrb, mXroot,
     &  maxorbout, nAcPar, iXsym, iAlphaBeta,
     &  iOverwr, iSUPSM, iCIrst, iPhName, nAcpr2, nOrbT, iClean,
     &  purify, iAdr15
      use general_data, only : nSym, mXSym,
     &  nDel, nBas, nOrb,
     &  nTot, nTot2, Invec, LuStartOrb, StartOrbFile, JobOld,
     &  JobIph, nSSH, maxbfn, mXAct

      implicit none

*     global data declarations
#include "output_ras.fh"
      Parameter (ROUTINE='READVC  ')
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "warnings.fh"
#include "wadr.fh"
#include "casvb.fh"
#include "raswfn.fh"
      Common /IDSXCI/ IDXCI(mxAct),IDXSX(mxAct)

      real*8, intent(in) :: CMO(*),OCC(*),D(*),DS(*),P(*),PA(*)

      logical :: found, changed
      integer :: iPrlev, nData, ifvb,
     &    i, j, iTIND, NNwOrd, iSym,
     &    LNEWORD, LTMPXSYM, iErr, IAD19, iJOB,
     &    lll, lJobH, ldJobH, lscr, iDisk,
     &    jRoot, kRoot, iDXsX, idXCI,
     &    iDummy(1), IADR19(30), iAD15, lEne, nTmp(8)
      real*8, allocatable :: CMOO(:)
      real*8 :: Dummy(1), Scal
#ifdef _HDF5_
      integer mh5id
      character(Len=maxbfn) typestring
#endif
      character(LENIN8*mxOrb) :: lJobH1
      character(2*72) :: lJobH2
      character(72) :: JobTit(mxTit)
      character(80) :: VecTit
      character(4) :: Label

      interface
        integer function isfreeunit(seed)
          integer, intent(in) :: seed
        end function

        integer function ip_of_Work_i(A)
          integer :: A
        end function
      end interface

*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call qEnter('ReadVc')
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
        Call GetMem('TIND','Allo','Inte',iTIND,maxbfn)
        CALL RDVEC(StartOrbFile,LUStartOrb,Label,NSYM,NBAS,NBAS,
     &          CMO, OCC, Dummy,iWork(iTIND), VECTIT, 0, iErr)
* If the typeindex array is used to resort orbitals, then if also
* a supersymmetry array is used, it has to be changed.
* The supersymmtry array is IXSYM().
* VECSORT is a utility that does not know about supersymmetry.
* So changing any orbital indices in IXSYM (or potentially any
* other orbital indices -- what about ALTER??) must be done HERE
* immediately togather with the VecSort, but not *inside* VecSort.

* But VecSort does not return any indicing information -- how are
* we to know how to change IXSYM?
* VecSort changed to include a reindexing array!
        NNwOrd=0
        Do ISym=1,NSym
         NNwOrd=NNwOrd+NBas(ISym)
        End Do
        Call GetMem('NewOrd','Allo','Inte',LNEWORD,NNwOrd)
*        Call VecSort(NSYM,NBAS,NBAS,CMO,OCC,iWork(iTIND),iErr)
        Call VecSort(NSYM,NBAS,NBAS,CMO,OCC,iWork(iTIND),
     &                              NNwOrd,iWork(lNewOrd),iErr)
* If there is a supersymmetry array, use the orbital mapping:
      If (iSUPSM.ne.0) Then
       Call GetMem('TmpXSym','Allo','Inte',LTMPXSYM,NNwOrd)
       Do I=1,NNwOrd
        J=iWork(lNewOrd-1+I)
        iWork(lTmpXSym-1+I)=IXSYM(J)
       End Do
       Call ICopy(NNwOrd,iWork(lTmpXSym),1,IXSYM,1)
       Call GetMem('TmpXSym','Free','Inte',LTMPXSYM,NNwOrd)
      End If

        Call GetMem('NewOrd','Free','Inte',LNEWORD,NNwOrd)
        Call GetMem('TIND','Free','Inte',iTIND,maxbfn)
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
        Write(LF,'(6X,A)') StartOrbFile
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
              Write(LF,*) '  File JOBOLD not found -- use JOBIPH.'
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
        lll = 1
        lll = MAX(lll,mxSym)
        lll = MAX(lll,mxOrb)
        lll = MAX(lll,RtoI)
        lll = MAX(lll,RtoI*mxRoot)
        CALL GETMEM('JOBOLD','ALLO','INTEGER',lJobH,lll)
        ldJobH=ip_of_Work_i(iWork(lJobH))
        iAd19=iAdr19(1)
        CALL WR_RASSCF_Info(JobOld,2,iAd19,
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      mxSym,
     &                      lJobH1,LENIN8*mxOrb,iWork(lJobH),
     &                      lJobH2,2*72,JobTit,72*mxTit,
     &                      Work(ldJobH),iWork(lJobH),
     &                      iWork(lJobH),iWork(lJobH),mxRoot,
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      Work(ldJobH))
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
        CALL GETMEM('JOBOLD','FREE','INTEGER',lJobH,lll)
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
         Call GetMem('Scr','Allo','Real',lscr,NACPR2)
         iDisk = IADR19(3)
         Do jRoot = 1,lRoots
           Scal = 0.0d0
           Do kRoot = 1,nRoots
             If ( iRoot(kRoot).eq.jRoot ) then
               Scal = Weight(kRoot)
             End If
           End Do
           Call DDaFile(JOBOLD,2,Work(lscr),NACPAR,iDisk)
           call daxpy_(NACPAR,Scal,Work(lscr),1,D,1)
           Call DDaFile(JOBOLD,2,Work(lscr),NACPAR,iDisk)
           call daxpy_(NACPAR,Scal,Work(lscr),1,DS,1)
           Call DDaFile(JOBOLD,2,Work(lscr),NACPR2,iDisk)
           call daxpy_(NACPR2,Scal,Work(lscr),1,P,1)
           Call DDaFile(JOBOLD,2,Work(lscr),NACPR2,iDisk)
           call daxpy_(NACPR2,Scal,Work(lscr),1,PA,1)
         End Do
         Call GetMem('Scr','Free','Real',lscr,NACPR2)
        End If
CSVC: read the L2ACT and LEVEL arrays from the jobiph file
         IAD19=IADR19(18)
         IF (IAD19.NE.0) THEN
           CALL IDAFILE(JOBOLD,2,IDXSX,mxAct,IAD19)
           CALL IDAFILE(JOBOLD,2,IDXCI,mxAct,IAD19)
         END IF
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
          Write(LF,'(6X,A)') StartOrbFile
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
          Call GetMem('TInd','Allo','Inte',iTInd,NNWOrd)
          Call GetMem('NewOrd','Allo','Inte',LNewOrd,NNwOrd)
          Call tpstr2tpidx(typestring,iWork(iTInd),NNWOrd)
          Call VecSort(NSYM,NBAS,NBAS,CMO,OCC,iWork(iTInd),
     &                                NNwOrd,iWork(lNewOrd),iErr)
* If there is a supersymmetry array, use the orbital mapping:
          If (iSUPSM.ne.0) Then
            Call GetMem('TmpXSym','Allo','Inte',LTmpXSym,NNwOrd)
            Do i=1,NNwOrd
              j=iWork(lNewOrd-1+i)
              iWork(lTmpXSym-1+i)=iXSym(j)
            End Do
            Call iCopy(NNwOrd,iWork(lTmpXSym),1,iXSym,1)
            Call GetMem('TmpXSym','Free','Inte',LTmpXSym,NNwOrd)
          End If
          Call GetMem('NewOrd','Free','Inte',LNewOrd,NNwOrd)
          Call GetMem('TInd','Free','Inte',iTInd,maxbfn)
        End If
#else
        write (6,*) 'Orbitals requested from HDF5, but this'
        write (6,*) 'installation does not support that, abort!'
        call abend
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
        CALL GETMEM('DumE','Allo','Real',LENE,nTot)
        CALL DCOPY_(nTot,[0.0D0],0,WORK(LENE),1)
        CALL PRIMO_RASSCF('Input orbitals',WORK(LENE),OCC,CMO)
        CALL GETMEM('DumE','Free','Real',LENE,nTot)
      END IF

*     cleaning orbitals for high symmetry cases

      If(iClean.ne.0) Call ClnMO(CMO)
      If(PURIFY(1:6).eq.'LINEAR') CALL LINPUR(CMO)
      If(PURIFY(1:4).eq.'ATOM') CALL SPHPUR(CMO)

*     orthogonalize the molecular orbitals
      call mma_allocate(CMOO, nTot2)
      CMOO(:nTot2) = CMO(:nTot2)
      CALL ONCMO(CMOO, CMO)
      call mma_deallocate(CMOO)

*     save start orbitals

      IAD15=IADR15(2)
      CALL DDAFILE(JOBIPH,1,CMO,NTOT2,IAD15)
      CALL DDAFILE(JOBIPH,1,OCC,NTOT,IAD15)

*     exit

      CALL QEXIT('READVC')
      RETURN
      END
