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

#ifdef _HDF5_
      use mh5, only: mh5_open_file_r, mh5_fetch_dset, mh5_close_file
#endif
      use sxci_pdft, only: idxci, idxsx
      use mcpdft_output, only: terse, verbose, debug, lf, iPrGlb, iPrLoc

      Implicit Real*8 (A-H,O-Z)

*     global data declarations

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "warnings.h"
#include "wadr.fh"
*     calling arguments

      Dimension CMO(*),OCC(*),D(*),DS(*),P(*),PA(*)

*     local data declarations

      Character*72 JobTit(mxTit)
      DIMENSION IADR19(30)
      Character*80 VecTit
      Logical Found
      Logical Changed
      Integer nTmp(8)
      Character*(LENIN8*mxOrb) lJobH1
      Character*(2*72) lJobH2
#ifdef _HDF5_
      integer mh5id
#endif

*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
C Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV >= DEBUG) THEN
        WRITE(LF,*)' Entering READVC'
      END IF
*----------------------------------------------------------------------*
* Do we use default orbitals?                                          *
*----------------------------------------------------------------------*
      If(InVec == 0) Then
        Call qpg_darray('RASSCF orbitals',Found,nData)
        If(Found) Then
          InVec=6
          IF(IPRLEV >= TERSE) THEN
            Write(6,'(6x,a)') 'Orbitals from runfile: rasscf orbitals'
          END IF
        End If
      End If
      Call Check_InVec_m(InVec)
      If(InVec == 0) Then
         Call qpg_darray('SCF orbitals',Found,nData)
         If(Found) Then
            InVec=7
            IF(IPRLEV >= TERSE) THEN
            Write(6,'(6x,a)') 'Orbitals from runfile: scf orbitals'
            END IF
         End If
      End If
      Call Check_InVec_m(InVec)
      If(InVec == 0) Then
         Call qpg_darray('Guessorb',Found,nData)
         If(Found) Then
            InVec=5
            IF(IPRLEV >= TERSE) THEN
            Write(6,'(6x,a)') 'Orbitals from runfile: guessorb orbitals'
            END IF
         End If
      End If
      Call Check_InVec_m(InVec)
      If(Invec == 0) Then
        Write(LF,'(6X,2A)')
     &  "MC-PDFT shouldn't be guessing orbitals...something wrong",
     &  'with your calculation or this module..'
        call abend()
      End If
*----------------------------------------------------------------------*
* read from unit formatted ascii file with starting orbitals

* Note: Inside RDVEC, the file StartOrbFile is opened, but uses blindly
* the unit number provided here. So that should better be a usable
* number, or else!

*     read from unit JOBOLD (binary file)
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
              Call DaName(JOBOLD,IPHNAME)
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
        lll = 10+RtoI
        lll = MAX(lll,mxSym)
        lll = MAX(lll,mxOrb)
        lll = MAX(lll,RtoI*mxRoot)
        CALL GETMEM('JOBOLD','ALLO','INTEGER',lJobH,lll)
        ldJobH=ip_of_Work_i(iWork(lJobH+10))
        iAd19=iAdr19(1)
        CALL WR_RASSCF_Info(JobOld,2,iAd19,
     &                      iWork(lJobH),iWork(lJobH+1),iWork(lJobH+2),
     &                      iWork(lJobH+3),iWork(lJobH),iWork(lJobH),
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      mxSym,
     &                      lJobH1,LENIN8*mxOrb,iWork(lJobH+4),
     &                      lJobH2,2*72,JobTit,72*mxTit,
     &                      Work(ldJobH),iWork(lJobH+5),
     &                      iWork(lJobH+6),iWork(lJobH),mxRoot,
     &                      iWork(lJobH),iWork(lJobH),iWork(lJobH),
     &                      iWork(lJobH+7),iWork(lJobH+8),
     &                      iWork(lJobH+9),
     &                      Work(ldJobH))
        IF(IPRLEV >= TERSE) THEN
         If (iJOB == 1) Then
            Write(LF,'(6X,A)')
     &      'The MO-coefficients are taken from the file:'
            Write(LF,'(6X,A)') 'JOBOLD'
         Else
            Write(LF,'(6X,A)')
     &      'The MO-coefficients are taken from the file:'
            Write(LF,'(6X,A)') trim(IPHNAME)
         End If
         Write(VecTit(1:72),'(A72)') JobTit(1)
         Write(LF,'(6X,2A)') 'Title:',VecTit(1:72)
        END IF
        CALL GETMEM('JOBOLD','FREE','INTEGER',lJobH,lll)
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
     &       ' file '//trim(IPHNAME)//
     &       ' and weighted together.'
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
        call mh5_fetch_dset(mh5id, 'MO_VECTORS', CMO)
        call mh5_close_file(mh5id)
#else
        write (6,*) 'Orbitals requested from HDF5, but this'
        write (6,*) 'installation does not support that, abort!'
        call abend
#endif

*     guess MO-coefficients

      Else If (InVec == 5) then
         IF(IPRLEV >= VERBOSE) Write(LF,'(6x,a)')
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
            IF(IPRLEV >= TERSE) THEN
             Write(LF,'(6X,A)')
     &       'The MO-coefficients are taken from guessorb on runfile'
            END IF
         End If
      Else If (InVec == 6) then
         IF(IPRLEV >= VERBOSE) Write(LF,'(6x,a)')
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

      End If
*     print start orbitals
      IF(IPRLEV >= DEBUG) THEN
        CALL GETMEM('DumE','Allo','Real',LENE,nTot)
        CALL DCOPY_(nTot,[0.0D0],0,WORK(LENE),1)
        CALL PRIMO_RASSCF_m('Input orbitals',WORK(LENE),OCC,CMO)
        CALL GETMEM('DumE','Free','Real',LENE,nTot)
      END IF

*     orthogonalize the molecular orbitals
* New orthonormalization routine, with additional deletion of
* linear dependence.
      CALL GETMEM('CMOO','ALLO','REAL',LCMOO,NTOT2)
      CALL DCOPY_(NTOT2,CMO,1,WORK(LCMOO),1)
      CALL GETMEM('CMOO','FREE','REAL',LCMOO,NTOT2)

      RETURN
      END
