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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine CStart_CI_Util(C,h0,TUVX,iSel,ExplE,ExplV,iFinal)
************************************************************************
*                                                                      *
*     Find initial CI-vectors                                          *
*                                                                      *
*     calling arguments:                                               *
*     C       : array of real*8                                        *
*               CI vector                                              *
*     h0      : array of real*8                                        *
*               one-electron integrals                                 *
*     TUVX    : array of real*8                                        *
*               two-electron integrals                                 *
*     iSel    : array of integer                                       *
*               list of configuration included in the explicit Hamilt. *
*     ExplE   : array of real*8                                        *
*               eigenvalues of the explicit Hamiltonian                *
*     ExplV   : array of real*8                                        *
*               eigenvectors of the explicit Hamiltonian               *
*     iFinal  : integer                                                *
*               status of optimization process                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension C(*) , h0(*) , TUVX(*), iSel(*), ExplE(*), ExplV(*)

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "csfbas.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      PARAMETER (ROUTINE='CSTART  ')
#include "SysDef.fh"
#ifdef _HDF5_
#  include "mh5.fh"
#endif

c      Dimension iToc(15)
      Character*80 String
      Logical Exist

      Call qEnter('CStart')
      IPRLEV=IPRLOC(3)

* special case: nConf=1

      If ( nConf.eq.1 .and. nAc.eq.0 ) then
        iDisk = IADR15(4)
C$$$        ExplE(1) = C(1)  ! Commented out by Jesper
C$$$        ExplV(1) = 1.0d0 ! Commented out by Jesper
        C(1) = 1.0d0
        Call Save_tmp_CI_vec(1,nConf,C,LuDavid)
        Call qExit('CStart')
        Return
      End If

* compute the explicit Hamiltonian
* (needs to be out here as we want to redetermine
*  dynamically the selection vector)

      Call ExplH2(C,h0,TUVX,iSel,ExplE,ExplV)

* special case: nSel=nConf

      If ( nSel.eq.nConf ) then
        If ( IPRLEV.ge. DEBUG ) Write (6,*)
     &                  ' Initial CI-vectors are obtained by',
     &                  ' diagonalizing the explicit Hamiltonian'
        iDisk = IADR15(4)
        Do i = 1,lRoots
          Call dCopy_(nConf,0.0D0,0,C,1)
          Do j = 1,nSel
            k = iSel(j)
            C(k) = ExplV(j+(i-1)*nSel)
          End Do
          Call Save_tmp_CI_vec(i,nConf,C,LuDavid)
          If ( IPRLEV.ge. INSANE ) then
            Write (String,'(A,I2)') 'CI vector of root',i
            Write (String,'(A,I4,A)') '(max. ',nSel,' elements)'
            l = Min(nSel,nConf)
            Call dVcPrt(String,' ',ExplV(1+(i-1)*nSel),l)
          End If
        End Do
        Call qExit('CStart')
        Return
      End If

      If (Start_Vectors) then
        Start_Vectors=.False.
        If ( ICIRST.ne.0 ) then
#ifdef _HDF5_
          Call f_Inquire(StartOrbFile,Exist)
          if (Exist) then
            if (mh5_is_hdf5(StartOrbFile)) then
              If (IPRLEV.ge.TERSE) Then
                write (6,'(1x,a)')
     $                  'reading initial CI vectors from '//StartOrbFile
              End If
              mh5id = mh5_open_file_r(StartOrbFile)

              Call GetMem('Scr1','Allo','Real',iTmp1,nConf)
              call GetMem('kcnf','allo','inte',ivkcnf,nactel)
              Do i = 1,lRoots
                call mh5_fetch_dset_array_real(mh5id,'CI_VECTORS',
     $                  Work(iTmp1),[nconf,1],[0,i-1])
                Call Reord2(NAC,NACTEL,LSYM,1,
     &                  iWork(KICONF(1)),iWork(KCFTP),
     &                  Work(iTmp1),C,iwork(ivkcnf))
                Call Save_CI_vec(i,nConf,C,LuDavid)
              End Do
              call GetMem('kcnf','free','inte',ivkcnf,nactel)
              Call GetMem('Scr1','Free','Real',iTmp1,nConf)

              call mh5_close_file(mh5id)
            else
              Exist = .false.
            end if
          end if
          If (.not.Exist) Then
#endif
* get start vectors by reading JOBOLD
          iJOB=0
          Call f_Inquire('JOBOLD',Exist)
          If (Exist) iJOB=1
          If (iJOB.eq.1) Then
            If (IPRLEV.ge.TERSE) Then
              write (6,'(1x,a)')
     $                'reading initial CI vectors from JOBOLD'
            End If
            If (JOBOLD.le.0) Then
              JOBOLD=20
              Call DaName(JOBOLD,'JOBOLD')
            End If
          Else
            If (IPRLEV.ge.TERSE) Then
              write (6,'(1x,a)')
     $                'reading initial CI vectors from JOBIPH'
            End If
            JOBOLD=JOBIPH
          End If

          iDisk = 0
          Call IDafile(JOBOLD,2,iToc,15,iDisk)
          iDisk = iToc(4)
          Call GetMem('Scr1','Allo','Real',iTmp1,nConf)
          Do i = 1,lRoots
            Call DDafile(JOBOLD,2,Work(iTmp1),nConf,iDisk)
            call GetMem('kcnf','allo','inte',ivkcnf,nactel)
            Call Reord2(NAC,NACTEL,LSYM,1,
     &              iWork(KICONF(1)),iWork(KCFTP),
     &              Work(iTmp1),C,iwork(ivkcnf))
            call GetMem('kcnf','free','inte',ivkcnf,nactel)
            Call Save_CI_vec(i,nConf,C,LuDavid)
            If ( IPRLEV.ge. INSANE ) then
              Write (String,'(A,I2)') 'Start vector of root',i
              Write (String,'(A,I4,A)') '(max. ',nSel,' elements)'
              l = Min(nSel,nConf)
              Call dVcPrt(String,' ',C,l)
            End If
          End Do
          Call GetMem('Scr1','Free','Real',iTmp1,nConf)
          If (iJOB.eq.1) Then
            If(JOBOLD.gt.0.and.JOBOLD.NE.JOBIPH) Then
              Call DaClos(JOBOLD)
              JOBOLD=-1
            Else If(JOBOLD.gt.0) Then
              JOBOLD=-1
            End If
          End If
#ifdef _HDF5_
          End If
#endif

        Else
* no CI restart, get start vectors by diagonalizing the explicit Hamiltonian
          If ( IPRLEV.ge. DEBUG ) Write (6,*)
     &            ' Initial CI-vectors are obtained by',
     &            ' diagonalizing the explicit Hamiltonian'
          Do i = 1,lRoots
            Call dCopy_(nConf,0.0d0,0,C,1)
            Do j = 1,nSel
              k = iSel(j)
              C(k) = ExplV(j+(i-1)*nSel)
            End Do
            Call Save_CI_vec(i,nConf,C,LuDavid)
            If ( IPRLEV.ge. INSANE ) then
              Write (String,'(A,I2)') 'Start vector of root',i
              Write (String,'(A,I4,A)') '(max. ',nSel,' elements)'
              l = Min(nSel,nConf)
              Call dVcPrt(String,' ',C,l)
            End If
          End Do

        End If

      Else
* no external start vector needed, get start vectors by reading JOBIPH
        If ( iFinal.eq.2 ) then
          If ( IPRLEV.ge. DEBUG ) Write (6,*)
     &            ' Initial CI-vectors are identical to the',
     &            ' transformed CI-vectors of the previous ',
     &            ' RASSCF iteration'
        Else
          If ( IPRLEV.ge.DEBUG ) Write (6,*)
     &            ' Initial CI-vectors are identical to the',
     &            ' CI-vectors of the previous RASSCF iteration'
        End If
        iDisk = IADR15(4)
        Do i = 1,lRoots-hRoots
          Call DDafile(JOBIPH,2,C,nConf,iDisk)
          Call Save_CI_vec(i,nConf,C,LuDavid)
          If ( IPRLEV.gt.10 ) then
            Write (String,'(A,I2)') 'Start vector of root',i
            Write (String,'(A,I4,A)') '(max. ',nSel,' elements)'
            l = Min(nSel,nConf)
            Call dVcPrt(String,' ',C,l)
          End If
        End Do
*MGD simple guess for missing ones : explV
*dangerous if linear dependence with converged states
        Do i= lRoots-hRoots+1,lRoots
           Call dCopy_(nConf,0.0d0,0,C,1)
           Do j = 1,nSel
             k = iSel(j)
             C(k) = ExplV(j+(i-1)*nSel)
           End Do
          Call Save_CI_vec(i,nConf,C,LuDavid)
        End Do

      End If

      Call qExit('CStart')

      Return
      End
