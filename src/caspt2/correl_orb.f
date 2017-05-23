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
      subroutine correlating_orbitals
************************************************************************
*
* Alter the correlating orbital space by various freeze-delete schemes.
*
************************************************************************
      use inputdata
      implicit none
#include "stdalloc.fh"
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"

      Real*8, Allocatable :: CMO(:), DPQ(:)
      Integer IDISK
      Integer ntri, NDPQ
      Real*8 Dummy(1)
      Integer I, iRC, iSkp, iSym

* memory to store MOs
      NCMO=NBSQT
      CALL MMA_ALLOCATE(CMO,NCMO)

* Read the MOs from the LUONEM file
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,CMO,NCMO,IDISK)

* AFreeze CASPT2 calculation, available only with
* Cholesky or RI type integral representation.
      If(Input % AFreeze) then
        If (.not.IfChol) Then
          Call WarningMessage(2,'AFreeze needs Cholesky/RI.')
          Call Quit_OnUserError
        End If
        Write(6,'(A)') ' Additional orbitals will be frozen or deleted'
        Write(6,'(A,18A4)') ' Selected atoms:  ',
     &    (Input%namfro(i),i=1,Input%lnfro)
        Write(6,'(A,8I4)')
     &  ' Frozen orbitals before selection:    ',(nfro(i),i=1,nsym)
        Write(6,'(A,8I4)')
     &  ' Inactive orbitals before selection:  ',(nish(i),i=1,nsym)
        Write(6,'(A,8I4)')
     &  ' Secondary orbitals before selection: ',(nssh(i),i=1,nsym)
        Write(6,'(A,8I4)')
     &  ' Deleted orbitals before selection:   ',(ndel(i),i=1,nsym)
        ntri=0
        Do isym=1,nsym
          ntri=ntri+(nbas(isym)+nbas(isym)**2)/2
        Enddo
        NDPQ=ntri
        CALL MMA_ALLOCATE(DPQ,NDPQ)
        Call AFreez(NSYM,NBAS,NFRO,NISH,NASH,NSSH,NDEL,NAME,
     &    INPUT%NAMFRO,INPUT%LNFRO,DPQ,
     &    Input%THRFR,Input%THRDE,IFQCAN,CMO,NCMO)
        CALL MMA_DEALLOCATE(DPQ)
        Write(6,'(A,8I4)')
     &  ' Frozen orbitals after selection     ',(nfro(i),i=1,nsym)
        Write(6,'(A,8I4)')
     &  ' Inactive orbitals after selection   ',(nish(i),i=1,nsym)
        Write(6,'(A,8I4)')
     &  ' Secondary orbitals after selection: ',(nssh(i),i=1,nsym)
        Write(6,'(A,8I4)')
     &  ' Deleted orbitals after selection:   ',(ndel(i),i=1,nsym)
      Endif

* LOV CASPT2 calculation, available only with
* Cholesky or RI type integral representation.
      If (Input % LovCASPT2) then
        If (.not.IfChol) Then
          Call WarningMessage(2,'LOV-CASPT2 needs Cholesky/RI.')
          Call Quit_OnUserError
        End If
        IF (IFQCAN.EQ.0) Then
          Call WarningMessage(2,'LOV-CASPT2 needs Canonical Orbitals.')
          Call Quit_OnUserError
        EndIf
        If (Input%thr_atm.lt.0.0d0 .or. Input%thr_atm.ge.1.0d0) Then
          write(6,*)' Threshold out of range! Must be in [0,1[ '
          Call Quit_OnUserError
        End If

        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A)') ' Start LovCASPT2 section '
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A,8I4)')
        Write(6,'(A,8I4)')
     &  ' Frozen orbitals before selection:   ',(nFro(i),i=1,nSym)
        Write(6,'(A,8I4)')
     &  ' Inactive orbitals before selection: ',(nIsh(i),i=1,nSym)
        Write(6,'(A,8I4)')
     &  ' Secondary orbitals before selection:',(nSsh(i),i=1,nSym)
        Write(6,'(A,8I4)')
     &  ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

        EMP2=0.0d0
        Call Lov_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &    NAME,nUniqAt,Input%thr_atm,IFQCAN,
     &    Input%DoMP2,Input%DoEnv,Input%VIRA,EMP2,CMO,NCMO)

        If (irc.ne.0) Then
          write(6,*) 'LovCASPT2 returned rc= ',irc
          Call abend()
        EndIf
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A)') ' End LovCASPT2 section '
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A,8I4)')
        iSkp=0
        Do iSym=1,nSym
          iSkp=Max(iSkp,nAsh(iSym))
        End Do
        If (iSkp.lt.1) Call xquit(0)
        Write(6,'(A,8I4)') ' Going to perform CASPT2 calculation '//
     &                      'on the active region only.'
        Write(6,'(A,8I4)')
      Endif

* Frozen Natural Orbital CASPT2 calculation, available only with
* Cholesky or RI type integral representation.
      If (Input % FnoCASPT2) then
        If (.not.IfChol) Then
          Call WarningMessage(2,'FNO-CASPT2 needs Cholesky/RI.')
          Call Quit_OnUserError
        End If
        If (Input%vFrac.le.0.0d0 .or. Input%vFrac.gt.1.0d0) Then
          Call WarningMessage(2,'FNO-CASPT2 fraction out of range.')
          Write(6,*)' Requested fraction of virtual space must be'
          Write(6,*)' between 0.0 and 1.0.'
          Call Quit_OnUserError
        End If

        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A)') ' Start FNO-CASPT2 section '
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A,8I4)')
        Write(6,'(A,I3,A)') ' NOs specified: ',
     &        int(Input%vfrac*1.0D2),'% of the total virtual space'
        Write(6,'(A,8I4)')
     &  ' Secondary orbitals before selection:',(nSsh(i),i=1,nSym)
        Write(6,'(A,8I4)')
     &  ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

        EMP2=0.0d0
        Call FNO_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &           Input%vfrac,IFQCAN,Input%DoMP2,EMP2,CMO,NCMO)

        If (irc.ne.0) Then
          write(6,*) 'FNO_CASPT2 returned rc= ',irc
          Call abend()
        EndIf
        Write(6,'(A,8I4)')
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A)') ' End FNO-CASPT2 section '
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A,8I4)')
        Write(6,'(A,8I4)')
      Endif

      If (Input % GhostDelete) then
        If (Input%ThrGD.lt.0.0d0 .or. Input%ThrGD.ge.1.0d0) Then
          Write(6,*)' GHOST threshold out of range! Must be in [0,1[ '
          Call Quit_OnUserError
        End If

        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A)') ' GHOST virtual space removal'
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A,8I4)')
        Write(6,'(A,8I4)')
     &  ' Secondary orbitals before selection:',(nSsh(i),i=1,nSym)
        Write(6,'(A,8I4)')
     &  ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

        Call Delete_Ghosts(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,
     &          NAME,nUniqAt,Input%ThrGD,.True.,CMO,Dummy)

        If (irc.ne.0) Then
          write(6,*) 'Delete_GHOSTS returned rc= ',irc
          Call abend()
        EndIf
        Write(6,'(A,8I4)')
        Write(6,'(A)')
     &  '-------------------------------------------------------'
        Write(6,'(A,8I4)')
        Write(6,'(A,8I4)')
      Endif

* Store the MOs on the LUONEM file
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,1,CMO,NCMO,IDISK)

      CALL MMA_DEALLOCATE(CMO)

* we need to force recanonicalization of the orbitals later
      IFQCAN=0

      return
      end
