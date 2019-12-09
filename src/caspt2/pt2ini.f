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
      SUBROUTINE PT2INI
      USE INPUTDATA
      USE REFWFN
      USE PT2WFN
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "compiler_features.h"

      INTEGER LuSpool
C     Cholesky
      Integer iSym, iRC

      character(128) StartFile
      COMMON /datafiles/StartFile

      CALL QENTER('PT2INI')
*
* Probe the RunFile for some basic information
*
      Call Get_cArray('Seward Title',Header,144)
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iScalar('Unique atoms',nUniqAt)
      nbast=sum(nbas(1:nsym))
      nbsqt=sum(nbas(1:nsym)**2)
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nbast)
      Call DecideOnCholesky(IfChol)
* PAM Feb 2008: The following statement was moved here from
* prpctl, in case prpctl call is bypassed by keyword ''NOPROP''.
      Call Put_cArray('Relax Method','CASPT2  ',8)
*
*     Allocate Input struct, read and process the Input
*
#ifdef ALLOC_SCAL
      ALLOCATE(Input)
#endif
      LuSpool=21
      Call SpoolInp(LuSpool)
      CALL READIN_CASPT2(LuSpool,nSym)
      Call Close_LuSpool(LuSpool)
* Initialize scratch files.
      CALL OPNFLS_CASPT2
* Initialize the reference wavefunction file and read basic info
      Call refwfn_init(Input%FILE)
      Call refwfn_info
* the input processing needs some data from the reference wavefunction
* to be set, hence this phase occurs after reading that data.
      Call ProcInp_Caspt2
*
* Allocate some arrays that will stay globally allocated:
*

* Finally read the MO and CI data from the refwfn file, and close it as
* we have no more need for it and the same filename might be reused for
* the pt2wfn file. We do this after input processing because we need to
* know which roots to pick up. The MOs are stored on LUONEM, at adress
* IAD1M(1), and the CI arrays on LUCIEX at IDCIEX.
      Call refwfn_data
      Call refwfn_close

* If necessary, make modifications to the inactive/virtual orbitals. The
* new MOs, if any, will overwrite the originals on LUONEM at IAD1M(1).
      If (Input%modify_correlating_MOs) Then
        Call correlating_orbitals
      End If

* After possible reconfiguration of inactives/virtuals, the computation
* of total sizes of orbital spaces is done here, before any other code
* (that might rely on these to be correctly set!)
      Call wfnsizes

* Create the PT2 wavefunction file (formerly JOBMIX). The reference file
* should not be active, as it might be the same file (in which case it
* is overwritten).
      call pt2wfn_init
      call pt2wfn_data
*
* Global allocations active throughout the program (NOTRI was computed
* during the call to wfnsizes, so this is done here.)
*
* The total fock matrix (sum of inactive and active contrib.)
      NFIFA=NOTRI
      CALL GETMEM('LFIFA','ALLO','REAL',LFIFA,NFIFA)
* The one-electron Hamiltonian
      NHONE=NOTRI
      CALL GETMEM('LHONE','ALLO','REAL',LHONE,NHONE)
* The fock matrix with contributions from inactive orbitals, only.
      NFIMO=NOTRI
      CALL GETMEM('LFIMO','ALLO','REAL',LFIMO,NFIMO)
* The fock matrix with contributions from active orbitals, only.
      NFAMO=NOTRI
      CALL GETMEM('LFAMO','ALLO','REAL',LFAMO,NFAMO)
* Density matrices, active indices.
      CALL GETMEM('LDREF','ALLO','REAL',LDREF,NDREF)
      CALL GETMEM('LPREF','ALLO','REAL',LPREF,NPREF)
* All the density matrices kept in memory
      CALL GETMEM('LDMIX','ALLO','REAL',LDMIX,NDREF*NSTATE)
      CALL GETMEM('LDWGT','ALLO','REAL',LDWGT,NSTATE*NSTATE)

* Print input data
      CALL PRINP_CASPT2

C INITIALIZE SPLIT-GRAPH UGA TABLES.
      CALL POLY0

C Initialize superindex tables. Check memory needs.
      CALL SIZES

C Initialize sizes, offsets etc used in equation solver.
      CALL EQCTL1

      If (IfChol) then
* initialize Cholesky information
        Call Cho_X_init(irc,0.0d0)
        if (irc.ne.0) then
          write(6,*) 'CASPT2: Non-zero rc in Cho_X_init'
          CALL QUIT(irc)
        endif
* import_ch transfers some values from cholesky.fh to chocaspt2.fh
        call import_cho(numcho_pt2,infvec_n2_pt2,maxvec_pt2)
        Call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'Allo')
* get unit numbers for Cholesky MO vectors
        Do iSym=1,nSym
          Call Cho_Caspt2_OpenF(0,1,iSym,nIsplit(iSym))
          Call Cho_Caspt2_OpenF(0,2,iSym,nAsplit(iSym))
        End Do
* set up size info for cholesky vector transformation in tracho2
        call trachosz
      End If

* Allocate global orbital arrays:
      CALL GETMEM('LCMOPT2','ALLO','REAL',LCMOPT2,NCMO)
* Allocate global orbital transformation arrays:
      CALL GETMEM('TORB','ALLO','REAL',LTORB,NTORB)
      CALL GETMEM('TAT','ALLO','REAL',LTAT,NTAT)

      CALL QEXIT('PT2INI')
      END

      SUBROUTINE PT2CLS
      USE SUPERINDEX
      USE INPUTDATA
      USE PT2WFN
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "ofembed.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"

      Integer iSym
#ifndef ALLOC_SCAL
      Integer iGrp
#endif
C     Cholesky return code
      INTEGER irc
C     size of idsct array
      INTEGER NIDSCT

      If (IfChol) then
*---  Finalize Cholesky information
         Call Cho_X_Final(irc)
         if (irc.ne.0) then
            write(6,*) 'CASPT2: Non-zero rc in Cho_X_Final'
            CALL QUIT(irc)
         endif
*     Delete files of MO Cholesky vectors.
         Do iSym=1,nSym
            Call Cho_Caspt2_OpenF(3,1,iSym,nIsplit(iSym)) !close+delete
            Call Cho_Caspt2_OpenF(3,2,iSym,nAsplit(iSym)) !close+delete
         End Do
*     Deallocate memory used as index arrays in the setup section
         Call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'Free')
* NOT TESTED
#if 0
         If (Done_OFemb) Then
            Call Free_Work(ipFMaux)
         EndIf
#endif
         ! deallocate chovec_io arrays
         call trachosz_free
      End If

* Deallocate SGUGA tables:
      CALL PCLOSE

C     Deallocate MAGEB, etc, superindex tables:
      CALL SUPFREE
* Deallocate global array for Fock matrix, etc:
      CALL GETMEM('LFIFA','FREE','REAL',LFIFA,NFIFA)
      CALL GETMEM('LHONE','FREE','REAL',LHONE,NHONE)
      CALL GETMEM('LFIMO','FREE','REAL',LFIMO,NFIMO)
      CALL GETMEM('LFAMO','FREE','REAL',LFAMO,NFAMO)
      CALL GETMEM('LDREF','FREE','REAL',LDREF,NDREF)
      CALL GETMEM('LPREF','FREE','REAL',LPREF,NPREF)
      CALL GETMEM('LDMIX','FREE','REAL',LDMIX,NDREF*NSTATE)
      CALL GETMEM('LDWGT','FREE','REAL',LDWGT,NSTATE*NSTATE)
* Deallocate global orbital transformation arrays:
      CALL GETMEM('TORB','FREE','REAL',LTORB,NTORB)
      CALL GETMEM('TAT','FREE','REAL',LTAT,NTAT)
* Deallocate global orbital arrays:
      CALL GETMEM('LCMOPT2','FREE','REAL',LCMOPT2,NCMO)
* Deallocate global RHS disk offsets (allocated in eqctl1):
      NIDSCT=MXSCT*8*MXCASE*MXVEC
      CALL GETMEM('IDSCT','FREE','INTE',LIDSCT,NIDSCT)

      call pt2wfn_close
C     Close all files:
      CALL CLSFLS_CASPT2

C free input struct
#ifdef ALLOC_SCAL
      DEALLOCATE(Input)
#else
      IF (ALLOCATED(Input%MultGroup%State))
     &   DEALLOCATE(Input%MultGroup%State)
      IF (ALLOCATED(Input%nXMulState)) DEALLOCATE(Input%nXMulState)
      IF (ALLOCATED(Input%XMulGroup)) THEN
        DO iGrp = 1, SIZE(Input%XMulGroup)
          IF (ALLOCATED(Input%XMulGroup(iGrp)%State))
     &       DEALLOCATE(Input%XMulGroup(iGrp)%State)
          DEALLOCATE(Input%XMulGroup)
        END DO
      END IF
      IF (ALLOCATED(Input%nFro)) DEALLOCATE(Input%nFro)
      IF (ALLOCATED(Input%nDel)) DEALLOCATE(Input%nDel)
      IF (ALLOCATED(Input%NamFro)) DEALLOCATE(Input%NamFro)
      IF (ALLOCATED(Input%HEff)) DEALLOCATE(Input%HEff)
#endif
      End
