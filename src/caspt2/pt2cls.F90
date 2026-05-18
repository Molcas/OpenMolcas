!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE PT2CLS()
      USE INPUTDATA, ONLY: CLEANUP_INPUT
      USE SUPERINDEX, ONLY: SUPFREE
      use PT2WFN, ONLY: PT2WFN_CLOSE
      use sguga, only: SGS, CIS, EXS, SG_Free
      use caspt2_global, only: FIMO, FIFA, DREF, PREF, DMIX,            &
     &                         DWGT, CMOPT2, TAT, TORB, IDSCT, Weight,  &
     &                         IDCIEX, IDTCEX
      use stdalloc, only: mma_deallocate
#ifdef _DMRG_
      use qcmaquis_interface, only:qcmaquis_interface_deinit
      use qcmaquis_interface_cfg, only:dmrg_file
      use qcmaquis_info, only: qcmaquis_info_deinit
      use caspt2_module, only: DMRG
#endif
! NOT TESTED
#if 0
      use OFembed, only: FMaux
#endif
      use ChoCASPT2, only: NASplit,NISplit,NumCho_PT2
      use caspt2_module, only: IfChol, nAsh, nIsh, nSsh, nSym
      use definitions, only: iwp, u6
      IMPLICIT NONE

      Integer(kind=iwp) iSym
!     Cholesky return code
      INTEGER(kind=iwp) irc
!     size of idsct array

      If (IfChol) then
!---  Finalize Cholesky information
         Call Cho_X_Final(irc)
         if (irc.ne.0) then
            Write(u6,*) 'CASPT2: Non-zero rc in Cho_X_Final'
            CALL QUIT(irc)
         endif
!     Delete files of MO Cholesky vectors.
         Do iSym=1,nSym
            Call Cho_Caspt2_OpenF(3,1,iSym,nIsplit(iSym)) !close+delete
            Call Cho_Caspt2_OpenF(3,2,iSym,nAsplit(iSym)) !close+delete
         End Do
!     Deallocate memory used as index arrays in the setup section
         Call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'Free')
! NOT TESTED
#if 0
         Call mma_deallocate(FMaux,safe='*')
#endif
         ! deallocate chovec_io arrays
         call trachosz_free()
      End If

! Deallocate SGUGA tables:
      CALL SG_Free(SGS,CIS,EXS)

! dealloacte DMRG stuff
#ifdef _DMRG_
      if (DMRG) then
        call mma_deallocate(dmrg_file%qcmaquis_checkpoint_file)
        call qcmaquis_info_deinit()
        call qcmaquis_interface_deinit()
      end if
#endif


!     Deallocate MAGEB, etc, superindex tables:
      CALL SUPFREE()
! Deallocate global array for Fock matrix, etc:
      CALL mma_deallocate(IDCIEX)
      CALL mma_deallocate(IDTCEX)
      CALL mma_deallocate(FIFA)
      CALL mma_deallocate(FIMO)
      CALL mma_deallocate(DREF)
      CALL mma_deallocate(PREF)
      CALL mma_deallocate(DMIX)
      CALL mma_deallocate(DWGT)
! Deallocate global orbital transformation arrays:
      CALL mma_deallocate(TORB)
      CALL mma_deallocate(TAT)
! Deallocate global orbital arrays:
      CALL mma_deallocate(CMOPT2)
! Deallocate global RHS disk offsets (allocated in eqctl1):
      CALL mma_deallocate(IDSCT)

      call pt2wfn_close()
!     Close all files:
      CALL CLSFLS_CASPT2()

! Weight array is allocated in refwfn_info
      call mma_deallocate(Weight)

! free input struct
      CALL CleanUp_Input()
      End SUBROUTINE PT2CLS
