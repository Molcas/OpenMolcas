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
Module LUCIA_INTERFACE
Private
Public Lucia_Util

Contains
!***********************************************************************
!  Lucia_Util
!
!> @brief
!>   Wrapper for using LUCIA utils in MOLCAS.
!> @author Jesper Wisborg Krogh
!>
!> @details
!> By using the LUCIA utils through this wrapper it is guaranteed
!> that all common blocks used have a common parent routine.
!>
!> @param[in] Module Identifier
!***********************************************************************
      Subroutine Lucia_Util(Module, iSym, iDisk, LU, Array, RVec)
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
      use strbas
      use rasscf_lucia
#include "implicit.fh"
      Character(LEN=*) Module
      Integer, Optional:: iSym
      Integer, Optional:: iDisk
      Integer, Optional:: LU
      Real*8, Optional:: Array(:)
      Real*8, Optional:: RVEC(:)

      Parameter(MxpLnc = 72)
      Character(LEN=MxpLnc) Module_
!
! Include all LUCIA include files to make sure
! they are available during the calculation.
!
#include "mxpdim.fh"
#include "cands.fh"

#include "cecore.fh"
#include "cgas.fh"

#include "cicisp.fh"
#include "cintfo.fh"
#include "clunit.fh"

#include "cprnt.fh"
#include "crun.fh"

#include "csm.fh"
#include "cstate.fh"



#include "gasstr.fh"
#include "intform.fh"
#include "irat.fh"

#include "lucinp.fh"

#include "oper.fh"
#include "orbinp.fh"
#include "spinfo_lucia.fh"
#include "stinf.fh"
#include "strinp.fh"
      Integer, Allocatable:: lVec(:)
!
#ifdef _DEBUGPRINT_
      INTEGER, SAVE :: COUNTER = 0
      COUNTER=COUNTER+1
      WRITE(6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)')
     & '[',COUNTER,']','ENTRY LUCIA_UTIL','(',Module,')'
#endif
!
! Make sure the Module variable is in upper case.
!
      Module_ = Module
      Call UppCas(Module_,MxpLnc)
!
! Call the appropriate routines according to Module
!
      If (Module_(1:4) .eq. 'DIAG') Then
         Call Diag_Master
      Else If (Module_(1:9) .eq. 'SIGMA_CVB') Then
!        iSym_LI is the symmetry to be used.
         Call Sigma_Master_CVB(C_POINTER,SIZE(C_POINTER),iSym)
      Else If (Module_(1:5) .eq. 'SIGMA') Then
!        write(6,*) 'blubbbbbbhc'
         Call Sigma_Master(C_POINTER,SIZE(C_POINTER))
      Else If (Module_(1:5) .eq. 'TRACI') Then
!        write(6,*) 'blubbbbbbtraci'
!        iDisk is the initial disk address (for read/write of JOBIPH)
!        Lu is the file unit for JOBIPH
!        Array is the transformation matrix (not sorted as LUCIA needs it).
         Call mma_allocate(lVec,MXNTTS,Label='lVec')
         Call Traci_Master(iDisk,LU,Array,lVec)
         Call mma_deallocate(lVec)
      Else If (Module_(1:5) .eq. 'DENSI') Then
         If (Present(RVEC)) Then
            Call Densi_Master(C_POINTER,SIZE(C_POINTER),RVEC=RVEC(:))
         Else
            Call Densi_Master(C_POINTER,SIZE(C_POINTER))
         End If
      Else If (Module_(1:3) .eq. 'INI') Then
         Call Lucia_Ini()
         Call DetCtl_Gas()
      Else If (Module_(1:5) .eq. 'CLOSE') Then
         Call DetCtl_Free()
         Call Lucia_Close()
      Else
         Write(6,*) 'Unknown module requested in Lucia_Util.'
         Write(6,*) 'Module = ',Module
         Write(6,*) 'Known modules are:'
         Write(6,*) 'Diag, Sigma, Sigma_CVB, Densi, DetCtl, Ini'
         Call Abend
      End If

#ifdef _DEBUGPRINT_
      WRITE(6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)')
     & '[',COUNTER,']','EXIT LUCIA_UTIL','(',Module,')'
#endif
      End Subroutine Lucia_Util





      SUBROUTINE densi_master(CIVec,nCIVec,RVec)
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
      use rasscf_lucia, only: kvec3_length, iSigma_on_Disk, PAtmp, Ptmp, DSTmp, Dtmp
!
! Controls the calculation of the densities, when Lucia is called
! from Molcas Rasscf.
!
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "crun.fh"
#include "cicisp.fh"
#include "clunit.fh"
#include "orbinp.fh"
#include "lucinp.fh"
#include "spinfo_lucia.fh"
#include "cstate.fh"
#include "io_util.fh"
      Integer nCIVec
      Real*8 CIVec(nCIVEC)
      Real*8, Optional:: RVec(:)

      logical iPack,tdm
      dimension dummy(1)
      Real*8, Allocatable:: VEC1(:), VEC2(:)
      Integer, Allocatable:: lVec(:)
      Real*8, Allocatable, Target:: SCR1(:), SCR3(:)
      Real*8, Allocatable:: SCR2(:), SCR4(:)
      Integer NSD, NCSF

!
! Put CI-vector from RASSCF on luc
!
!     if rvec is associated, it should be a pointer to a second CI vector
!     and a one-particle transition density matrix will be computed
      tdm = Present(rVec)

      NSD = NSD_PER_SYM(IREFSM)
      NCSF= NCSF_PER_SYM(IREFSM)
      CALL mma_allocate(SCR1,NSD,Label='SCR1')
      CALL mma_allocate(SCR2,NSD,Label='SCR2')

      CALL COPVEC(CIVEC,SCR1,NCSF)

      Call mma_allocate(lVec,MXNTTS,Label='lVec')
      IF (tdm) THEN
         CALL mma_allocate(SCR3,NSD,Label='SCR3')
         CALL mma_allocate(SCR4,NSD,Label='SCR4')

         CALL COPVEC(rvec,SCR3,NCSF)
         CALL CSDTVC(SCR3,SCR4,1,DTOC,SDREO, IREFSM, 1)
         CALL CPCIVC(SCR3,NSD,LUHC, MXNTTS, IREFSM, 1,lVec)
         CALL mma_deallocate(SCR3)
         CALL mma_deallocate(SCR4)
      END IF
      CALL CSDTVC(SCR1,SCR2,1,DTOC,SDREO, IREFSM, 1)
      CALL CPCIVC(SCR1,NSD,LUC, MXNTTS, IREFSM, 1,lVec)
      Call mma_deallocate(lVec)

!
! Determine length of arrays VEC1 and VEC2
!
!      IF(ISIMSYM.EQ.0) THEN
         LBLOCK = MXSOOB
!      ELSE
!         LBLOCK = MXSOOB_AS
!      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
!PAM06      LBLOCK = MAX(XISPSM(IREFSM,1),DBLE(MXSOOB))
      LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
      IF(PSSIGN.NE.0.0D0) LBLOCK = 2*INT(XISPSM(IREFSM,1))
!
! Allocate arrays
!
      Call mma_allocate(VEC1,LBLOCK,Label='VEC1')
      Call mma_allocate(VEC3,kvec3_length,Label='VEC3')
!
! Copy Sigma-vector from disc to core
!
      Call mma_allocate(VEC2,LBLOCK,Label='VEC2')
       IF (iSigma_on_disk .ne. 0) THEN
          Call mma_allocate(lVec,MXNTTS,Label='lVec')
          CALL cpsivc(lusc34, mxntts, vec2,lVec)
          Call mma_deallocate(lVec)
       ELSE
          vec2(:) = 0.0d0
       ENDIF
!
! Information needed on file handling
!
      LBLK = - 1
!
! Copy vector on file LUC to LUSC1 and LUHC
!
      IDISK(LUC)=0
      IDISK(LUSC1)=0
      CALL COPVCD(LUC,LUSC1,VEC1,0,LBLK)
      IF (.not.tdm) CALL COPVCD(LUSC1,LUHC,VEC1,1,LBLK)
!
! Calculate one- and two-body densities
!
      IPACK = .TRUE.
      DUMMY = 0.0D0
      IF (tdm) THEN
         CALL densi2_lucia(1,Dtmp,dummy,dummy,dummy,vec1,vec2,lusc1,luhc,exps2,1,DStmp,IPACK)
      ELSE
         CALL densi2_lucia(2,rho1,dummy,Ptmp,PAtmp,vec1,vec2,lusc1,luhc,exps2,1,srho1,IPACK)
      END IF

!
! Explanation of calling parameters
!
!      2      : DONE!!! - Calculate both one and two body densities.
!      rho1  : DONE!!! - Output - include in module glbbas
!      krho2  : DONE!!! - Output - include in moduke glbbas
!      vec1  : DONE!!! - CI-vector
!      vec2  : DONE!!! - Sigma-vector
!      lusc1  : DONE!!! - file pointer
!      luhc   : DONE!!! - file pointer
!      exps2  : DONE!!! - Output - expectation value of S**2.
!      1      : DONE!!! - Calculate spin density
!      srho1 : DONE!!! - Comming with module glbbas
!
      IF (.not.tdm) THEN
!        Save densities in trigonal format for use in Molcas
!
         CALL TriPak(rho1, Dtmp, 1, ntoob, ntoob)
         CALL TriPak(srho1, DStmp, 1, ntoob, ntoob)
      END IF
!
      CALL CSDTVC(scr1,scr2,2,dtoc,SDREO, iRefSm, 1)
!
      CALL mma_deallocate(SCR1)
      CALL mma_deallocate(SCR2)
      Call mma_deallocate(VEC1)
      Call mma_deallocate(VEC2)
      Call mma_deallocate(VEC3)
!
      END SUBROUTINE densi_master
!
End Module LUCIA_INTERFACE
