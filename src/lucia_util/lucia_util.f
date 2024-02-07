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
* Copyright (C) Jesper Wisborg Krogh                                   *
************************************************************************
*  Lucia_Util
*
*> @brief
*>   Wrapper for using LUCIA utils in MOLCAS.
*> @author Jesper Wisborg Krogh
*>
*> @details
*> By using the LUCIA utils through this wrapper it is guaranteed
*> that all common blocks used have a common parent routine.
*>
*> @param[in] Module Identifier
*> @param[in] iArg2   Argument to LUCIA
*> @param[in] Array1 Argument to LUCIA
************************************************************************
      Subroutine Lucia_Util(Module)
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS
      use strbas
      use rasscf_lucia
      use Lucia_Interface, only: iSym_LI, iDisk_LI, Lu_LI, Array_LI
#include "implicit.fh"

      Parameter(MxpLnc = 72)
      Character(LEN=*) Module
      Character(LEN=MxpLnc) Module_
*
* Include all LUCIA include files to make sure
* they are available during the calculation.
*
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
*
#ifdef _DEBUGPRINT_
      INTEGER, SAVE :: COUNTER = 0
      COUNTER=COUNTER+1
      WRITE(6,'(1X,A1,I6,A1,1X,A,1X,A1,A,A1)')
     & '[',COUNTER,']','ENTRY LUCIA_UTIL','(',Module,')'
#endif
*
* Make sure the Module variable is in upper case.
*
      Module_ = Module
      Call UppCas(Module_,MxpLnc)
*
* Call the appropriate routines according to Module
*
      If (Module_(1:4) .eq. 'DIAG') Then
         Call Diag_Master
      Else If (Module_(1:9) .eq. 'SIGMA_CVB') Then
*        iSym_LI is the symmetry to be used.
         Call Sigma_Master_CVB(iSym_LI)
      Else If (Module_(1:5) .eq. 'SIGMA') Then
!        write(6,*) 'blubbbbbbhc'
         Call Sigma_Master(C_POINTER,SIZE(C_POINTER))
      Else If (Module_(1:5) .eq. 'TRACI') Then
!        write(6,*) 'blubbbbbbtraci'
*        iDisk_LI is the initial disk address (for read/write of JOBIPH)
*        Lu_LI is the file unit for JOBIPH
*        Array_LI is the transformation matrix (not sorted as LUCIA needs it).
         Call mma_allocate(lVec,MXNTTS,Label='lVec')
         Call Traci_Master(iDisk_LI,LU_LI,Array_LI,lVec)
         Call mma_deallocate(lVec)
      Else If (Module_(1:5) .eq. 'DENSI') Then
         Call Densi_Master(C_POINTER,SIZE(C_POINTER))
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
      Return
      End
