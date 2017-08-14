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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2003-2005, Valera Veryazov                             *
*               2017, Roland Lindh                                     *
************************************************************************
      subroutine SCF_DMET(ireturn)
************************************************************************
*                                                                      *
*     purpose: perform RHF calculations                                *
*                                                                      *
*     calls to: OpnFls_SCF,ReadIn_SCF,SOrb,WrInp_SCF,WfCtl_SCF,Final,  *
*               ClsFls_SCF                                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, P. Borowski, M. Schuetz, and R. Lindh              *
*     University of Lund, Sweden, 1992-1998                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use SCF_Arrays
      Implicit Real*8 (a-h,o-z)
*
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "stdalloc.fh"
#include "twoswi.fh"
#include "file.fh"
#include "hflda.fh"
#include "warnings.fh"
      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
*
      Character*8 EMILOOP
      Logical FstItr, Semi_Direct
      Real*8 SIntTh
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      write(6,*) "scf 0"
      Call CWTime(TCPU1,TWall1)
      write(6,*) "scf 1"
      Call qEnter('SCF')
      write(6,*) "scf 2"
      HFLDA=0.0
      Call SCF_Init()
      write(6,*) "scf 3"
      iTerm=0
*
      Call OpnFls_SCF
      write(6,*) "scf 4"
      Call ReadIn_SCF(SIntTh)
      write(6,*) "scf 5"
      LuOrb=LuInp
*
      Semi_Direct = DSCF .and. (nDisc.ne.0 .or. (nCore.ne.0
     &      .and. nDisc.eq.0))
      If (Semi_Direct) Then
         Call GetMem('MaxMem','Max','Real',iDum,MemSew)
         MemLow=Min(MemSew/2,1024*1024)
         MemSew=Max(MemSew/10,MemLow)
         Call xsetmem_ints(MemSew)
      End If
*
      Call Init_SCF()
      write(6,*) "scf 6"
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up the starting orbitals and set the occupation numbers.
*
      LuOrb=LuInp
      nD = iUHF + 1
      Call SOrb(LuOrb,SIntTh,iTerm)
      write(6,*) "scf 7"
      Call OccDef(OccNo,nnB,nD,CMO,nBB)
      write(6,*) "scf 8"
*
      Call mma_deallocate(HDiag)
      If (Aufb) Then
         lthH = nBB
      Else
         lthH = nOV
      End If
      Call mma_allocate(HDiag,lthH,nD,Label='HDiag')
*                                                                      *
************************************************************************
*                                                                      *
      Call WrInp_SCF(SIntTh)
      write(6,*) "scf 9"

      Call Cre_SCFWfn
      write(6,*) "scf 10"

      FstItr=.True.

      If(.not.OnlyProp) Then
         Call WfCtl_SCF_DMET(iTerm,KSDFT,FstItr,SIntTh)
            write(6,*) "scf after WfCtl_SCF_0"
      End If

*     so that iPsLst is right in case of nIter==0
            write(6,*) "scf after WfCtl_SCF_1"
      If (nIter(nIterP).eq.0) iter0=-1
            write(6,*) "scf after WfCtl_SCF_2"
      Call Final_DMET()
            write(6,*) "scf after WfCtl_SCF_3"
      If (DSCF) Call Free_TLists
            write(6,*) "scf after WfCtl_SCF_4"
*
      Call CWTime(TCPU2,TWall2)
            write(6,*) "scf after WfCtl_SCF_5"
      Call SavTim(4,TCPU2-TCPU1,TWall2-TWall1)
            write(6,*) "scf after WfCtl_SCF_6"
*
      Call GMFree_DMET()
            write(6,*) "scf after WfCtl_SCF_7"
      Call ClsFls_SCF
            write(6,*) "scf after WfCtl_SCF_8"
      If (Semi_Direct) Call xRlsMem_Ints
            write(6,*) "scf after WfCtl_SCF_9"
*
*     Call MolDen Interface
*
      If(iUHF.eq.0) Then
            write(6,*) "scf after WfCtl_SCF_10"
         Call Molden_Interface(iUHF,'SCFORB','MD_SCF',AddFragments)
c         Call grid_driver(-1,'SCF','SCFORB',iRc)
      Else
            write(6,*) "scf after WfCtl_SCF_11"
         Call Molden_Interface(iUHF,'UHFORB','MD_SCF',AddFragments)
c         Call grid_driver(-1,'SCF','UNAORB',iRc)
            write(6,*) "scf after WfCtl_SCF_12"
      End If
            write(6,*) "scf after WfCtl_SCF_13"
      Call qExit('SCF')
            write(6,*) "scf after WfCtl_SCF_14"
      if(iStatPRN.gt.0) then
            write(6,*) "scf after WfCtl_SCF_15"
       Call qStat(' ')
       Call FastIO('STATUS')
      endif
*
*.... Everything has to come to an end...
*
      iReturn=iTerm
*
      If (Do_OFemb) Then
         Call GetEnvF('EMIL_InLoop',EMILOOP)
         If (EMILOOP.eq.' ') EMILOOP='0'
         If (EMILOOP(1:1).ne.'0') Then
            If (iReturn.ne._RC_ALL_IS_WELL_) Then
               Call WarningMessage(1,'SCF: non-zero return code.')
            EndIf
            iReturn=_RC_CONTINUE_LOOP_
            Call Check_FThaw(iReturn)
         EndIf
      EndIf
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      return
*
      End
************************************************************************
