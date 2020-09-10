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
      subroutine SCF(ireturn)
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
#include "interfaces_scf.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call CWTime(TCPU1,TWall1)
      Call qEnter('SCF')
      HFLDA=0.0
      Call SCF_Init()
      iTerm=0
*
      Call OpnFls_SCF
      Call ReadIn_SCF(SIntTh)
      LuOrb=LuInp
*
      Semi_Direct = DSCF .and. (nDisc.ne.0 .or. (nCore.ne.0
     &      .and. nDisc.eq.0))
      If (Semi_Direct) Then
         Call mma_MaxDBLE(MemSew)
         MemLow=Min(MemSew/2,1024*1024)
         MemSew=Max(MemSew/10,MemLow)
         Call xsetmem_ints(MemSew)
      End If
*
      Call Init_SCF()
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up the starting orbitals and set the occupation numbers.
*
      LuOrb=LuInp
      nD = iUHF + 1
      Call SOrb(LuOrb,SIntTh,iTerm)
      Call OccDef(OccNo,nnB,nD,CMO,nBB)
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

      Call Cre_SCFWfn

      FstItr=.True.

      If(.not.OnlyProp) Then
         Call WfCtl_SCF(iTerm,KSDFT,FstItr,SIntTh)
      End If

*     so that iPsLst is right in case of nIter==0
      If (nIter(nIterP).eq.0) iter0=-1
      Call Final()
      If (DSCF) Call Free_TLists
*
      Call CWTime(TCPU2,TWall2)
      Call SavTim(4,TCPU2-TCPU1,TWall2-TWall1)
*
      Call GMFree()
      Call ClsFls_SCF
      If (Semi_Direct) Call xRlsMem_Ints
*
*     Call MolDen Interface
*
      If(iUHF.eq.0) Then
         Call Molden_Interface(iUHF,'SCFORB','MD_SCF')
      Else
         Call Molden_Interface(iUHF,'UHFORB','MD_SCF')
      End If
      Call qExit('SCF')
      if(iStatPRN.gt.0) then
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
      SubRoutine IniLLs
*     initialize the diverse linked lists
      Implicit Real*8 (a-h,o-z)

#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "llists.fh"
#include "lnklst.fh"
*
#ifdef _DEBUG_
      Call QEnter('IniLLs')
#endif
*
*     MemRsv set tentatively to the size of six density matrices
c     MemRsv=6*nBT
      LLlist=0
      LLGrad=0
      Call IniLst(LLGrad,20)
      Call IniLst(LLDgrd,20)
      Call IniLst(LLDelt,20)
      Call IniLst(LLy,20)
      Call IniLst(LLx,MxOptm)
      Init_LLs=1
*
#ifdef _DEBUG_
      Call QExit('IniLLs')
#endif
      Return
      End
*----------------------------------------------------------------------*
#ifdef _NOTUSED_
      Subroutine StatLLS()
      Implicit Real*8 (a-h,o-z)
#include "llists.fh"
      If (Init_LLs.eq.1) Then
         Call StlLst(LLGrad)
         Call StlLst(LLDgrd)
         Call StlLst(LLDelt)
         Call StlLst(LLy)
         Call StlLst(LLx)
      Else
         Write (6,*) '****** W A R N I N G ! ******'
         Write (6,*) ' Linked lists are not there!'
      End If
      Return
      End
#endif
*----------------------------------------------------------------------*
      SubRoutine KiLLs
*     dispose the diverse linked lists
      Implicit Real*8 (a-h,o-z)
#include "llists.fh"
      If (Init_LLs.eq.1) Then
         Call KilLst(LLGrad)
         Call KilLst(LLDgrd)
         Call KilLst(LLDelt)
         Call KilLst(LLy)
         Call KilLst(LLx)
         Init_LLs=-1
      Else
         Write (6,*) '****** W A R N I N G ! ******'
         Write (6,*) ' Linked list already killed!'
      End If
      Return
      End
*----------------------------------------------------------------------*
      Subroutine RclLLs(iDskPt)
      Implicit Real*8 (a-h,o-z)
#include "infso.fh"
#include "file.fh"
#include "llists.fh"
      Integer iDskPt(5)
      Call RclLst(LLGrad,LuGrd,iDskPt(1),MemRsv)
      Call RclLst(LLDgrd,LuDGd,iDskPt(2),MemRsv)
      Call RclLst(LLDelt,LuDel,iDskPt(3),MemRsv)
      Call RclLst(LLy   ,Lux  ,iDskPt(4),MemRsv)
      Call RclLst(LLx   ,Luy  ,iDskPt(5),MemRsv)
      Init_LLs=1
*     Call StatLLs
      Return
      End
*----------------------------------------------------------------------*
      Subroutine DmpLLs(iDskPt)
      Implicit Real*8 (a-h,o-z)
#include "file.fh"
#include "llists.fh"
      Integer iDskPt(5)
      If (Init_LLs.eq.1) Then
*        Call StatLLs
         Call DmpLst(LLGrad,LuGrd,iDskPt(1))
         Call DmpLst(LLDgrd,LuDGd,iDskPt(2))
         Call DmpLst(LLDelt,LuDel,iDskPt(3))
         Call DmpLst(LLy   ,Lux  ,iDskPt(4))
         Call DmpLst(LLx   ,Luy  ,iDskPt(5))
      Else
         Write (6,*) '****** W A R N I N G ! ******'
         Write (6,*) ' Linked list already killed!'
      End If
      Return
      End
*----------------------------------------------------------------------*
      Subroutine StlLst(LLink)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      return
       Write (6,*)
       Write (6,*) '*********** Status of Linked List *************'
       Write (6,*)
       Write (6,*) ' LLink:',LLink
       Write (6,*)
       Write (6,*) ' CNOD data'
       Write (6,*) 'Error code:                       ',iWork(LLink  )
       Write (6,*) 'Pointer to first NODE in the list:',iWork(LLink+1)
       Write (6,*) 'Actual length of list:            ',iWork(LLink+2)
       Write (6,*) '# of vectors in core:             ',iWork(LLink+3)
       Write (6,*)
       iRoot=iWork(LLink+1)
       Do while (iRoot.ne.0)
          Write (6,*) ' NODE data'
          Write (6,*) 'NODE @:                         ',iRoot
          Write (6,*) 'Pointer to next NODE:           ',iWork(iRoot  )
          Write (6,*) 'Pointer to stored vector:       ',iWork(iRoot+1)
          If (iWork(iRoot+5).ge.1) Then
             Write (6,*) 'Vector status:                  in Core'
          Else
             Write (6,*) 'Vector status:                  on Disk'
          End If
          Write (6,*) 'Next free position:             ',iWork(iRoot+2)
          Write (6,*) 'Length of vector:               ',iWork(iRoot+3)
          Write (6,*) 'Iteration number:               ',iWork(iRoot+4)
          Write (6,*)
          iRoot=iWork(iRoot)
       End Do
       Write (6,*) '************ End of Status Report *************'
       Write (6,*)
      Return
      End
#ifdef _NOTUSED_
*----------------------------------------------------------------------*
      Subroutine Init_TLists
      Implicit Real*8 (a-h,o-z)
*     Include 'mxdm.fh'
#include <mxdm.fh>

      Logical Triangular
*
      If (DSCF) Then
         Triangular=.True.
         Call Init_TList(Triangular)
         Call Init_PPList
         Call Init_GTList
      End If
*
      Return
      End
#endif
*----------------------------------------------------------------------*
      Subroutine Free_TLists
      Implicit Real*8 (a-h,o-z)

#include "mxdm.fh"
#include "infscf.fh"
*
      If (DSCF) Then
         Call Free_TList
         Call Free_PPList
         Call Free_GTList
      End If
*
      Return
      End
*----------------------------------------------------------------------*
      Subroutine Reduce_Thresholds(EThr_,SIntTh)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
*
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
      Common /Save/ SIntTh_old, EThr_old, DThr_old, DltNTh_old,
     &              FThr_old, ThrInt_Old
*
      Write (6,*)
      Write (6,*) 'Temporary increase of thresholds...'
      Write (6,*)
      SIntTh_old=SIntTh
      EThr_old=EThr
      DThr_old=DThr
      DltNTh_old=DltNTh
      FThr_old=FThr
*
*---- Get threshold used in connection of products of integrals and
*     densities
*
      ThrInt_Old=Get_ThrInt()
*
      EThr=EThr_
      If (EThr_old.eq.Zero) Then
         Relax=One
      Else
         Relax=EThr/EThr_old
      End If
      SIntTh=SIntTh*Relax
      DThr=DThr*Relax
      DltNTh=100.0D0*EThr
      FThr=FThr*Relax
      Call xSet_ThrInt(ThrInt_Old*Relax)
*
      Return
      End
      Subroutine Reset_Thresholds
      Implicit Real*8 (a-h,o-z)
*
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
      Common /Save/ SIntTh_old, EThr_old, DThr_old, DltNTh_old,
     &              FThr_old, ThrInt_Old
*
      Write (6,*)
      Write (6,*) 'Restore thresholds...'
      Write (6,*)
      SIntTh=SIntTh_old
      EThr=EThr_old
      DThr=DThr_old
      DltNTh=DltNTh_old
      FThr=FThr_old
      Call xSet_ThrInt(ThrInt_Old)
*
      Return
      End
