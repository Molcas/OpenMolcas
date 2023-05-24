!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003-2005, Valera Veryazov                             *
!               2017, Roland Lindh                                     *
!***********************************************************************
      subroutine SCF(ireturn)
!***********************************************************************
!                                                                      *
!     purpose: perform RHF and UHF calculations                        *
!                                                                      *
!                                                                      *
!***********************************************************************
      Use SCF_Arrays, only: CMO, HDiag, OccNo
      Use Interfaces_SCF, Only: OccDef
      use OFembed, only: Do_OFemb
      use InfSCF, only: DSCF, nDisc, nCore, nD, AufB, nBB, mOV, OnlyProp, iStatPrn, Atom, KSDFT, Name, nnB, Type
      use stdalloc, only: mma_allocate, mma_deallocate
      use Files, only: LuInp
      Implicit None
!
#include "twoswi.fh"
#include "warnings.h"
!
      Integer iReturn

      Character(LEN=8) EMILOOP
      Logical FstItr, Semi_Direct
      Real*8 SIntTh,TCPU1, TCPU2, TWALL1, TWALL2
      Integer iTerm, LUOrb, MemLow, MemSew, LthH
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
      Call CWTime(TCPU1,TWall1)
      Call SCF_Init()
      iTerm=0
!
      Call OpnFls_SCF()
      Call ReadIn_SCF(SIntTh)
      LuOrb=LuInp
!
      Semi_Direct = DSCF .and. (nDisc.ne.0 .or. (nCore.ne.0 .and. nDisc.eq.0))
      If (Semi_Direct) Then
         Call mma_MaxDBLE(MemSew)
         MemLow=Min(MemSew/2,1024*1024)
         MemSew=Max(MemSew/10,MemLow)
         Call xsetmem_ints(MemSew)
      End If
!
      Call Init_SCF()
!                                                                      *
!***********************************************************************
!                                                                      *
!     Pick up the starting orbitals and set the occupation numbers.
!
      LuOrb=LuInp
      Call SOrb(LuOrb,SIntTh,iTerm)
      Call OccDef(OccNo,nnB,nD,CMO,nBB)
!
      Call mma_deallocate(HDiag)
      If (Aufb) Then
         lthH = nBB*nD
      Else
         lthH = mOV
      End If
      Call mma_allocate(HDiag,lthH,Label='HDiag')
!                                                                      *
!***********************************************************************
!                                                                      *
      Call WrInp_SCF(SIntTh)

      Call Cre_SCFWfn()

      FstItr=.True.

      If(.not.OnlyProp) Call WfCtl_SCF(iTerm,KSDFT,FstItr,SIntTh)

      Call Final()
      If (DSCF) Call Free_TLists()
      Call mma_deallocate(Type)
      Call mma_deallocate(Atom)
      Call mma_deallocate(Name)
!
      Call CWTime(TCPU2,TWall2)
!
      Call GMFree()
      Call ClsFls_SCF()
      If (Semi_Direct) Call xRlsMem_Ints()
!
!     Call MolDen Interface
!
      If(nD==1) Then
         Call Molden_Interface(nD-1,'SCFORB','MD_SCF')
      Else
         Call Molden_Interface(nD-1,'UHFORB','MD_SCF')
      End If
      if(iStatPRN.gt.0) then
       Call FastIO('STATUS')
      endif
!
!.... Everything has to come to an end...
!
      iReturn=iTerm
!
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
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
      return
!
      End subroutine SCF
!***********************************************************************
      SubRoutine IniLLs()
!     initialize the diverse linked lists
      use LnkLst, only: LLlist,LLGrad,LLdGrd,LLDelt,LLy,LLx,Init_LLs
      use MxDM, only: MxOptm
      Implicit None

      LLlist=0
      LLGrad=0
      Call IniLst(LLGrad,20)
      Call IniLst(LLDgrd,20)
      Call IniLst(LLDelt,20)
      Call IniLst(LLy,20)
      Call IniLst(LLx,MxOptm)
      Init_LLs=.True.
!
      End subroutine IniLLs
!----------------------------------------------------------------------*
#define _NOTUSED_
#ifdef _NOTUSED_
      Subroutine StatLLS()
      use LnkLst, only: LLGrad,LLdGrd,LLDelt,LLy,LLx,Init_LLs
      Implicit None

      If (Init_LLs) Then
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
      End Subroutine StatLLS
#endif
!----------------------------------------------------------------------*
      SubRoutine KiLLs()
!     dispose the diverse linked lists
      use LnkLst, only: LLGrad,LLdGrd,LLDelt,LLy,LLx,Init_LLs
      Implicit None
      If (Init_LLs) Then
         Call KilLst(LLGrad)
         Call KilLst(LLDgrd)
         Call KilLst(LLDelt)
         Call KilLst(LLy)
         Call KilLst(LLx)
         Init_LLs=.False.
      Else
         Write (6,*) '****** W A R N I N G ! ******'
         Write (6,*) ' Linked list already killed!'
      End If
      Return
      End SubRoutine KiLLs
!----------------------------------------------------------------------*
      Subroutine RclLLs(iDskPt)
      use InfSO, only: MemRsv
      use LnkLst, only: LLGrad,LLdGrd,LLDelt,LLy,LLx,Init_LLs
      use Files, only: LuDel, LuDgd, LuGrd, Lux, Luy
      Implicit None
      Integer iDskPt(5)
      Call RclLst(LLGrad,LuGrd,iDskPt(1),MemRsv)
      Call RclLst(LLDgrd,LuDGd,iDskPt(2),MemRsv)
      Call RclLst(LLDelt,LuDel,iDskPt(3),MemRsv)
      Call RclLst(LLy   ,Lux  ,iDskPt(4),MemRsv)
      Call RclLst(LLx   ,Luy  ,iDskPt(5),MemRsv)
      Init_LLs=.True.
!     Call StatLLs()
      Return
      End Subroutine RclLLs
!----------------------------------------------------------------------*
      Subroutine DmpLLs(iDskPt)
      use LnkLst, only: LLGrad,LLdGrd,LLDelt,LLy,LLx,Init_LLs
      use Files, only: LuDel, LuDgd, LuGrd, Lux, Luy
      Implicit None
      Integer iDskPt(5)
      If (Init_LLs) Then
!        Call StatLLs()
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
      End Subroutine DmpLLs
!----------------------------------------------------------------------*
      Subroutine StlLst(LLink)
      use LnkLst, only: nLList
      Implicit None
      Integer LLink, iRoot
!     return
       Write (6,*)
       Write (6,*) '*********** Status of Linked List *************'
       Write (6,*)
       Write (6,*) ' LLink:',LLink
       Write (6,*)
       Write (6,*) ' CNOD data'
       Write (6,*) 'Error code:                       ',nLList(LLink,0)
       Write (6,*) 'Pointer to first NODE in the list:',nLList(LLink,1)
       Write (6,*) 'Actual length of list:            ',nLList(LLink,2)
       Write (6,*) '# of vectors in core:             ',nLList(LLink,3)
       Write (6,*)
       iRoot=nLList(LLink,1)
       Do while (iRoot.ne.0)
          Write (6,*) ' NODE data'
          Write (6,*) 'NODE @:                         ',iRoot
          Write (6,*) 'Pointer to next NODE:           ',nLList(iRoot,0)
          Write (6,*) 'Pointer to stored vector:       ',nLList(iRoot,1)
          If (nLList(iRoot,5).ge.1) Then
             Write (6,*) 'Vector status:                  in Core'
          Else
             Write (6,*) 'Vector status:                  on Disk'
          End If
          Write (6,*) 'Next free position:             ',nLList(iRoot,2)
          Write (6,*) 'Length of vector:               ',nLList(iRoot,3)
          Write (6,*) 'Iteration number:               ',nLList(iRoot,4)
          Write (6,*)
          iRoot=nLList(iRoot,0)
       End Do
       Write (6,*) '************ End of Status Report *************'
       Write (6,*)
      Return
      End Subroutine StlLst
!----------------------------------------------------------------------*
      Subroutine Free_TLists()
      use InfSCF, Only: DSCF
      Implicit None
!
      If (DSCF) Then
         Call Free_TList()
         Call Free_PPList()
         Call Free_GTList()
      End If
!
      Return
      End Subroutine Free_TLists
!----------------------------------------------------------------------*
      Subroutine Reduce_Thresholds(EThr_,SIntTh)
      use InfSO, only: DltNTh
      use InfSCF, only: EThr, DThr, FThr
      use Save_Stuff, only: DltNTh_old, DThr_Old, EThr_old, FThr_Old, SIntTh_old, ThrInt_old
      use Constants, only: Zero, One, Ten
      Implicit None
      Real*8 EThr_, SIntTh, Relax
      Real*8, External:: Get_ThrInt
!
      Write (6,*)
      Write (6,*) 'Temporary increase of thresholds...'
      Write (6,*)
      SIntTh_old=SIntTh
      EThr_old=EThr
      DThr_old=DThr
      DltNTh_old=DltNTh
      FThr_old=FThr
!
!---- Get threshold used in connection of products of integrals and
!     densities
!
      ThrInt_Old=Get_ThrInt()
!
      EThr=EThr_
      If (EThr_old.eq.Zero) Then
         Relax=One
      Else
         Relax=EThr/EThr_old
      End If
      SIntTh=SIntTh*Relax
      DThr=DThr*Relax
      DltNTh=Ten**2*EThr
      FThr=FThr*Relax
      Call xSet_ThrInt(ThrInt_Old*Relax)
!
      Return
      End Subroutine Reduce_Thresholds
      Subroutine Reset_Thresholds()
      use InfSO, only: DltNTh
      use InfSCF, only: EThr, DThr, FThr
      use Save_Stuff, only: DltNTh_old, DThr_Old, EThr_old, FThr_Old, ThrInt_Old
      Implicit None
!
      Write (6,*)
      Write (6,*) 'Restore thresholds...'
      Write (6,*)
      EThr=EThr_old
      DThr=DThr_old
      DltNTh=DltNTh_old
      FThr=FThr_old
      Call xSet_ThrInt(ThrInt_Old)
!
      Return
      End Subroutine Reset_Thresholds
