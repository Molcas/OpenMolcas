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
      SUBROUTINE Cho_VecTransp(Vec,Jin,Jfi,iSym,iRed,iPass)
      Use Para_Info, Only: MyRank, nProcs
      Implicit Real*8 (a-h,o-z)
      Real*8   Vec(*)
      Integer  Jin, Jfi, iSym, iRed, iPass

      Character*13 SecNam
      Parameter (SecNam = 'Cho_VecTransp')

#if defined (_MOLCAS_MPP_)
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "cholq.fh"
#include "choglob.fh"
#include "WrkSpc.fh"
#include "mafdecls.fh"

      Logical LocDbg
#if defined (_DEBUGPRINT_)
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif

      Parameter (N2 = InfVec_N2)

      External  ga_create_irreg, ga_destroy
      Logical   ga_create_irreg, ga_destroy, ok
      Integer   g_a
CVVP:2014 DGA is here
#ifndef _GA_
      Logical   ga_create_local
      Integer   ga_local_woff,nelm,iGAL
      External  ga_local_woff,ga_create_local
#endif
***************************************************************
      map(i) = iWork(ip_map-1+i)
      iDV(i) = iWork(ip_iVecR-1+i)
      nRSL(i) = iWork(ip_nRSL-1+i)
      iAdrLG(i,j) = iWork(ip_iAdrLG-1+MxRSL*(j-1)+i)
      IndRed(i,j) = iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      InfVec_G(i,j,k)=
     & iWork(ip_InfVec_G-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      iL2G(i) = iWork(ip_iL2G-1+i)
***************************************************************

      If (.not.Cho_Real_Par) Then
         If (LocDbg) Then
            Write(6,'(A,A,A)') 'Illegal call to ',SecNam,':'
            Write(6,*)
     &      'Should only be called in parallel, but Cho_Real_Par = ',
     &      Cho_Real_Par
         End If
         Call Cho_Quit('Illegal call to '//SecNam,103)
      End If

      If (iRed .eq. 2) Then
         jRed=3
      Else If (iRed .eq. 3) Then
         jRed=2
      Else
         Call Cho_Quit('iRed must be 2 or 3 in '//SecNam,104)
      End If
      nRS_l = nnBstR(iSym,iRed)   ! local  red set dimension
      nRS_g = nnBstR_G(iSym,iRed) ! global red set dimension
      nV = Jfi - Jin + 1
      nVR = 0

      l_iVecR = nV
      Call GetMem('iVecR','Allo','Inte',ip_iVecR,l_iVecR)
      Call cho_p_distrib_vec(Jin,Jfi,iWork(ip_iVecR),nVR)
      l_VecR = nRS_g*(nVR+1)
      Call GetMem('VecR','Allo','Real',ip_VecR,l_VecR)

      l_nRSL=nProcs
      Call GetMem('RSL','Allo','Inte',ip_nRSL,l_nRSL)
      Call iZero(iWork(ip_nRSL),nProcs)
      iWork(ip_nRSL+MyRank) = nRS_l  ! MyRank starts from 0
      Call Cho_GAIGOP(iWork(ip_nRSL),nProcs,'+')

      MxRSL=nRSL(1)
      Do i=2,nProcs
         MxRSL=max(MxRSL,nRSL(i))
      End Do
      l_iAdrLG=MxRSL*nProcs
      Call GetMem('iAdrLG','Allo','Inte',ip_iAdrLG,l_iAdrLG)

      l_Map = nProcs
      Call GetMem('Map','Allo','Inte',ip_Map,l_Map)
      nProcs_eff = 0
      iStart = 1
      myStart = 0
      Do i=1,nProcs
         If (nRSL(i) .gt. 0) Then
            iWork(ip_Map+nProcs_eff) = iStart
            If ((i-1) .eq. myRank) myStart = iStart
            iStart = iStart + nRSL(i)
            nProcs_eff = nProcs_eff + 1
         End If
      End Do

      If (LocDbg) Then
         Write(LuPri,*)
         Write(LuPri,*) SecNam,': debug info.'
         Write(LuPri,*) '#nodes: ',nProcs,'  myRank: ',myRank
         Write(LuPri,*) '#contributing nodes: ',nProcs_eff
         Write(LuPri,*) 'Symmetry block: ',iSym
         Write(LuPri,*) 'On this node:'
         Write(LuPri,*) 'Vector dimension : ',nRS_l
         Write(LuPri,*) 'Number of vectors: ',nV,' (',Jin,'-',Jfi,')'
         Write(LuPri,*) 'Global vector dimension : ',nRS_g
         Write(LuPri,*) 'Number of global vectors: ',nVR
         Write(Lupri,*) 'MAP:'
         Write(LuPri,*) (map(i),i=1,nProcs_eff)
      End If
CVVP:2014 Local rather than Global
#ifdef _GA_
      ok = ga_create_irreg(mt_dbl,nRS_g,nV,'Ga_Vec',iWork(ip_Map),
     &                     nProcs_eff,1,1,g_a)
#else
      ok = ga_create_local(mt_dbl,nRS_g,nV,'Ga_Vec',g_a)
#endif
      If (.not. ok) Call Cho_Quit(SecNam//': ga_create_irreg error',101)

      If (nRS_l .gt. 0) Then
         myEnd = myStart + nRS_l - 1
#ifdef _GA_
         Call ga_put(g_a,myStart,myEnd,1,nV,Vec,nRS_l)
#else
CVVP:2014 the minimal latency and scalable putC call
         Call ga_putc(g_a,myStart,myEnd,1,nV,Vec,nRS_l)
#endif
      End If
#ifndef _GA_
      nelm=nRS_g*nV
      iGAL=ga_local_woff(g_a)
      Call Cho_GAdGOP(Work(iGAL),nelm,'+')
#else
      Call Cho_GASync()
#endif
      Jin0 = Jin - 1
      Do i=1,nVR
         iv=ip_VecR+nRS_g*(i-1)
         jv=iDV(i) - Jin0
#ifdef _GA_
         Call ga_get(g_a,1,nRS_g,jv,jv,Work(iv),nRS_g)
#else
CVVP:2014 the minimal latency and scalable getC call
         Call ga_getc(g_a,1,nRS_g,jv,jv,Work(iv),nRS_g)
#endif
      End Do

      ok = ga_destroy(g_a)
      If (.not. ok) Call Cho_Quit(SecNam//': ga_destroy error',101)

C --- write the reordered vec on disk

      Call Cho_P_IndxSwp()
      irc=-1
      Call Cho_X_RSCopy(irc,1,jRed)
      If (irc .ne. 0) Then
         Call Cho_Quit(SecNam//
     &                 ': Non-zero return code from Cho_X_RSCopy',
     &                 104)
      End If
      l_mapRS2RS = nnBstR(iSym,1)
      Call GetMem('mapRS2RS','Allo','Inte',ip_mapRS2RS,l_mapRS2RS)
      Call Cho_RS2RS(iWork(ip_mapRS2RS),l_mapRS2RS,jRed,iRed,iPass,iSym)
      Call Cho_P_IndxSwp()

      Call iZero(iWork(ip_iAdrLG),l_iAdrLG)
      Do i = 1,nRS_l
         i1 = IndRed(iiBstR(iSym,iRed)+i,iRed) ! addr in local rs1
         j1 = iL2G(i1) ! addr in global rs1
         j = iWork(ip_mapRS2RS-1+j1-iiBstR_G(iSym,1)) ! addr in glob. rs
         iWork(ip_iAdrLG-1+MxRSL*myRank+i) = j
      End Do
      Call Cho_GAIGOP(iWork(ip_iAdrLG),l_iAdrLG,'+')

      Call GetMem('mapRS2RS','Free','Inte',ip_mapRS2RS,l_mapRS2RS)

      If (LocDbg) Then
         iCount=0
         Do iNode=1,nProcs
            Do iRSL=1,nRSL(iNode)
               iCount=iCount+1
            End Do
         End Do
         If (iCount .ne. nRS_g) Then
            Call Cho_Quit('nRSL error in '//SecNam,104)
         End If
      End If

      iScr=ip_VecR+NRS_g*nVR
      Do j=1,nVR
         iv=ip_VecR+nRS_g*(j-1)
         Call dCopy_(nRS_g,Work(iv),1,Work(iScr),1)
         iCount=0
         iv=iv-1
         Do iNode=1,nProcs
            Do iRSL=1,nRSL(iNode)
               Work(iv+iAdrLG(iRSL,iNode))=Work(iScr+iCount)
               iCount=iCount+1
            End Do
         End Do
      End Do

      If (CHO_ADRVEC .ne. 1) THEN ! only WA files!!
         Call Cho_Quit('CHO_ADRVEC error in '//SecNam,102)
      Else ! write to disk and update InfVec_G(*,3,iSym)
         iVec1=myNumCho(iSym)+1
         lTot=nRS_g*nVR
         If (lTot .gt. 0) Then
            iOpt=1
            iAdr=InfVec_G(iVec1,3,iSym)
            Call dDAfile(LuCho_G(iSym),iOpt,Work(ip_VecR),lTot,iAdr)
         End If
         Do iVec = 1,nVR
            jVec = iVec1 + iVec - 1
            If (jVec .lt. MaxVec) Then
               iWork(ip_InfVec_G+MaxVec*N2*(iSym-1)+MaxVec*2+jVec) =
     &         iWork(ip_InfVec_G+MaxVec*N2*(iSym-1)+MaxVec*2+jVec-1)
     &         + nRS_g
            End If
         End Do
         LastV = myNumCho(iSym) + nVR
         If (LastV .gt. MaxVec) Then
            Call Cho_Quit('Max. number of vectors exceeded in '//SecNam,
     &                    104)
         End If
      End If
      myNumCho(iSym) = myNumCho(iSym) + nVR

C --- deallocations

      Call GetMem('Map','Free','Inte',ip_Map,l_Map)
      Call GetMem('iAdrLG','Free','Inte',ip_iAdrLG,l_iAdrLG)
      Call GetMem('RSL','Free','Inte',ip_nRSL,l_nRSL)
      Call GetMem('VecR','Free','Real',ip_VecR,l_VecR)
      Call GetMem('iVecR','Free','Inte',ip_iVecR,l_iVecR)
#else
      Call Cho_Quit(SecNam//
     &              ' should never be called in serial installation',
     &              103)
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Vec)
         Call Unused_integer(Jin)
         Call Unused_integer(Jfi)
         Call Unused_integer(iSym)
         Call Unused_integer(iRed)
         Call Unused_integer(iPass)
      End If
#endif

      End
