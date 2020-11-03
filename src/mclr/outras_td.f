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
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
       SubRoutine OutRAS_td(iKapDisp,iCiDisp)
********************************************************************
*                                                                  *
* Writes the response to a permenent file MCKINT - Instead of      *
* just a temporary one.                                            *
* Contracts the response coefficient to the hessian                *
*                                                                  *
* Input                                                            *
*       iKapDisp : Disk locations of solutions to respons equation *
*       iCIDisp  : Disk locations of CI Soulutions to response     *
*                                                                  *
* Author: Anders Bernhardsson, 1996                                *
*         Theoretical Chemistry, University of Lund                *
********************************************************************
       Implicit Real*8 (a-h,o-z)
#include "detdim.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "cicisp_mclr.fh"
#include "WrkSpc.fh"
       Character*8 Label
       Integer Pstate_sym
       Integer iKapDisp(nDisp),iCiDisp(nDisp)
       Logical CI
*
*-------------------------------------------------------------------*
*
* Ok construct hessian
*
*-------------------------------------------------------------------*
*
*       Call RecPrt('iKapDisp',' ',Work(iKapDisp),,1)
*
       Write(6,*)
       Write(6,*) '      Writing response to disk in Split guga '//
     &     'GUGA format'
       Write(6,*)
       idisp=0
       Do 100 iSym=1,nSym
          rsuM=0.0d0
          Call Setup_MCLR(iSym)
          PState_SYM=iEor(State_Sym-1,iSym-1)+1
          nconfM=ncsf(PState_Sym)
          nconf1=ncsf(PState_Sym)
          CI=.false.
          If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.
          If (CI.and.nconf1.eq.1.and.isym.eq.1) CI=.false.
          if (Timedep) nconfM=nconfM*2
*
*    Allocate areas for scratch and state variables
*
          Call GetMem('kappa1','Allo','Real',ipKap1,nDens2)
          Call GetMem('kappa2','Allo','Real',ipKap2,nDens2)
          Call GetMem('kappa3','Allo','Real',ipKap3,nDens2)
          If (CI) Then
           Call GetMem('CI1','Allo','Real',ipcip1,nconfM)
           call InCSFSD(Pstate_sym,State_sym,.true.)
          End If
        Do 110 jDisp=1,lDisp(iSym)
         iDisp=iDisp+1
          kdisp=DspVec(idisp)
          iDisk=iKapDisp(iDisp)
          Len=nDensC
*---------------------------------------------------------
* LuTemp temp file in wfctl where the response is written
*---------------------------------------------------------
          Call dDaFile(LuTemp,2,Work(ipKap1),Len,iDisk)
*          if (ndensc.ne.0)
*     &      Call RecPrt('K',' ',Work(ipkap1),ndensc,1)
          Call Uncompress(Work(ipKap1),Work(ipKap3),isym)
          If (CI) Then
           ilen=nconfM
           idis=iCIDisp(iDisp)
           Call dDaFile(LuTemp,2,Work(ipCIp1),iLen,iDis)
*           Call RecPrt(' ',' ',Work(ipCIp1),nconfM,1)
          End If
          Call TCMO(Work(ipKap3),isym,-1)
          irc=ndens2
          Label='KAPPA   '
          iopt=128
          isyml=2**(isym-1)
          ipert=kdisp
          Call dWrMCk(iRC,iOpt,Label,ipert,Work(ipKap3),isyml)
          if (irc.ne.0) Call Abend()
          irc=nconfM
          iopt=128
          Label='CI      '
          isyml=2**(isym-1)
          ipert=kdisp

          If (iAnd(kprint,8).eq.8) Write(6,*) 'Perturbation ',ipert

          If (Timedep.and.CI) then
            Call GetMem('TEMPCI','ALLO','REAL',ipCITmp,nconf1)
            call dcopy_(nconf1,Work(ipCIp1+nconf1),1,Work(ipCITmp),1)
            Call Guganew(Work(ipCITmp),0,pstate_sym)
            Call DSCAL_(nconf1,-1.0d0,Work(ipCITmp),1)
*           Call RecPrt(' ',' ',Work(ipCItmp),nconf1,1)
          End If

          If (CI) call Guganew(Work(ipCIp1),0,pstate_sym)
          If (Timedep) then
            If (CI) Then
*              Call RecPrt(' ',' ',Work(ipCIp1),nconf1,1)
               Call GetMem('TEMPCI2','ALLO','REAL',ipCITmp2,nconfM)

               call dcopy_(nconf1,Work(ipcip1),1,Work(ipCITmp2),1)
               call dcopy_(nconf1,Work(ipcitmp),1,
     &                     Work(ipCITmp2+nconf1),1)
               Call GetMem('TEMPCI','FREE','REAL',ipCITmp,nconf1)
              If (imethod.eq.2.and.(.not.CI).and.nconfM.eq.1)  Then
                  Work(ipcip1)=0.0d0
                  Work(ipcip1+1)=0.0d0
              End If
*             Call RecPrt('M',' ',Work(ipCItmp2),nconfM,1)
              Call dWrMCk(iRC,iOpt,Label,ipert,Work(ipcitmp2),isyml)
              Call GetMem('TEMPCI2','FREE','REAL',ipCITmp2,nconfM)
            End If
          Else
            If (imethod.eq.2.and.(.not.CI).and.nconf1.eq.1)
     &            Work(ipcip1)=0.0d0
            Call dWrMCk(iRC,iOpt,Label,ipert,Work(ipcip1),isyml)
            if (irc.ne.0) Call Abend()
          End If
**********************************************************************
*
 110     Continue
*
*    Free areas for scratch and state variables
*
*
*
          If (CI)  Then
          Call GetMem('CI1','Free','Real',ipcip1,nconfM)
          End If
          Call Getmem('rkappa3','FREE','Real',ipkap3,nDensC)
          Call Getmem('rkappa2','FREE','Real',ipkap2,nDensC)
          Call Getmem('rkappa1','FREE','Real',ipkap1,nDensC)
 100  Continue

*
      Return
      End
