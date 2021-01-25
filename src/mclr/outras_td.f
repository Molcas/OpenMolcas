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
#include "stdalloc.fh"
       Character(LEN=8) Label
       Integer Pstate_sym
       Integer iKapDisp(nDisp),iCiDisp(nDisp)
       Logical CI
       Real*8, Allocatable:: Kap1(:), Kap2(:), Kap3(:), CIp1(:,:)
*
*-------------------------------------------------------------------*
*
* Ok construct hessian
*
*-------------------------------------------------------------------*
*
       Write(6,*)
       Write(6,*) '      Writing response to disk in Split guga '//
     &     'GUGA format'
       Write(6,*)
       idisp=0
       Do 100 iSym=1,nSym
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
          Call mma_allocate(Kap1,nDens2,Label='Kap1')
          Call mma_allocate(Kap2,nDens2,Label='Kap2')
          Call mma_allocate(Kap3,nDens2,Label='Kap3')
          If (CI) Then
             If (TimeDep) Then
                Call mma_allocate(CIp1,nconf1,2,Label='CIp1')
             Else
                Call mma_allocate(CIp1,nconf1,1,Label='CIp1')
             End If
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
          Call dDaFile(LuTemp,2,Kap1,Len,iDisk)
*         if (ndensc.ne.0) Call RecPrt('K',' ',Kap1,ndensc,1)
          Call Uncompress(Kap1,Kap3,isym)
          If (CI) Then
             ilen=nconfM
             idis=iCIDisp(iDisp)
             Call dDaFile(LuTemp,2,CIp1,iLen,iDis)
*            Call RecPrt(' ',' ',CIp1,nconfM,1)
          End If
          Call TCMO(Kap3,isym,-1)
          irc=ndens2
          Label='KAPPA   '
          iopt=128
          isyml=2**(isym-1)
          ipert=kdisp
          Call dWrMCk(iRC,iOpt,Label,ipert,Kap3,isyml)
          if (irc.ne.0) Call Abend()
          irc=nconfM
          iopt=128
          Label='CI      '
          isyml=2**(isym-1)
          ipert=kdisp

          If (iAnd(kprint,8).eq.8) Write(6,*) 'Perturbation ',ipert

          If (Timedep.and.CI) then
            Call Guganew(CIp1(:,2),0,pstate_sym)
            Call DSCAL_(nconf1,-1.0d0,CIp1(:,2),1)
          End If

          If (CI) call Guganew(CIp1(:,1),0,pstate_sym)
          If (Timedep) then
            If (CI) Then
*              Call RecPrt(' ',' ',CIp1,nconfM,1)
               Call dWrMCk(iRC,iOpt,Label,ipert,CIp1,isyml)
            End If
          Else
            If (imethod.eq.2.and.(.not.CI).and.nconf1.eq.1)
     &         CIp1(1,1)=0.0d0
            Call dWrMCk(iRC,iOpt,Label,ipert,CIp1,isyml)
            if (irc.ne.0) Call Abend()
          End If
**********************************************************************
*
 110     Continue
*
*    Free areas for scratch and state variables
*
          If (CI) Call mma_deallocate(CIp1)
          Call mma_deallocate(Kap3)
          Call mma_deallocate(Kap2)
          Call mma_deallocate(Kap1)
 100  Continue

*
      Return
      End
