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
      SubRoutine OutRAS(iKapDisp,iCiDisp)
********************************************************************
*                                                                  *
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
      Real*8, Allocatable:: Kap1(:), Kap2(:), Kap3(:), CIp1(:)
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
      Do iSym=1,nSym
         rsuM=0.0d0
         Call Setup_MCLR(iSym)
         PState_SYM=iEor(State_Sym-1,iSym-1)+1
         nconfM=ncsf(PState_Sym)
         nconf1=ncsf(PState_Sym)
         CI=.false.
         If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.
         If (CI.and.nconf1.eq.1.and.isym.eq.1) CI=.false.
*
*    Allocate areas for scratch and state variables
*
         Call mma_allocate(Kap1,nDens2,Label='Kap1')
         Call mma_allocate(Kap2,nDens2,Label='Kap2')
         Call mma_allocate(Kap3,nDens2,Label='Kap3')
         If (CI) Then
            Call mma_allocate(CIp1,nconfM,Label='CIp1')
            call InCSFSD(Pstate_sym,State_sym,.true.)
         End If
         Do jDisp=1,lDisp(iSym)
            iDisp=iDisp+1
            If (iAnd(ntpert(idisp),2**4).eq.16) Then
               kdisp=DspVec(idisp)
*
               iDisk=iKapDisp(iDisp)
               If (iDisk.ne.-1) Then
                  Len=nDensC
                  Call dDaFile(LuTemp,2,Kap1,Len,iDisk)
                  Call Uncompress(Kap1,Kap3,isym)
                  If (CI) Then
                     ilen=nconfM
                     idis=iCIDisp(iDisp)
                     Call dDaFile(LuTemp,2,CIp1,iLen,iDis)
                  End If
                  Call GASync()
               Else
                  Call GASync()
                  Len=nDensC
                  Call FZero(Kap1,Len)
                  Call GADSum(Kap1,Len)
                  If (CI) Then
                     len=nconfM
                     Call FZero(CIp1,Len)
                     Call GADSum(CIp1,Len)
                  End If
               End If
               Call GASync()
               Call TCMO(Kap3,isym,-1)
               irc=ndens2
               Label='KAPPA   '
               iopt=128
               isyml=2**(isym-1)
               ipert=kdisp
               write(6,'(A,I5," jDisp: ",I5," and iSym:",I5)')
     &           "Writing KAPPA and CI in mclr for iDisp:",
     &           iDisp, jDisp, iSym
               Call dWrMCk(iRC,iOpt,Label,ipert,Kap3,isyml)
               if (irc.ne.0) Call SysAbendMsg('outras','Error in wrmck',
     &              'label=KAPPA')
               irc=nconfM
               iopt=128
               Label='CI      '
               isyml=2**(isym-1)
               ipert=kdisp

               If (iAnd(kprint,8).eq.8)
     &              Write(6,*) 'Perturbation ',ipert
               If (CI) call Guganew(CIp1,0,pstate_sym)
               If (imethod.eq.2.and.(.not.CI).and.nconfM.eq.1)
     &              CIp1(1)=0.0d0
               Call dWrMCk(iRC,iOpt,Label,ipert,CIp1,isyml)
               if (irc.ne.0) Call SysAbendMsg('outras','Error in wrmck',
     &              ' ')
            End If
**********************************************************************
*
         End Do
*
*     Free areas for scratch and state variables
*
*
         If (CI) Call mma_deallocate(CIp1)
         Call mma_deallocate(Kap3)
         Call mma_deallocate(Kap2)
         Call mma_deallocate(Kap1)
      End Do
*
      Return
      End
