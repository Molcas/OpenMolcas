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
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Drvetc(ngrad)
************************************************************************
*                                                                      *
* Object: driver for computation of gradients of one-electron matrices.*
*                                                                      *
* Called from: Mckinley                                                *
*                                                                      *
* Calling    :                                                         *
*              GetMem                                                  *
*              Cnt1El2                                                 *
*                                                                      *
*             Written by Anders Bernhardsson for electric field        *
*             Gradients                                                *
*             October  97                                              *
************************************************************************
      use Basis_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External ElGrd,elgrddot
      External ElMem
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "disp.fh"
      Character*8 Lbl
      Real*8 Ccoor(3)
      Call QEnter('DrvEtc')
      iRc=-1
      iOpt=1
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
      Call Get_D1ao_Var(ipD0,Length)
      If ( length.ne.nDens ) Then
         Write (6,*) 'DrvEtc: length.ne.nDens'
         Write (6,*) 'Length,nDens=',Length,nDens
         Call QTrace()
         Call Abend()
      End If
*     Write(*,*) Ddot_(ndens,Work(ipD0),1,Work(ipD0),1)
      call dcopy_(3,[0.0d0],0,CCOOR,1)
      ncomp=1
      loper=0
      Call GetMem('ELEGRD','ALLO','REAL',ipEG,3*ngrad)
      Call Dot1El2(ElGrddot,ElMem,Work(ipEG),3*nGrad,
     &                 .true.,CCoor,
     &                 Work(ipD0),1)
      CALL DSCAL_(3*ngrad,-1.0d0,Work(ipEG),1)
      Call GetMem('D0  ','FREE','Real',ipD0,nDens)
      Call GetMem('TEMP','ALLO','Real',ipTemp,3*ngrad)
      Call GetMem('TEMP','CHEC','Real',ipTemp,3*ngrad)
      call dcopy_(3*ngrad,[0.0d0],0,Work(ipTemp),1)
      Call Drvel1(Work(ipTemp))
      Call DaXpY_(3*ngrad,1.0d0,Work(ipTemp),1,Work(ipEG),1)
      Lbl='NUCELGR'
      idum=1
      iopt=128
      irc=3*ngrad
      Call dWrMCk(irc,iopt,LBL,idum,Work(ipTemp),idum)
      if(irc.ne.0)
     & Call SysAbendMsg('drvect','error during write in dwrmck',' ')
      idum=1
      iopt=128
      irc=3*ngrad
      Lbl='DOTELGR'
      Call dWrMCk(irc,iopt,LBL,idum,Work(ipEG),idum)
      if (irc.ne.0)
     & Call SysAbendMsg('drvect','error during write in dwrmck',' ')
      Call GetMem('ELGR','FREE','REAL',ipeg,3*ngrad)
      Call GetMem('TEMP','FREE','Real',ipTemp,3*ngrad)
* needed in RASSI
      Do iCar=1,3
      isym=irrfnc(2**(icar-1))! nropr(ichbas(1+iCar))
      Write(Lbl,'(a,i2)') 'ELEC ',iCar
      idcnt=0
      Do iCnttp=1,nCnttp
        Do iCnt=1,dbsc(iCnttp)%nCntr
          idcnt=idcnt+1
          Do idCar=1,3
            Call Cnt1El2(ELGRD,ELMEM,Lbl,idcnt,idcar,loper,
     &               1.0d0,.true.,
     &               lbl,0,isym,icar,1)
          End Do
        End Do
      End Do
      End Do
      Call QExit('DrvEtc')
      Return
      End
