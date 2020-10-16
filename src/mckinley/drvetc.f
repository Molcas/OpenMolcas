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
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
      Character*8 Lbl
      Real*8 Ccoor(3)
      Real*8, Allocatable:: D0(:), EG(:), Temp(:)

      Ccoor(:)=Zero

      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do

      Call mma_Allocate(D0,nDens,Label='D0')
      Call Get_D1ao_Var(D0,nDens)

      Call mma_allocate(EG,3*nGrad,Label='EG')
      Call Dot1El2(ElGrddot,ElMem,EG,3*nGrad,.true.,CCoor,D0,1)
      Call mma_deallocate(D0)

      EG(:) = -EG(:)

      Call mma_allocate(Temp,3*nGrad,Label='Temp')
      Temp(:)=Zero

      Call Drvel1(Temp)

      EG(:) = EG(:) + Temp(:)

      Lbl='NUCELGR'
      idum=1
      iopt=128
      irc=3*ngrad
      Call dWrMCk(irc,iopt,LBL,idum,Temp,idum)
      If(irc.ne.0) Call SysAbendMsg('drvect',
     &                           'error during write in dwrmck',' ')

      idum=1
      iopt=128
      irc=3*ngrad
      Lbl='DOTELGR'
      Call dWrMCk(irc,iopt,LBL,idum,EG,idum)
      If (irc.ne.0) Call SysAbendMsg('drvect',
     &                           'error during write in dwrmck',' ')
      Call mma_deallocate(EG)
      Call mma_deallocate(Temp)

* needed in RASSI
      loper=0
      Do iCar=1,3

         isym=irrfnc(2**(icar-1))      ! nropr(ichbas(1+iCar))
         Write(Lbl,'(a,i2)') 'ELEC ',iCar
         idcnt=0
         Do iCnttp=1,nCnttp
           Do iCnt=1,dbsc(iCnttp)%nCntr
             idcnt=idcnt+1
             Do idCar=1,3
               Call Cnt1El2(ELGRD,ELMEM,Lbl,idcnt,idcar,loper,
     &                      One,.true.,lbl,0,isym,icar,1)
             End Do
           End Do
         End Do

      End Do

      Return
      End
