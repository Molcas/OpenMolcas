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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      SubRoutine GetDens(FName,Density,iPrint)
************************************************************************
*                                                                      *
* Object: to get the 1 particle density from file INPORB               *
*                                                                      *
* Called from: Drv1el                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OneEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chem. Phys.                       *
*             University of Lund, SWEDEN                               *
*             January 2000                                             *
************************************************************************
      use PrpPnt
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "stdalloc.fh"
      Character Line*80
      Character*(*) FName
      Logical Density
*
      nDen=n2Tri(1)
      nVec=n2Tot
      nOcc=nDim
      If (Density) Call mma_allocate(Den,nDen,label='Den')
      iadDen=1
      Call mma_allocate(Vec,nVec,label='Vec')
      Call mma_allocate(Occ,nOcc,label='Occ')
      iadDen=1
      iadVec=1
      iadOcc=1
      iOpt = 1
      LuVec = 19
      Call RdVec(FName,LuVec,'CO',nIrrep,nBas,nBas,
     &           Vec, Occ, Dummy, iDummy, Line,0,iErr)
      Write (6,*)
      Write (6,'(A)') ' Header from vector file:'
      Write (6,*)
      Write (6,'(A)') Line(:mylen(Line))
      Write (6,*)
*
      If (Density) Then
*
*        Build the density matrix.
*
         call dcopy_(nDen,Zero,0,Den,1)
*
         ictv=iadVec
         icto=iadOcc
         ictd=iadDen
         Do iIrrep=0,nIrrep-1
            Do iBas=1,nBas(iIrrep)
               icv=ictv
               ico=icto
               icd=ictd
               Do j1=0,nBas(iIrrep)-1
                  Do j2=0,j1-1
                     Den(icd)=Den(icd)+Occ(ico)*Vec(icv+j1)
     &                                   *Two  *Vec(icv+j2)
                    icd=icd+1
                  End Do
                  Den(icd)=Den(icd)+Occ(ico)*Vec(icv+j1)
     &                                      *Vec(icv+j1)
                  icd=icd+1
               End Do
               ictv=ictv+nBas(iIrrep)
               icto=icto+1
            End Do
            ictd=ictd+nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
         Call mma_deallocate(Occ)
         Call mma_deallocate(Vec)
         iadVec=iadDen
         iadOcc=iadDen
         nOcc = nDen
         nVec = nDen
         If (iPrint.ge.10) Call PrMtrx(' Density matrix',1,1,iadDen,Den)
*
      End If
*
      Return
      End
