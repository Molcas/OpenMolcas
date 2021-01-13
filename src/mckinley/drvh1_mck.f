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
* Copyright (C) 1991, Roland Lindh                                     *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Drvh1_mck(nGrad,Nona)
************************************************************************
*                                                                      *
* Object: driver for computation of gradients of one-electron matrices.*
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
*             Modified by Anders Bernhardsson for Gradients            *
*             May 95                                                   *
************************************************************************
      Use Basis_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External OvrGrd_mck,KneGrd_mck,nagrd_mck,prjgrd_mck,m1grd_mck ,
     &         srogrd_mck, nona2
      External OvrMem_mck,KneMem_mck,namem_mck,prjmm1,m1mm1, na2mem,
     &         sromm1
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Character*8 Label
      Logical Nona, lECP
      Real*8, Allocatable:: Fock(:), D0(:)
*
      If (show) Then
         nFock = 0
         nDens = 0
         Do 1 iIrrep = 0, nIrrep - 1
            nFock = nFock + nBas(iIrrep)*(nBas(iIrrep)+1)/2
            nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
 1       Continue
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
         Call mma_allocate(D0,nDens,Label='D0')
         Call Get_D1ao_Var(D0,nDens)
*...  Read the generalized Fock matrix
*...  Fock matrix in AO/SO basis
         Call mma_allocate(Fock,nFock,Label='Fock')
         Call Get_Fock_Occ(Fock,nFock)
      Else
         nFock = 1
         nDens = 1
         Call mma_allocate(Fock,nFock,Label='Fock')
         Call mma_allocate(D0,nDens,Label='D0')
         Fock(1)=Zero
         D0(1)=Zero
      End If
      If (Nona) Then
************************************************************************
*0a)                                                                   *
*                                                                      *
*     Antisymmetric gradient of Overlap matrix                         *
*                                                                      *
************************************************************************
          Label='OVRGRDA '
          idcnt=0
          Do iCnttp=1,nCnttp
             Do iCnt=1,dbsc(iCnttp)%nCntr
                idcnt=idcnt+1
                Do idCar=1,3
            Call Cnt1El(OvrGrd_mck,OvrMem_mck,Label,idcnt,idcar,loper,
     &              -One,.false.,Fock,'OVRGRDA ',0)
                End Do
             End Do
          End Do
************************************************************************
*0b)                                                                   *
*                                                                      *
*     Non-adiabatic second derivative integrals                        *
*                                                                      *
************************************************************************
          Label='NONA2   '
          idcnt=0
          Do iCnttp=1,nCnttp
             Do iCnt=1,dbsc(iCnttp)%nCntr
                idcnt=idcnt+1
                Do idCar=1,3
            Call Cnt1El(NONA2,NA2Mem,Label,idcnt,idcar,loper,
     &                  One,.false.,Fock,'NONA2   ',0)
                End Do
             End Do
          End Do
*
      End If
*
************************************************************************
*1)                                                                    *
*                                                                      *
*     Gradient of Overlap matrix                                       *
*                                                                      *
************************************************************************
      Label='OVRGRD  '
      idcnt=0
      Do iCnttp=1,nCnttp
        Do iCnt=1,dbsc(iCnttp)%nCntr
          idcnt=idcnt+1
          Do idCar=1,3
            Call Cnt1El(OvrGrd_mck,OvrMem_mck,Label,idcnt,idcar,loper,
     &               One,.false.,Fock,'OVRGRD  ',0)
          End Do
        End Do
      End Do
*
************************************************************************
*2)                                                                    *
*                                                                      *
*     Gradient of Kinetic operator                                     *
*                                                                      *
*                                                                      *
************************************************************************
      Label='KNEGRD  '
      idcnt=0
      Do iCnttp=1,nCnttp
        Do iCnt=1,dbsc(iCnttp)%nCntr
          idcnt=idcnt+1
          Do idCar=1,3
            Call Cnt1El(KneGrd_mck,KneMem_mck,Label,idcnt,idcar,loper,
     &                  One,.false.,D0,'ONEGRD  ',0)
          End Do
        End Do
      End Do
*
************************************************************************
*3)                                                                    *
*                                                                      *
*     Gradient of Nuclear attraction Operator                          *
*                                                                      *
*                                                                      *
************************************************************************
      Label='NAGRD   '
      idcnt=0
      Do iCnttp=1,nCnttp
        Do iCnt=1,dbsc(iCnttp)%nCntr
          idcnt=idcnt+1
          Do idCar=1,3
            Call Cnt1El(NaGrd_mck,NaMem_mck,Label,idcnt,idcar,loper,
     &                  One,.true.,D0,'ONEGRD  ',1)
          End Do
        End Do
      End Do
*
*
************************************************************************
*3)                                                                    *
*                                                                      *
*     Gradient of Nuclear attraction Operator  ECP-part                *
*                                                                      *
*                                                                      *
************************************************************************
      lECP = .False.
      DO i = 1, nCnttp
         lECP = lECP .or. dbsc(i)%ECP
      End Do
      If (lecp) Then
      idcnt=0
      Do iCnttp=1,nCnttp
        Do iCnt=1,dbsc(iCnttp)%nCntr
          idcnt=idcnt+1
          Do idCar=1,3
            Label='PRJGRD  '
            Call Cnt1El(Prjgrd_mck,PrjMm1,Label,idcnt,idcar,loper,
     &               One,.true.,D0,'ONEGRD  ',1)
            Label='M1GRD  '
            Call Cnt1El(m1grd_mck,m1Mm1,Label,idcnt,idcar,loper,
     &               One,.true.,D0,'ONEGRD  ',1)
            Label='SROGRD  '
            Call Cnt1El(Srogrd_mck,sroMm1,Label,idcnt,idcar,loper,
     &               One,.true.,D0,'ONEGRD  ',1)
          End Do
        End Do
      End Do
      End if
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(D0)
      Call mma_deallocate(Fock)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nGrad)
      End
