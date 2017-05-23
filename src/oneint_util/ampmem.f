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
* Copyright (C) 1996, Per Ake Malmqvist                                *
*               1996, Roland Lindh                                     *
************************************************************************
      Subroutine AMPMem(nHer,MemAMP,la,lb,lr)
C     Statement function for Cartesian index
      nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2

C Mem1: Workspace for MltPrm.
C Mem2: Tables Tpp,Tp,T0,Tm, and Tmm.
C Mem3: Result from AMPr.
      Mem1=0
      Call MltMmP(nOrder,Mem,la,lb+2,2)
      Mem1=max(Mem1,Mem)
      Mem2=6*nElem(la)*nElem(lb+2)
      nHer = nOrder
      Call MltMmP(nOrder,Mem,la,lb+1,1)
      Mem1=max(Mem1,Mem)
      Mem2=Mem2+3*nElem(la)*nElem(lb+1)
      Call MltMmP(nOrder,Mem,la,lb  ,2)
      Mem1=max(Mem1,Mem)
      Mem2=Mem2+6*nElem(la)*nElem(lb)
      If (lb.ge.1) Then
        Call MltMmP(nOrder,Mem,la,lb-1,1)
        Mem1=max(Mem1,Mem)
        Mem2=Mem2+3*nElem(la)*nElem(lb-1)
        If (lb.ge.2) Then
          Call MltMmP(nOrder,Mem,la,lb-2,2)
          Mem1=max(Mem1,Mem)
          Mem2=Mem2+6*nElem(la)*nElem(lb-2)
        End If
      End If
      Mem3=6*nElem(la)*nElem(lb)

      MemAMP=Mem1+Mem2+Mem3+1

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
