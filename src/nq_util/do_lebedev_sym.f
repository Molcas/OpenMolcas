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
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Do_Lebedev_Sym(L_Eff,mPt,ipR)
      Implicit None
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer, Intent(In) :: L_Eff
      Integer, Intent(Out) :: mPt, ipR
      Integer :: mPt_, i, j
      Real*8, Parameter :: Thr = 1.0D-16
      Real*8, Allocatable:: R(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
         Subroutine Do_Lebedev(L_Eff,nPoints,R)
         Implicit None
         Integer L_Eff, nPoints
         Real*8, Allocatable:: R(:,:)
         End Subroutine Do_Lebedev
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      Call Do_Lebedev(L_Eff,mPt_,R)
      mPt=0
      outer: Do i=1,mPt_
        Do j=1,i-1
          If (All(Abs(R(1:3,j)+R(1:3,j)).lt.Thr)) Then
            R(4,i)=0.0D0
            Cycle outer
          End If
        End Do
        mPt=mPt+1
      End Do outer

      Call GetMem('AngRW','Allo','Real',ipR,4*mPt)

      j=1
      Do i=1,mPt_
        If (R(4,i).ne.0.0D0) Then
          Call DCopy_(4,R(:,i),1,Work(ipR+(j-1)*4),1)
          j=j+1
        End if
      End Do
      Call mma_deallocate(R)
      End Subroutine Do_Lebedev_Sym
