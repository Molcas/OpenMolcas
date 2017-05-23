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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
      Subroutine GetGap(Eorb,nData,nAufb,Gap,Efermi)
************************************************************************
*                                                                      *
* This routine figure out the homo lumo gap.                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Integer nData
      Integer nAufb
      Real*8  Eorb(nData)
      Real*8  Gap
      Real*8  Efermi
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer i,j,k
      Real*8  tmp
*----------------------------------------------------------------------*
* Setup                                                                *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Sort array                                                           *
*----------------------------------------------------------------------*
*     Write(6,*) 'Unsorted orbitals energies'
*     Write(6,'(10F12.6)') Eorb
      Do i=1,nData-1
         k=i
         Do j=i+1,nData
            If(Eorb(j).lt.Eorb(k)) k=j
         End Do
         tmp=Eorb(k)
         Eorb(k)=Eorb(i)
         Eorb(i)=tmp
      End Do
*     Write(6,*) 'Sorted orbitals energies'
*     Write(6,'(10F12.6)') Eorb
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      If(nAufb.le.0) Then
         Gap=1.0d3
         Efermi=Eorb(1)
      Else If(nAufb.ge.nData) Then
         Gap=1.0d3
         Efermi=Eorb(nData)+1.0d-3
      Else
         Gap=Eorb(nAufb+1)-Eorb(nAufb)
         Efermi=0.5d0*(Eorb(nAufb+1)+Eorb(nAufb))
      End If
*     Write(6,*) 'Gap:',Gap
*     Write(6,*) 'Efermi:',Efermi
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
