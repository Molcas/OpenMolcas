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
************************************************************************
*                                                                      *
* This routine dumps characters in xml format.                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Subroutine xml_cDump(Name,Appear,Units,Level,Data,nx,ny)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Character*(*) Name
      Character*(*) Appear
      Character*(*) Units
      Character*(*) Data(*)
      Integer       nx
      Integer       ny
      Integer       Level
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer ix
      Integer iy
      Integer ind
      Integer opt
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      If(ny.eq.1.and.nx.lt.5) Then
         Call xml_cDumpa(Name,Len(Name),
     &                   Appear,Len(Appear),
     &                   Units,Len(Units),Level,
     &                   nx,ny,0)
         Do ix=1,nx
            Call xml_cDumpb(Data(ix),Len(Data(ix)),0)
         End Do

      Else
         Call xml_cDumpa(Name,Len(Name),
     &                   Appear,Len(Appear),
     &                   Units,Len(Units),Level,
     &                   nx,ny,1)
         Do iy=1,ny
            Do ix=1,nx
               ind=(ix-1)*ny+iy
               If((ix/10)*10.eq.ix.or.ix.eq.nx) Then
                  opt=1
               Else
                  opt=0
               End If
               Call xml_cDumpb(Data(ind),Len(Data(ind)),opt)
            End Do
         End Do
      End If
      Call xml_cDumpc(Name,Len(Name))
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*

      Return
      End
