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
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine UpOrb(n,o,w,Orb,Lab)
      Implicit Real*8 (a-h,o-z)
      Dimension Orb(n)
#include "parm.fh"
#include "common.fh"
      Character*(LENIN8) Lab(n)
      dimension indx((MxLqn+1)*(MxLqn+1))
      dimension jndx((MxLqn+1)*(MxLqn+1))
#include "symlab.fh"
      Do 100 i=1,(MxLqn+1)*(MxLqn+1)
         indx(i)=0
         jndx(i)=0
100   Continue
      Do 200 iBas=1,n
*        Write(*,'(a,i5)') '     iBas:',iBas
         If(Lab(iBas)(1:LENIN).eq.Center) Then
            iBlk=0
            Do 210 i=1,(MxLqn+1)*(MxLqn+1)
               If(type(i).eq.Lab(iBas)(LENIN1:LENIN8)) iBlk=i
210         Continue
*           Write(*,'(a,i5)') '       iBlk:',iBlk
            If(iBlk.eq.0) Then
               Write(6,*) 'Unknown basis function: ',
     &                    Lab(iBas)(LENIN1:LENIN8)
               Call Quit_OnUserError()
            End If
            indx(iBlk)=indx(iBlk)+1
*           Write(*,'(a,i5)') '     indx:',indx(iBlk)
            jndx(iBlk)=0
            Do 220 jBas=1,iBas
*              Write(*,'(a,i5)') '       jBas:',jBas
               If(Lab(jBas)(1:LENIN).eq.Center) Then
                  If(Lab(iBas)(LENIN1:LENIN8).eq.
     &               Lab(jBas)(LENIN1:LENIN8)) Then
                     jndx(iBlk)=jndx(iBlk)+1
*                    Write(*,'(a,i5)') '     jndx:',jndx(iBlk)
                     ind=jndx(iBlk)
     &                  +indx(iBlk)*(indx(iBlk)-1)/2
     &                  +iSymBk(iBlk)-1
*                    add=w*o*Cmo(iBas)*Cmo(jBas)
                     add=w*o*Orb(iBas)*Orb(jBas)
                     pDsym(ind)=pDsym(ind)+add
*                    Write(*,'(a,f12.6,a,i5)') 'add',add,' to',ind
                  End If
               End If
220         Continue
         End If
200   Continue
      Return
      End
