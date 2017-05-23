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
      Subroutine ExpHinvv(rdia,v,u,alpha,beta)
*
*     Preconditioning of the state transfer part
*     of the  electronic hessian with an subunit
*     described with the Explicit hessian and
*     the rest with the diagonal
*
*                              -1
*  |u> = alpha|u> + beta  (H-E ) |v>
*                           0 0
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
#include "incdia.fh"
      Real*8 v(*),u(*),rdia(*)
*
      If (nexp.ne.0) Then
      Call GetMem('Tmp1','ALLO','REAL',ipTmp1,nExp)
      Call GetMem('Tmp4','ALLO','REAL',ipTmp4,nExp)
*
      Do i=0,nExp-1
       j=iwork(iplst+i)
       Work(ipTmp1+i)=v(j)
       Work(ipTmp4+i)=u(j)
      End Do
*
      irc=0
      call dgetrs_('N',NEXP,1,Work(iphx),nexp,
     &                iwork(ipvt),Work(ipTmp1),nexp,irc)
      If (irc.ne.0) then
       write(6,*) 'Error in DGETRS called from exphinvv'
       Call Abend
      endif
*
      If (alpha.eq.0.0d0.and.beta.eq.1.0d0) Then
      Call DVEM(nConf1,v,1,rdia,1,u,1)
      Else If (alpha.eq.0.0d0) Then
      Do i=1,nConf1
        u(i)=beta*rDia(i)*v(i)
      End Do
      Else If (alpha.eq.1.0d0) Then
      Do i=1,nConf1
        u(i)=u(i)+beta*rDia(i)*v(i)
      End Do
      else
      Do i=1,nConf1
        u(i)=alpha*u(i)+beta*rDia(i)*v(i)
      End Do
      End If
*
      Do i=0,nExp-1
       j=iwork(iplst+i)
       u(j)=alpha*Work(iptmp4+i)+beta*Work(ipTmp1+i)
      End Do
      Call GetMem('Tmp1','FREE','REAL',ipTmp1,nExp)
      Call GetMem('Tmp4','FREE','REAL',iptmp4,nExp)

      Else
      If (alpha.eq.0.0d0.and.beta.eq.1.0d0) Then
      Call DVEM(nConf1,v,1,rdia,1,u,1)
      Else If (alpha.eq.0.0d0) Then
      Do i=1,nConf1
        u(i)=beta*rDia(i)*v(i)
      End Do
      Else If (alpha.eq.1.0d0) Then
      Do i=1,nConf1
        u(i)=u(i)+beta*rDia(i)*v(i)
      End Do
      else
      Do i=1,nConf1
        u(i)=alpha*u(i)+beta*rDia(i)*v(i)
      End Do
      End If

      End If
*
      Return
      End
