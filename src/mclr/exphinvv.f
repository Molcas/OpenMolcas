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
      use Exp, only: H0S, H0F, SBIDT, nExp
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
#include "stdalloc.fh"
#include "incdia.fh"
      Real*8 v(*),u(*),rdia(*)
      Real*8, Allocatable:: Tmp1(:), Tmp4(:)
*
      If (nExp.ne.0) Then
      Call mma_allocate(Tmp1,nExp,Label='Tmp1')
      Call mma_allocate(Tmp4,nExp,Label='Tmp4')
*
      Do i=1,nExp
       j=SBIDT(i)
       Tmp1(i)=v(j)
       Tmp4(i)=u(j)
      End Do
*
      irc=0
      call dgetrs_('N',NEXP,1,H0S,nexp,H0F,Tmp1,nexp,irc)
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
      Do i=1,nExp
       j=SBIDT(i)
       u(j)=alpha*Tmp4(i)+beta*Tmp1(i)
      End Do
      Call mma_deallocate(Tmp1)
      Call mma_deallocate(Tmp4)

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
