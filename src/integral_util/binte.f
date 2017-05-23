************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine binte (k,alfa,beta,r0,a,ggrin,nz)
      Implicit Real*8(A-H,O-Z)
#include "real.fh"
#include "welcom.fh"
      Real*8 grint(0:kmax,kmax), ggrin(nz,0:k,k/2+1,k/4+1),
     &       alfa(nz), a(nz)
*
*     iq = 1
      Call qEnter('Binte')
*     Call RecPrt(' In Binte: Alfa',' ',alfa,nz,1)
*     Call RecPrt(' In Binte: A   ',' ',a   ,nz,1)
      call dcopy_(nz*(k+1)*(k/2+1)*(k/4+1),Zero,0,ggrin,1)
      Do 110 iz=1,nz
         call dcopy_((kMax+1)*kMax,Zero,0,grint,1)
         Call rrint(k,alfa(iz),a(iz),beta,r0,grint,kmax)
         Do 10 i=0,k
            Do 11 j=0,i,2
               j2=j/2+1
               ggrin(iz,i,j2,1)=Zero
               Do 12 l=j,i
                  If (i-l.eq.0) Then
                     ggrin(iz,i,j2,1)=ggrin(iz,i,j2,1)+
     &                    grint(l,j2)*binom(i-j,l-j)
                  Else
                     ggrin(iz,i,j2,1)=ggrin(iz,i,j2,1)+
     &                    grint(l,j2)*a(iz)**(i-l)*binom(i-j,l-j)
                  End If
12             Continue
               ind=1
               Do 13 k2=2,j2-1,2
                  tal=fiint((j-k2)/2,k2/2)/fiint(j/2,0)
                  ind=ind+1
                  ggrin(iz,i,j2,ind)=ggrin(iz,i,j2,1)*tal
13             Continue
11          Continue
10       Continue
110   Continue
*
*     Call RecPrt(' In Binte: Ggrin',' ',
*    &            Ggrin,nz,(k+1)*(k/2+1)*(k/4+1))
      Call qExit('Binte')
      Return
      End
