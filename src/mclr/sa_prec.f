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
      Subroutine SA_PREC(S,rdia)
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
      Real*8 S(nroots**2,nroots),rdia(*)

      Do i=1,nroots
!        write(*,*)"rdia,S,CI,ENE"
        Call SA_PREC2(rdia,
     &           S(1,i),Work(ipin(ipci)),ERASSCF(i))
      End Do
      Return
      End

      Subroutine SA_PREC2(rdia,S,CI,ENE)
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "negpre.fh"
#include "WrkSpc.fh"
#include "incdia.fh"
      Real*8 rdia(*),CI(*),S(nroots,nroots)
      Do i=0,nroots-1
       Do j=0,nroots-1
         S(i+1,j+1)=0.0d0
         Do k=1,ncsf(State_Sym)
          dnum=rdia(k)-Ene
          dnum=Sign(Max(Abs(dnum),1.0d-16),dnum)
          S(i+1,j+1)=S(i+1,j+1)+
     &    CI(i*ncsf(State_Sym)+k)*CI(j*ncsf(State_Sym)+k)/dnum
         End Do
       End Do
      End Do
!      Call RecPrt('S',' ',S,nroots,nroots)
      Call MatInvert(S,nroots)
!      Call RecPrt('Inverted S',' ',S,nroots,nroots)
      Return
      end
