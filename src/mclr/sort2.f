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
      Subroutine MMSort2(A,B,P,iel)
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "disp_mclr.fh"
      Real*8 A(*),B(*),P(*)
      integer iel(3)
      logical geomi,geomj
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      ijD=0
      iG=0
      ijG=0
      ijP=0
      iii=0
      Call icopy(3,[0],0,iel,1)
      Do iSym=1,nsym
       Do idisp=1,ldisp(isym)
        geomi=iand(ntpert(idisp+iii),16).eq.16
        if (.not.geomi) Then
         iG=iG+1
         iel(ig)=isym
         Do jdisp=1,ldisp(isym)
          geomj=iand(ntpert(jdisp+iii),16).eq.16
          if (geomj) Then
           ijg=ijg+1
           ijd1=ijD+itri(idisp,jdisp)
           B(ijG)=A(ijD1)
           Else if (idisp.le.jdisp) Then
            ijP=itri(dspvec(jdisp+iii),dspvec(idisp+iii))
            P(ijP)=A(ijD+itri(idisp,jdisp))
          End If
         End do
        End If
       End do
       ijD=ijD+ldisp(isym)*(ldisp(isym)+1)/2
       iii=iii+ldisp(isym)
      End do
      return
      end
