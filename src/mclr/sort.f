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
      Subroutine MMSort(A,B,ldisp1)
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "disp_mclr.fh"
      Real*8 A(*),B(*)
      integer ldisp1(nsym)
      logical geomi,geomj
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      ijG=0
      ijD=0
      iii=0
      Call icopy(nsym,[0],0,ldisp1,1)
      Do iSym=1,nsym
       iG=0
       Do idisp=1,ldisp(isym)
        geomi=iand(ntpert(iii+idisp),16).eq.16
        if (geomi) Then
         ldisp1(isym)=ldisp1(isym)+1
         iG=iG+1
         jG=0
         Do jdisp=1,idisp
          geomj=iand(ntpert(jdisp+iii),16).eq.16
          if (geomj) Then
           jG=jG+1
           ijg1=ijG+itri(ig,jg)
           ijd1=ijD+itri(idisp,jdisp)
           B(ijG1)=A(ijD1)
          End If
         End do
        End If
       End do
       ijG=ijG+iG*(iG+1)/2
       ijD=ijD+ldisp(isym)*(ldisp(isym)+1)/2
       iii=iii+ldisp(isym)
      End do
      return
      end
