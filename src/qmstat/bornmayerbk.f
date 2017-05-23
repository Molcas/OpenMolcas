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
      Subroutine BornMayerBK(iQ_Atoms,BoMaH,BoMaO)
*
*-- With the Brdarski-Karlstrom scheme, construct the Born-Mayer
*   parameters.
*
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      Dimension BoMaH(MxAt),BoMaO(MxAt)
      Dimension rBdi(MxCen),rBdiQ(MxAt)
      Parameter (cjhr=0.1734d0)

*
*-- The solvent part.
*
      Do 3, i=1,2
        rdi2=0
        Do 4, j=1,3
          rdi2=rdi2+quadi(j,i)
4       Continue
        rBdi(i)=Sqrt(rdi2/charDi(i))
3     Continue

*
*-- The solute part.
*
      Do 5, i=1,iQ_Atoms
        rdi2=0
        Do 6, j=1,3
          rdi2=rdi2+QuadiQ(j,i)
6       Continue
        rBdiQ(i)=Sqrt(rdi2/charDiQ(i))
5     Continue

*
*-- Put together.
*
      Do 7,i=1,iQ_Atoms
        BoMaH(i)=1/(cjhr*(RBdiQ(i)+rbdi(1)))
        BoMaO(i)=1/(cjhr*(RBdiQ(i)+rbdi(2)))
        If(iPrint.ge.8) then
          Write(6,*)'   Born-Mayer parameters.'
          Write(6,'(A,i2,A,2(f12.4))')'    Atom ',i,' (H/O):',BoMaH(i)
     &,BoMaO(i)
        Endif
7     Continue

      Return
      End
