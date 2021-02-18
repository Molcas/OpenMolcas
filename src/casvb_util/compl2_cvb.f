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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine compl2_cvb(a,nvec,n,awrk,bwrk,cwrk)
      implicit real*8 (a-h,o-z)
      dimension a(n,n),awrk(n,nvec+n),bwrk(n,n),cwrk(n),dum(1)

      call fmove_cvb(a,awrk,n*nvec)
      call mxunit_cvb(awrk(1,1+nvec),n)
      call schmidt_cvb(awrk,nvec,dum,n,0)
      call schmidtd_cvb(awrk,nvec,awrk(1,nvec+1),n,dum,n,0)
c  Sort N vectors in order of decreasing norms
      do 100 i=1,n
      cwrk(i)=ddot_(n,awrk(1,i+nvec),1,awrk(1,i+nvec),1)
100   continue
      do 200 j=1,n
      cmx=cwrk(1)
      imx=1
      do 300 i=2,n
      if(cwrk(i).gt.cmx)then
        cmx=cwrk(i)
        imx=i
      endif
300   continue
      cwrk(imx)=-DBLE(j)
      call fmove_cvb(awrk(1,imx+nvec),bwrk(1,j),n)
200   continue
      call schmidt_cvb(bwrk,n,dum,n,0)
c  Extract N-NVEC remaining vectors with largest norms
      do 400 i=1,n
      cwrk(i)=ddot_(n,bwrk(1,i),1,bwrk(1,i),1)
400   continue
      do 500 j=1,n-nvec
      cmx=cwrk(1)
      imx=1
      do 600 i=2,n
      if(cwrk(i).gt.cmx)then
        cmx=cwrk(i)
        imx=i
      endif
600   continue
      cwrk(imx)=-DBLE(j)
      call fmove_cvb(bwrk(1,imx),a(1,nvec+j),n)
500   continue
      call nize_cvb(a(1,nvec+1),n-nvec,dum,n,0,0)
      return
      end
