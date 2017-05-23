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
      subroutine two2mean13(carteSO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb)
cbs   gives the two first contributions
cbs   < i M | j M >  with Malpha  and Mbeta
cbs   the other orbit parts cancel
      implicit real*8 (a-h,o-z)
#include "para.fh"
      dimension carteSO(ncontmf,ncontmf,norbsum,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(irun,jrun,icartleft,icartright)
      enddo
      enddo
      enddo
      enddo
c     write(6,*) 'effective integrals'
c     do jrun=1,ncontmf
c     write(6,'(4E20.14)') (onecart(irun,jrun),irun=1,ncontmf)
c     enddo
      return
      end
