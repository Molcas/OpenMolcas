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
************************************************************************
*                                                                      *
*  Subroutine Shell_MxDens:   returns max density values for each      *
*                             shell pair...                            *
*                                                                      *
************************************************************************
      subroutine Shell_MxDens(Dens,DMax,nSkal)
c----------------------------------------------------------------------
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "shinf.fh"
#include "WrkSpc.fh"
      dimension dmax(nskal,nskal),dens(*)
      ijoff=0
      call fzero(dmax,nskal*nskal)
      Do irp=0,nirrep-1
        ie=0
        Do iSh = 1, nSkal
          n=nbfshl(ish,irp)
          ia=ie+1
          ie=ie+n
          je=0
          Do jSh=1,iSh
            m=nbfshl(jsh,irp)
            ja=je+1
            je=je+m
            Do i=ia,ie
              ij=i*(i-1)/2+ja+ijoff
              Do j=ja,min(i,je)
                dmax(jsh,ish)=max(dmax(jsh,ish),abs(dens(ij)))
                ij=ij+1
              End Do
            End Do
          dmax(ish,jsh)=dmax(jsh,ish)
          End Do
        End Do
        ijoff=ijoff+ie*(ie+1)/2
      End Do
      return
      end
