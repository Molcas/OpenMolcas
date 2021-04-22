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
       subroutine unpackk_ic_3 (i,vint,ndimvi,Vic)
c
c     this routine vint(j,k,l) = <i,j|k,l>
c     for given i from incore (reduced) expanded block Vic
c     ie. symp=symq=symr=syms
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimvi - (norb(symi)) (I)
c     Vic    - incore expanded block of integrals (I)
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimvi
       real*8 vint(1:ndimvi,1:ndimvi,1:ndimvi)
       real*8 Vic(1:(ndimvi*(ndimvi+1)/2)*(1+ndimvi*(ndimvi+1)/2)/2)
c
c     help variables
c
      integer j,k,l,ik,jl,ikjl
c
c
        do k=1,ndimvi
c
c        def ik
        if (i.ge.k) then
          ik=i*(i-1)/2+k
        else
          ik=k*(k-1)/2+i
        end if
c
          jl=0
        do j=1,ndimvi
          do l=1,j
c
c         def jl
          jl=jl+1
c
c           def ikjl
            if (ik.ge.jl) then
            ikjl=ik*(ik-1)/2+jl
            else
            ikjl=jl*(jl-1)/2+ik
            end if
c
            vint(j,k,l)=Vic(ikjl)
            vint(l,k,j)=Vic(ikjl)
c
          end do
        end do
        end do
c
       return
       end
