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
       subroutine unpackk_ic_2 (i,vint,ndimvi,ndimvj,Vic)
c
c     this routine vint(j,k,l) = <i,j|k,l>
c     for given i from incore (reduced) expanded block Vic
c     ie. symp=symr,symq=syms
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimvi - (norb(symi)) (I)
c     ndimvj - (norb(symj)) (I)
c     Vic    - incore expanded block of integrals (I)
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimvi,ndimvj
       real*8 vint(1:ndimvj,1:ndimvi,1:ndimvj)
       real*8 Vic(1:(ndimvi*(ndimvi+1)/2),1:(ndimvj*(ndimvj+1)/2))
c
c     help variables
c
      integer j,k,l,ik,jl
c
c
        do k=1,ndimvi
        if (i.ge.k) then
          ik=i*(i-1)/2+k
        else
          ik=k*(k-1)/2+i
        end if
          jl=0
          do j=1,ndimvj
          do l=1,j
            jl=jl+1
            vint(j,k,l)=Vic(ik,jl)
            vint(l,k,j)=Vic(ik,jl)
          end do
          end do
        end do
c
       return
       end
