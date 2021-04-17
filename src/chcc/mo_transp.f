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
        subroutine mo_transp(cmo,cmo_t,no,nv,ndel,nbas)
c
c CMO(p,alpha) <- CMO_t(alpha,p+del),  p=o+v
c
        integer no,nv,nbas,ndel
        integer i,j
        real*8 cmo(1:(no+nv),1:nbas)
        real*8 cmo_t(1:nbas,1:(no+nv+ndel))
c
        do i=1,nbas
        do j=1,(no+nv)
c
        cmo(j,i)=cmo_t(i,j)
        end do
        end do
c
        return
        end
