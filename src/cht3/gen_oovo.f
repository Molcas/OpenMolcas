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
        subroutine gen_oovo (w,l0,l1,tmp)
c
c this routine genetates (ij,a,k) integrals from
c blocked MO cholesky vectors
c
c --------
c
c       L0(m,IJ)    L0vctr  I>=J
c       L1(m,I ,A') L1vcxx xx - Group of A'
c

        implicit none
c
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        real*8 tmp(*),l0(*),l1(*),w(*)
        integer a,dima,length,last
c
        integer a_tmp
c
c1        read tmp(m,IJ)
c
        length=nc*(no*(no+1))/2
c
cmp!        write (6,'(A,A6)') 'L0vcrt ','L0vcrt'
cmp!        write (6,*) 'length = ',length
cmp!        write (6,*) 'file size (ifort) = ',8+8*length
c
        call GetX_t3 (tmp,length,LunAux,'L0vctr',1,1)
c2        map l0(IJ,m)   <- tmp(m,IJ)
c
        call Map2_21_t3 (tmp,l0,nc,(no*(no+1)/2))
c
c3        loop over A'
c
        do a=1,NvGrp
cmp@@        dima=nv/NvGrp
        dima=DimGrpaR(a)
c
c4        read tmp(m,I,A')
c
cmp!        write (6,'(A,i3,2x,A6)') 'a,L1Name(a) ',a,L1Name(a)
c
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
c
cmp!         write (6,*) 'dima = ',dima
         length=nc*no*dima
cmp!         write (6,*) 'length = ',length
cmp!         write (6,*) 'file size (ifort) = ',8+8*length
c
        call GetX_t3 (tmp,length,LunAux,L1Name(a),1,1)
c
c5        grow l1(m,I,A)
c
cmp        last=(a-1)*(nv/NvGrp)
c
          last=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          last=last+DimGrpaR(a_tmp)
          end do
        end if
c
        call grow_l1(l1,tmp,dima,nc,no,nv,last)
c
c6        end loop over A'
c
        end do
c
c7        map tmp(m,A,I) <- l1(m,I,A)
c
        call Map3_132_t3 (l1,tmp,nc,no,nv)
c
c7.1        zero w
c
        call zeroma (w,1,((no*(no+1))/2)*nv*no)
c
c8        mult w(IJ,A,I)  <- l0(IJ,m) . tmp(m,A,I)
c
        call mc0c1a3b (
     & (no*(no+1))/2,nc,
     &  nc,nv*no,
     & (no*(no+1))/2,nv*no,
     & (no*(no+1))/2,nc,nv*no,l0,tmp,w)
c
        return
        end
