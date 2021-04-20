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
        subroutine gen_vvvo(occ_ind,w3,l1_1,l2_1,tmp)
c
c this routine do
c
c regenerate VVVo integrals from cholesky vectors
c
c -------------------
c
c structure of the cholesky vector files :
c
c       L1(m,I ,A') L1vcxx xx - Group of A'
c
c       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
c                          yy - Group of B'
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer a,b,c,dima,dimb,dimc,occ_ind
        integer length
        integer lasta,lastb,lastc
c
        real*8 w3(1:(nv*(nv+1))/2,1:nv)
        real*8 tmp(*),l1_1(*),l2_1(*)
c
        integer a_tmp,b_tmp,c_tmp
c
c algoritmus je dobry ak maxdim > no
c inak treba vymenit citanie L1 za L2
c
c dalo by sa to urobit podstatne lepsie, kedby
c dircc nevyzadoval VVV ako (ab,c) ale ako (a,b,c)
c
c mozno urob sort L1i (m,c') <- L1(m,i,c')
c ---
c
c1         loop over a'
c
        do a=1,NvGrp
c
c2        loop over b'
c
        do b=1,a
c
c2.1        read L2(m,a',b')
c
         if (a.eq.b) then  ! a=b
c open the pertinent file
c
cmp@        dima=nv/NvGrp
cmp@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
        dima=DimGrpaR(a)
c
cmp!         write (6,*) 'dima = ',dima
         dimb=dima
         length=(dima*(dima+1)*nc)/2
cmp!         write (6,*) 'length L2Name(a,b) = ',L2Name(a,b),length
cmp!     write (6,*) 'file size (g77) = ',16+length*8
cmp!         write (6,*) 'file size L2Name(a,b) (ifort) = ',
cmp!     & L2Name(a,b),8+length*8
c
        call GetX_t3 (tmp,length,LunAux,L2Name(a,b),1,1)
c
         else ! a>b
c open the pertinent file
c
cmp@@        dima=nv/NvGrp
cmp@@        dimb=nv/NvGrp
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
cmp@@         if (b.eq.NvGrp) dimb=nv-((NvGrp-1)*dimb)
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
c
cmp!         write (6,*) 'dima, dimb = ',dima,dimb
         length=dima*dimb*nc
cmp!         write (6,*) 'length L2Name(a,b) = ',L2Name(a,b),length
cmp!     write (6,*) 'file size (g77) = ',16+length*8
cmp!         write (6,*) 'file size L2Name(a,b) (ifort) = ',
cmp!     & L2Name(a,b),8+length*8
c
        call GetX_t3 (tmp,length,LunAux,L2Name(a,b),1,1)
        end if
c
c2.2        map  L2_1(a',b',m) <- tmp(m,a',b')
c
        if (a.eq.b) then ! expand and map
c expand and map l2_1 (a',b',m) <- tmp (m,ab')
        call exMap3_231 (tmp,l2_1,nc,dima)
        else
        call Map3_231_t3 (tmp,l2_1,nc,dima,dimb)
        end if
c
c3         loop over c'
c
        do c=1,NvGrp
c
c3.1        read L1(m,i,c')
c
cmp@@        dimc=nv/NvGrp
cmp@@         if (c.eq.NvGrp) dimc=nv-((NvGrp-1)*dimc)
        dimc=DimGrpaR(c)
c
cmp!         write (6,*) 'dimc = ',dimc
         length=nc*no*dimc
cmp!         write (6,*) 'length L1Name(c) = ',L1Name(c),length
cmp!         write (6,*) 'file size L1Name(c) (ifort) = ',
cmp!     & L1Name(c),8+8*length
c
        call GetX_t3 (tmp,length,LunAux,L1Name(c),1,1)
c
c3.2        extract l1_1 (m,c')_i <- tmp (m,i,c')
c toto by sa dalo nahradit mapovanim
c
        call ext_o_32 (tmp,l1_1,nc,no,dimc,occ_ind)
c
c3.2.1        zero tmp
c
        call zeroma(tmp,1,dima*dimb*dimc)
c
c3.3         mult tmp (a',b',c') <- L2_1 (a',b',m) l1_1 (m,c')
c
               call mc0c1a3b
     & (dima*dimb,nc,nc,dimc,dima*dimb,dimc,
     & dima*dimb,nc,dimc,l2_1,l1_1,tmp)
C
c3.4        add W(ab,c) <- tmp (a'b'c')
c
cmp@@        lasta=(a-1)*(nv/NvGrp)
cmp@@        lastb=(b-1)*(nv/NvGrp)
cmp@@        lastc=(c-1)*(nv/NvGrp)
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.1) then
          do b_tmp=1,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
          lastc=0
        if (c.gt.1) then
          do c_tmp=1,c-1
          lastc=lastc+DimGrpaR(c_tmp)
          end do
        end if
c
cmp!        write (6,'(A,3(i4),2x,3(i4))') 'lasta, lastb, lastc = ',
cmp!     & lasta,lastb,lastc,a,b,c

c sme v gen_vvvo
        call grow_w3 (w3,tmp,
     & nv,nv,dima,dimb,dimc,lasta,lastb,lastc)
c
c3.5        end loop over c'
        end do
c4        end loop over b'
        end do
c5        end loop over a'
        end do
c
        return
        end
