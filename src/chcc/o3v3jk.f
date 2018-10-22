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
c
        subroutine o3v3jk (wrk,wrksize,NvGrp,maxdim,LunAux)
c
c
c       This routine do:
c
c       --- Part I - generation of Q and K intermediates
c
c       Q(be,u,i,a) = 2 J(be,u,i,a)   - K(be,u,i,a)
c       Q(be,u,i,a) <-
c              Q1     + [ 2(u,be|a,i) - (u,i|a,be) ]
c              Q2     - [ 2(u,j |a,i) - (u,i|a,j ) ] . t(be,j)
c              Q4     + [ 2(a,i |b,j) - (a,j|b,i ) ] . D(be,b,j,u)
c
c       K(be,u,i,a) <-
c              K1     + [  (u,i |a,be)             ]
c              K2     - [  (u,i |a,j )             ] . t(be,j)
c              K4     - [  (b,i |a,j )             ] . T(be,b,j,u)
c
c       Hoo(i,u)  <<-
c       Hoo1        +  Foo(i,u)
c       Hoo2        +  sum(a,b,j) [2 (ai|bj) - (aj|bi)] Ta(a,b,u,j)
c        implemented as:
c       Hoo2        +  sum(a,b,j) [2 (bi|aj) - (bj|ai)] Ta(a,b,j,u)
c
c       Hvv(a,be) <<-
c       Hvv1        +  Fvv(a,be)
c       Hvv2        +  sum(b,i,j) [2 (ai|bj) - (aj|bi)] Ta(be,b,i,j)
c
c       Hvo(a,i)  <<-
c       Hvo1        +  Fvo(a,i)
c       Hvo2        +  sum(be,j)   [2 (be,j|ai) - (be,i|aj)] T1(be,j)
c
c        Aex(ij,u,v) <-
c       Aex1  + sum(a,b)  [  (ai|bj) . Ta(a,b,u,v)
c
c       A(ij,u,v) <-
c       A1    +   <-      [  (iu,jv)           ]
c       A2    + sum(a)    [  (aj|iu) . t1(a,v) ]
c       A3    + sum(a)    [  (ai|jv) . t1(a,u) ]
c       A4    +           [   Aex(ij,u,v)      ]
c
c        T1(be,u) <<-
c        T161 + sum(a,i)   [ (2(ai|be,u) . t(a,i) ]
c        T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(a,be,i,j) ]
c
c        T18 is implemented in the form:
c        T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(be,a,j,i) ]
c
c        where:
c       D(a,b,i,j) =  t(a,b,j,i) - T(a,b,i,j)
c       T(a,b,i,j) =  t(a,b,i,j)/2 + t(a,i) . t(b,j)
c        Ta(a,b,i,j)=  t(a,b,i,j)   + t(a,i) . t(b,j)
c
c
c1      intermediates Q and K will be temporarry stored in files
c       QFil, KFil as follows
c       do aGrp
c       do beGrp
c         Q(K) (be',u,i,a')
c       end do
c       end do
c
c2      Structure of files, where selected group of (pq|rs) are
c       stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
c
c       (A'I|JK)  I1inxx xx - Group of A'
c
c       (A'I|B'J) I2xxyy xx - Group of A'
c                        yy - Group of B'
c
c       (A'B'|IJ) I3xxyy xx - Group of A'
c                        yy - Group of B'
c
c
c3      Structure of Cholesky vector files
c
c       L1(m,I ,A')  L1vcxx xx - Group of A'
c
c       L2(m,A',B')  L2xxyy xx - Group of A'
c                           yy - Group of B'
c
c4      Structure of Amplitude file
c       t2(A',B',IJ)  T2xxyy xx - Group of A'
c                            yy - Group of B'
c
c
        implicit none
#include "chcc1.fh"
#include "chcc_parcc.fh"
#include "para_info.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
#include "wrk.fh"
#include "chcc_casy.fh"
c
        integer NvGrp,maxdim,LunAux

c       help variables
        integer dim1,dim2,dima,dimb,dimbe
        integer aGrp,bGrp,beGrp,adda,addb,addbe
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4,PossH5
        integer PossK,PossQ
c
        integer PossT

c
        character*6 LunName
c
c       --- introduction part ---
c
cx      Distribute memory
c        @@ yet distribution is the same for all three o3v3 drivers,
c        v buducnosti nosnost usposobit vsetky 3 vlastne
c
c
c
        PossT=PossFree
        call DistMemo3v3jk (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,PossV4,
     c       PossH1,PossH2,PossH3,PossH4,PossH5,
     c       PossK,PossQ,
     c       PossT)
c
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
cmp!           write (6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
       write (6,*) ' Not Enough memory in o3v3jk step!',
     & 'Increase large and/or small segmentation ',
     c                  (1.0d0*PossT)/(1.0d0*wrksize)
           call abend()
        end if
c
c       if (generkey.eq.1) then
cx      vytvorenie suboru T2       vektorov (docasne)
c       call UrobT2 (wrk(PossV1),NbeGrp2,NbGrp2,LunAux)
c       write (6,*) 'T2Vc    done'
c       end if
c
c
c
c##        Vanish Hoo,Hvo,Hvv, T1n, A, Aex (because of paralelization)
        dim1=no*no
        call mv0zero (dim1,dim1,wrk(PossHoo))
        dim1=no*nv
        call mv0zero (dim1,dim1,wrk(PossHvo))
        dim1=nv*nv
        call mv0zero (dim1,dim1,wrk(PossHvv))
        dim1=no*nv
        call mv0zero (dim1,dim1,wrk(PossT1n))
        dim1=no*no*no*(no+1)/2
        call mv0zero (dim1,dim1,wrk(PossA))
        if (intkey.eq.0) then
          call mv0zero (dim1,dim1,wrk(PossAex))
        end if
c
c
c*      cycle over be'
c
        addbe=0
c
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
c
c##        test, if something for this beGrp is planed to be run on this node
        dim1=0
        do dim2=1,NvGrp
          dim1=dim1+BeAID(myRank,beGrp,dim2)
        end do
c       dim1=dim1+BetaID(myRank,beGrp)
        if (dim1.eq.0) goto 11
        if (printkey.gt.1) then ! toto som si nie isty ...
        write (6,*) ' o3v3 JK - ID,beGrp',myRank,beGrp
        end if
cmp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
cmp
c
c
cT1G          vanish H4(be',u)
          dim1=dimbe*no
          call mv0zero (dim1,dim1,wrk(PossH4))
c
cG             Extract H1(be',J) <- T1(be,j)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
c
c
c
c*        cycle over a'
c
          adda=0
c
          do aGrp=1,NvGrp
          dima=DimGrpv(aGrp)
c
c          Test, if this Be' A' combination is to be run on this node
          if (BeAID(myRank,beGrp,aGrp).eq.0) goto 12
          if (printkey.ge.10) then
          write (6,*) ' o3v3 JK - ID,be,a',myRank,beGrp,aGrp
          end if
c
cG             Extract H2(a',J) <- T1(a,j)
          call ExtT1 (wrk(PossH2),wrk(PossT1o),dima,adda)
c
cGx         vanish Q(K)(be',u,i,a')

            dim1=no*no*dimbe*dima
            call mv0zero (dim1,dim1,wrk(PossQ))
            call mv0zero (dim1,dim1,wrk(PossK))
c
cQK2.1.1    read V2(a',o_a,JK) <- I2 (a',o_a|JK)
            LunName=I1name(aGrp)
            dim1=no*(no+1)*no*dima/2
            call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
cQK2.1.2    Expand V1(a',o_a,J,K) <- V2(a',o_a,JK)
            dim1=no*(no+1)/2
                   call Exp2 (wrk(PossV2),wrk(PossV1),dima,no,dim1,no)
c
cQK2.2      map  V4(j,u,i,a') <- V1(a',j,u,i)
            call Map4_4123 (wrk(PossV1),wrk(PossV4),dima,no,no,no)
c
c
            if (aGrp.eq.beGrp) then
c            term A23 only for a'=be'
c
cA23.1        Calc V1(I,JK,L) <- V2(T)(a',I,JK) . H2(a',L)
              dim1=no*no*no*(no+1)/2
              call mv0zero (dim1,dim1,wrk(PossV1))
              dim1=no*no*(no+1)/2
              call mc0c1at3b (dima,dim1,dima,no,dim1,no,
     c                        dim1,dima,no,
     c                        wrk(PossV2),wrk(PossH2),wrk(PossV1))
c
cA23.2f       Add A(ij,u,v) <<- V1(j,iu,v) + V1(i,jv,u)
              dim1=no*(no+1)/2
              call AdV_A23 (wrk(PossV1),wrk(PossA),dim1,no)
c
            end if
c
c
cQK1.1.1    read V2(bea',ui) = (be',a'|IJ)
            if (beGrp.gt.aGrp) then
              LunName=I3name(beGrp,aGrp)
              dim1=no*(no+1)*dima*dimbe/2
            else if (beGrp.eq.aGrp) then
              LunName=I3name(beGrp,aGrp)
              dim1=no*(no+1)*dima*(dima+1)/4
            else
              LunName=I3name(aGrp,beGrp)
              dim1=no*(no+1)*dima*dimbe/2
            end if
            call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
cQK1.1.2    Expand V1(be',a',u,i) <- V2(bea',ui)
            if (beGrp.gt.aGrp) then
              dim1=no*(no+1)/2
              call Exp2 (wrk(PossV2),wrk(PossV1),dimbe,dima,dim1,no)
            else if (beGrp.eq.aGrp) then
              dim1=dima*(dima+1)/2
              dim2=no*(no+1)/2
              call Exp4 (wrk(PossV2),wrk(PossV1),dim1,dima,dim2,no)
            else
              dim1=no*(no+1)/2
              call Exp2i (wrk(PossV2),wrk(PossV1),dima,dimbe,dim1,no)
            end if
c
cQK1.2      map V2(be'u,i,a') <- V1(be',a',u,i)
            call Map4_1423 (wrk(PossV1),wrk(PossV2),dimbe,dima,no,no)
c
cQ1.3       Q(be',u,i,a) <-  - V2(be',u,i,a')
            dim1=no*no*dimbe*dima
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossQ),1,-1.0d0)
c
cK1.3f      K(be',u,i,a) <-  V2(be',u,i,a')
            dim1=no*no*dimbe*dima
            call mv0v1u(dim1,wrk(PossV2),1,wrk(PossK),1,1.0d0)
c
cQ1.4       read V1(be',o_be,a',o_a) = (be'I|a'J)
            LunName=I2name(beGrp,aGrp)
            dim1=no*no*dima*dimbe
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
cQ1.5       Map V2(be',u,i,a') <- V1(be',u,a',i)
            call Map4_1243 (wrk(PossV1),wrk(PossV2),dimbe,no,dima,no)
c
cQ1.6f      Q(be',u,i,a) <-  2 V2(be',u,i,a')
            call mv0v1u(no*no*dimbe*dima,
     c                  wrk(PossV2),1,wrk(PossQ),1,2.0d0)
c
cT161.1     V2(be',u) <- V1(be',u,a',i) . H2(a',i)
            dim1=dimbe*no
            dim2=dima*no
            call mv0zero (dim1,dim1,wrk(PossV2))
               call mv0v1a3u (dim1,dim2,dim2,dim1,
     c                     dim1,dim2,1,1,
     c                     wrk(PossV1),wrk(PossH2),wrk(PossV2))
c
cT161.2f    Add H4(be',u) <<- 2 V2(be',u)
            dim1=dimbe*no
               call mv0v1u (dim1,wrk(PossV2),1,wrk(PossH4),1,2.0d0)
c
cHvo2.1            Make V2(a',i,be',j) <- [2 V1(be',j|a'i) - V1(be',i|a'j)]
            call MkV_Hvo2 (wrk(PossV1),wrk(PossV2),dimbe,dima,no)
c
cHvo2.2     V3(a',i) <- V2(a',i,be',j) . H1(be',j)
            call mv0zero (dima*no,dima*no,wrk(PossV3))
            call mv0v1a3u (dima*no,dimbe*no,dimbe*no,dima*no,
     c                     dima*no,dimbe*no,1,1,
     c                     wrk(PossV2),wrk(PossH1),wrk(PossV3))
c
cHvo2.3f    Add Hvo(a,i) <<- V3(a',i) (da sa spravit aj vlastna rutina)
            call AdT_T17 (wrk(PossHvo),wrk(PossV3),
     c                    nv,dima,no,adda,1.0d0)
c
c
cK2.3       map  V2(j,u,i,a') = -V4(j,u,i,a')
            dim1=no*no*no*dima
            call MkV_K22 (wrk(PossV2),wrk(PossV4),dim1)
c
cQ2.3       make V1(j,u,i,a') = -2V4(i,u,j,a')+V4(j,u,i,a') uz prepermutovanuo
            call MkV_Q22 (wrk(PossV4),wrk(PossV1),dima)
c
cQ2.4f      Q(be',u,i,a') <-  H1(be',j).V1(j,u,i,a')
            call mc0c1a3b (dimbe,no,no,no*no*dima,dimbe,no*no*dima,
     c                     dimbe,no,no*no*dima,
     c                     wrk(PossH1),wrk(PossV1),wrk(PossQ))
c
cK2.4f      K(be',u,i,a') <- H1(be',j).V2(j,u,i,a')
            call mc0c1a3b (dimbe,no,no,no*no*dima,dimbe,no*no*dima,
     c                     dimbe,no,no*no*dima,
     c                     wrk(PossH1),wrk(PossV2),wrk(PossK))
c
cHvv2.1            vanish H3(be',a')
            dim1=dimbe*dima
            call mv0zero(dim1,dim1,wrk(PossH3))
c
c*            cycle over b'
c
              addb=0
c
              do bGrp=1,NvGrp
              dimb=DimGrpv(bGrp)
c
cGx           Extract H5(b',i) <- T1(b,i)
              call ExtT1 (wrk(PossH5),wrk(PossT1o),dimb,addb)
c
c
cQK4.1.12     read V1(be'b',o_be,o_b) <- t2(be',b',I,J)
              LunName=T2name(beGrp,bGrp)
              if (bGrp.eq.beGrp) then
                dim1=dimbe*(dimbe+1)*no*no/2
                call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
cQK4.1.2        Expand V1(be',b',p,q) <- V2(be'b',p,q)
                dim1=dimb*(dimb+1)/2
                call ExpT2 (wrk(PossV2),wrk(PossV1),dimbe,dim1,no)
              else
                dim1=dimbe*dimb*no*no
                call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
              end if
c
c
c        T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(be,a,j,i) ]
              if (aGrp.eq.bGrp) then
c              term T18 only in the case a'=b'
c              preserve: V1(be',a',I,J) = T2(be',a',I,J)
c                        V4(J,K,L,a') = (a',J|K,L)
c              distroy : V2,V3
c
cT18.1                Make V2(a',j,i,u) <- - [2(ai|ju)-(aj|iu)] from V4
                call MkV_T18 (wrk(PossV2),wrk(PossV4),dima,no)
c
cT18.2                Set V3(be',a',j,i) <- V1(be',a',j,i)
                dim1=dimbe*dima*no*no
                call mv0u(dim1,wrk(PossV1),1,wrk(PossV3),1)
c
cT18.3                Make Tau in V3
                   call MkTau_chcc (wrk(PossV3),wrk(PossH1),wrk(PossH2),
     c                           dimbe,dima,no,1.0d0,1.0d0)
c
cT18.4f                Calc H4(be',u) <<- V3(be',a',j,i) . V2(a',j,i,u)
                   dim1=no*no*dima
                     call mc0c1a3b (dimbe,dim1,dim1,no,dimbe,no,
     c                         dimbe,dim1,no,
     c                         wrk(PossV3),wrk(PossV2),wrk(PossH4))
c
              end if
c
c
cQK4.2        read V2(b',o_b,a',o_a) = (b'I|a'J)
              LunName=I2Name(bGrp,aGrp)
              dim1=dimb*dima*no*no
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
              if (aGrp.eq.beGrp) then
c              term A4 only for a'=be' ,If Tau is needed instead
c             of T, A4 term can be joined with Hoo2 subpart
c
cAex1.1         Extract V3(a',b',ij) <- (ai|bj) from V2(b',J,a',I)
                dim1=no*(no+1)/2
                call MkV_A4 (wrk(PossV3),wrk(PossV2),dimb,dima,no,dim1)
c
cAex1.2         Make Tau in V1 (@@@ toto sa menilo oproti T1=0 - OK)
                call MkTau_chcc (wrk(PossV1),wrk(PossH2),wrk(PossH5),
     c                           dima,dimb,no,1.0d0,1.0d0)
c
                if (intkey.eq.0) then
c                cholesky generation of integrals
cAex1.3f          Aex(ij,u,v) <<- V3(T)(a',b',ij) . V1(a',b',u,v)
                  dim1=no*(no+1)/2
                  call mc0c1at3b (dimb*dima,dim1,dimb*dima,no*no,
     c                            dim1,no*no,
     c                            dim1,dima*dimb,no*no,
     c                            wrk(PossV3),wrk(PossV1),wrk(PossAex))
                else
c                W4 and W3 integrals from disc
cAex(A4).3f       A(ij,u,v) <<- V3(T)(a',b',ij) . V1(a',b',u,v)
                  dim1=no*(no+1)/2
                  call mc0c1at3b (dimb*dima,dim1,dimb*dima,no*no,
     c                            dim1,no*no,
     c                            dim1,dima*dimb,no*no,
     c                            wrk(PossV3),wrk(PossV1),wrk(PossA))
                end if
c
cAex1.4ff       nesuspended: V pripade Tau rekonstruovat naspat T2 vo V1
                call MkTau_chcc (wrk(PossV1),wrk(PossH2),wrk(PossH5),
     c                           dima,dimb,no,1.0d0,-1.0d0)
c
              end if
c
cHvv2.2       make Tau(be',b',o_be,o_b) in V1
              call MkTau_chcc (wrk(PossV1),wrk(PossH1),wrk(PossH5),
     c                         dimbe,dimb,no,1.0d0,1.0d0)
c
cHvv2.3       Make V3(b',i,j,a') <-2(a'i|b'j)-(a'j|b'i) from V2(b'J|a'I)
              call MkV_Hvv2 (wrk(PossV3),wrk(PossV2),dima,dimb,no)
c
cHvv2.4       Calc H3(be',a') <<- V1(be',b',i,j) . V3(b',i,j,a')
              dim1=dimb*no*no
              call mc0c1a3b (dimbe,dim1,dim1,dima,dimbe,dima,
     c                       dimbe,dim1,dima,
     c                       wrk(PossV1),wrk(PossV3),wrk(PossH3))
c
              if (aGrp.eq.beGrp) then
c              terms Hoo2 only for a'=be'
c              @ musim tu sachovat, lebo nemam dalsie V a tak
c                je tu jeden lacny perm navyse Hoo2.4 :-(, ak bude
c                vytvorene nove V5, tak  to treba prerobit @
c
cHoo2.1                Map V3(a',b',j,u) <- V1(a',b',u,j)
                call Map3_132 (wrk(PossV1),wrk(PossV3),dima*dimb,no,no)
c
cHoo2.2                Make V1(i,a',b',j) <- 2 V2(b'j|a'i) - V2(b'i|a'j)
                call MkV_Hoo2 (wrk(PossV1),wrk(PossV2),dima,dimb,no)
c
cHoo2.3f        Hoo(i,u) <<- + V1(i,a',b',j) . V3(a',b',j,u)
                dim1=dima*dimb*no
                call mc0c1a3b (no,dim1,dim1,no,no,no,
     c                         no,dim1,no,
     c                         wrk(PossV1),wrk(PossV3),wrk(PossHoo))
cHoo2.4post        Map V1(a',b',u,j) <- V3(a',b',j,u)
                call Map3_132 (wrk(PossV3),wrk(PossV1),dima*dimb,no,no)
c
              end if
c
cQK4.3        make T(be',b',o_be,o_b) in V1 (in V1 is Tau from Hvv2.2)
              call MkTau_chcc (wrk(PossV1),wrk(PossH1),wrk(PossH5),
     c                         dimbe,dimb,no,0.5d0,0.5d0)
c
cQK4.4        map V3(be',o_b,b',o_be) <-V1(be',b',o_be,o_b)(now T in V3)
              call Map4_1342 (wrk(PossV1),wrk(PossV3),dimbe,dimb,no,no)
c
cQK4.5        map V1(b',o_a,o_b,a') <- V2(b',o_b,a',o_a) (now I2 in V1)
              call Map4_1342 (wrk(PossV2),wrk(PossV1),dimb,no,dima,no)
c
cK4.6f        K(be',u,i,a') <<- - V3(be',u_b,b',j_be) . V1(b',j_a,i_b,a')
              call mc0c2a3b (dimbe*no,dimb*no,dimb*no,dima*no,
     c                       dimbe*no,dima*no,
     c                       dimbe*no,dimb*no,dima*no,
     c                       wrk(PossV3),wrk(PossV1),wrk(PossK))
c
cQ4.6         Make D: V2(be',u_b,b',j_be) (now D in V2)
cQ4.6c        from T - V3(be',u_b,b',j_be) and T1 - H1(be',j),H2(b',u)
c             velice specialna procedurka, ale da sa urobit fok :-)))
              call MkD_Q46
     c             (wrk(PossV2),wrk(PossV3),wrk(PossH1),wrk(PossH5),
     c              dimbe,dimb,no)
c
cQ4.7         Make V3(b',o_a,o_b,a') = 2(a',i |b',j) - (a',j|b',i )
cQ4.7c        from V1(b',o_a,o_b,a')
c             dalsia, o nieco menej korenista rutinka fok :-)))
              call MkI_Q47 (wrk(PossV3),wrk(PossV1),dimb,dima,no)
c
cQ4.8f        Q(be',u,i,a') <<- V2(be',u_b,b',j_be) . V3(b',j_a,i_b,a')
              call mc0c1a3b (dimbe*no,dimb*no,dimb*no,dima*no,
     c                       dimbe*no,dima*no,
     c                       dimbe*no,dimb*no,dima*no,
     c                       wrk(PossV2),wrk(PossV3),wrk(PossQ))
c
            addb=addb+dimb
            end do
c
c
cHvv2.5f    Add Hvv(a,be) <<- - H3(be',a')
              call AdH_Hvv2 (wrk(PossH3),wrk(PossHvv),
     c                     dima,dimbe,adda,addbe,nv)
c
cx          write Q and K submatrix to corresponding files
            LunName=Tmp1Name(beGrp,aGrp)
            dim1=dimbe*dima*no*no
            call SaveX (wrk(PossQ),dim1,LunAux,LunName,1,0)
            call SaveX (wrk(PossK),dim1,LunAux,LunName,0,1)
c
12        adda=adda+dima
          end do
c
cT1G        Add T1n(be,u) <<- H4(be',u)
        call AdT_T17(wrk(PossT1n),wrk(PossH4),nv,dimbe,no,addbe,1.0d0)
c
11      addbe=addbe+dimbe
        end do
c
c
#ifdef _MOLCAS_MPP_
c##        Synchronizacny bod:
c        Allreduce Hoo,Hvv,Hvo
c        Allreduce A,Aex
        dim1=no*no
        call gadgop (wrk(PossHoo),dim1,'+')
        dim1=nv*nv
        call gadgop (wrk(PossHvv),dim1,'+')
        dim1=no*nv
        call gadgop (wrk(PossHvo),dim1,'+')
        dim1=no*no*no*(no+1)/2
        call gadgop (wrk(PossA),dim1,'+')
        if (intkey.eq.0) then
          call gadgop (wrk(PossAex),dim1,'+')
        end if
#endif
c
cHoo1.1 Hoo(i,u) <<- Foo(i,u)
        dim1=no*nv
        call mv0v1u (dim1,wrk(PossFoo),1,wrk(PossHoo),1,1.0d0)
c
cHvv1.1 Hvv(a,be) <<- Fvv(a,be)
        dim1=nv*nv
        call mv0v1u (dim1,wrk(PossFvv),1,wrk(PossHvv),1,1.0d0)
c
cHvo1.1 Hvo(a,i) <<- Foo(a,i)
        dim1=no*nv
        call mv0v1u (dim1,wrk(PossFvo),1,wrk(PossHvo),1,1.0d0)
c
cA1.1   read V1(IJ,KL) <- I0(ij,kl)
        LunName=I0Name
        dim1=no*(no+1)*no*(no+1)/4
        call  GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
cA1.2f  Add A(ij,u,v) <<- (iu|jv) from V1(iu,jv)
        dim1=no*(no+1)/2
        call MkV_A1 (wrk(PossA),wrk(PossV1),dim1,no)
c
        if (intkey.eq.0) then
cA4.1f          Add A(ij,u,v) <<- Aex(ij,u,v)
          dim1=no*no*no*(no+1)/2
           call mv0v1u (dim1,wrk(PossAex),1,wrk(PossA),1,1.0d0)
        end if
c
c
c@@
c        call Chck_Hoo (wrk(PossHoo))
c        call Chck_Hvv (wrk(PossHvv))
c        call Chck_Hvo (wrk(PossHvo))
c        call Chck_A (wrk(PossA))
c@@
        return
        end
c
c        ---------------------
c
        subroutine Chck_Hvo (Hvo)
c
c        check Hvo
c
        implicit none
#include "chcc1.fh"
        real*8 Hvo(1:nv,1:no)
c
c        help var
        integer i,a,j,b,bad,tot
        real*8 s
c
        bad=0
        tot=0
c
        do i=1,no
        do a=1,nv
c
          s=0.0d0
c
          do j=1,no
          do b=1,nv
          s=s+(2.0d0*Q21(b,j,a,i)-Q21(b,i,a,j))*T1c(b,j)
          end do
          end do
c
          Hvoc(a,i)=s
c
          if (abs(Hvo(a,i)-s).gt.1.0d-10) then
          bad=bad+1
          end if
          tot=tot+1
c
        end do
        end do
c
        write (6,*) ' Hvo Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Hoo (Hoo)
c
c        check Hoo
c
        implicit none
#include "chcc1.fh"
        real*8 Hoo(1:no,1:no)
c
c        help var
        integer i,u,j,a,b,bad
        real*8 s
c
        bad=0
c
        do i=1,no
        do u=1,no
c
          s=0.0d0
          do j=1,no
          do a=1,nv
          do b=1,nv
          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*
     c        (T2c(a,b,u,j)+T1c(a,u)*T1c(b,j))
          end do
          end do
          end do
c
          Hooc(i,u)=s
c
          if (abs(Hoo(i,u)-s).gt.1.0d-10) then
          bad=bad+1
c          write (6,*) Hoo(i,u),s
          end if
c
        end do
        end do
c
        write (6,*) ' Hoo Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Hvv (Hvv)
c
c        check Hoo
c
        implicit none
#include "chcc1.fh"
        real*8 Hvv(1:nv,1:nv)
c
c        help var
        integer i,j,a,b,be,bad
        real*8 s
c
        bad=0
c
        do be=1,nv
        do a=1,nv
c
          s=0.0d0
          do i=1,no
          do j=1,no
          do b=1,nv
          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*
     c        (T2c(be,b,i,j)+T1c(be,i)*T1c(b,j))
          end do
          end do
          end do
          s=-s
c
          Hvvc(be,a)=s
c
          if (abs(Hvv(a,be)-s).gt.1.0d-10) then
          bad=bad+1
c          write (6,*) Hoo(i,u),s,
          end if
c
        end do
        end do
c
        write (6,*) ' Hvv Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_A (AA)
c
c        check AA(ij,u,v)
c
        implicit none
#include "chcc1.fh"
        real*8 AA(1:no*(no+1)/2,1:no,1:no)
c
c        help var
        integer i,j,ij,u,v,a,b,bad
        real*8 s
c
        bad=0
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
c
        do u=1,no
        do v=1,no
c
          s=Q0(u,i,v,j)
c
          do a=1,nv
           s=s+Q1(a,j,u,i)*T1c(a,v)
          end do
c
          do a=1,nv
           s=s+Q1(a,i,v,j)*T1c(a,u)
          end do
c
          do a=1,nv
          do b=1,nv
           s=s+Q21(a,i,b,j)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
          end do
c
          Ac(i,j,u,v)=s
c
c          write (6,99) i,j,u,v,AA(ij,u,v),s,AA(ij,u,v)-s
c99          format (4(i2,1x),3(f15.10,1x))
          if (abs(AA(ij,u,v)-s).gt.1.0d-10) then
          bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
c
        do u=1,no
        do v=1,no
        do i=2,no
        do j=1,i-1
          Ac(j,i,v,u)=Ac(i,j,u,v)
        end do
        end do
        end do
        end do
c
        write (6,*) ' A   Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Tjedna (T1)
c
c        check T1
c
        implicit none
#include "chcc1.fh"
        real*8 T1(1:nv,1:no)
c
c        help var
        integer u,a,bad
        real*8 s
c
        bad=0
c
        do u=1,no
        do a=1,nv
          s=T1c(a,u)
          if (abs(T1(a,u)-s).gt.1.0d-10) then
          bad=bad+1
          T1(a,u)=s
          end if
        end do
        end do
c
        write (6,*) ' Tjedna   Chck :',bad
c
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Th (T2)
c
c        check
c       T2 = T2(a,b,j_b,u_a)
c
        implicit none
#include "chcc1.fh"
        real*8 T2(1:nv*(nv+1)/2,1:no,1:no)
c
c        help var
        integer u,a,b,j,bad,ab
        real*8 s
c
        bad=0
c
        do u=1,no
        do j=1,no
        ab=0
        do a=1,nv
        do b=1,a
        ab=ab+1
c
          s=T2c(a,b,j,u)+T1c(a,j)*T1c(b,u)
c
          if (abs(T2(ab,j,u)-s).gt.1.0d-10) then
          T2(ab,j,u)=s
          bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' T2  Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_V (VV)
c
c        check  V
c
        implicit none
#include "chcc1.fh"
        real*8 VV(1:nv,1:no,1:no,1:no)
c
c        help var
        integer be,v,u,j,b,bad
        real*8 s
c
        bad=0
c
        do j=1,no
        do u=1,no
        do v=1,no
        do be=1,nv
c
          s=0.0d0
          do b=1,nv
          s=s+Q22(be,b,u,j)*T1c(b,v)
          end do
c
          if (abs(VV(be,v,u,j)-s).gt.1.0d-10) then
          VV(be,v,u,j)=s
          bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' V  Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_T17g (V,dima,adda,dimbe,addbe)
c
c        check V(a',u) -
c        sum(be',b,i)  (b,i|be',a') . [ 2 Ta(b,be',i,u) - Ta(b,be',u,i)]
c
        implicit none
#include "chcc1.fh"
        integer dima,adda,dimbe,addbe
        real*8 V(1:dima,1:no)
c
c        help var
        integer a,u,be,b,i,bad,tot
        real*8 s
c
        tot=0
        bad=0
c
        do u=1,no
        do a=adda+1,adda+dima
c
          s=0.0d0
          do b=1,nv
          do be=1+addbe,addbe+dimbe
          do i=1,no
c
           s=s+Q3(be,a,b,i)*(2.0d0*(T2c(be,b,u,i)+T1c(be,u)*T1c(b,i))
     c                         -(T2c(be,b,i,u)+T1c(be,i)*T1c(b,u)))
          end do
          end do
          end do
c
          if (abs(V(a-adda,u)-s).gt.1.0d-10) then
          bad=bad+1
          end if
          tot=tot+1
c
        end do
        end do
c
        write (6,*) ' T17 Chck :',bad,tot
c
        return
        end
c
