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
        subroutine o3v3chol (wrk,wrksize,NvGrp,maxdim,LunAux)
c
c       This routine do:
c
c       --- Part II - generation of Q3,K3,Gvv,Goo,T26-9 in X form
c
c       T(be,u)      <<-
c       T162 - sum(b,i)   [ (iu|be,b) . t(b,i) ]
c       T17  + sum(a,b,i) [ 2(a,i|b,be) - (b,i|a,be) ] . Ta(a,b,i,u)
c          = + sum(a,b,i)  (a,i|b,be) . [ 2 Ta(a,b,i,u) - Ta(a,b,u,i) ]
c
c       N.B. T17 implementede in expression:
c       T(a,u)      <<-
c       T17  + sum(be,b,i)  (b,i|be,a) . [ 2 Ta(b,be,i,u) - Ta(b,be,u,i) ]
c
c
c       Q(be,u,i,a) <-
c              Q3     +   [ 2(b,be|a,i) - (b,i|a,be) ] . t(b,u )
c
c       K(be,u,i,a) <-
c              K3     +   [  (b,i |a,be)             ] . t(b ,u)
c
c       Gvv(be,a)    <-
c       Gvv1  +           [ Hvv(a,be) ]
c       Gvv2  - sum(i)    [ Fvo(a,i) . t1(be,i) ]
c       Gvv3  + sum(b,i)  [ -(b,be|a,i) +2(b,i|a,be) ] . t(b,i )
c
c       Goo(i,u) <-
c       Goo1  +  <-       [ Hoo(i,u) ]
c       Goo2  + sum(a)    [ Fvo(a,i) . t1(a,u) ]
c       Goo3  + sum(a,j)  [ (2(aj|iu) - (ai|ju)) . t1(a,j) ]
c
c       terms Goo2,3 implemented as:
c       Goo2  + sum(be)   [ Fvo(be,i) . t1(be,u) ]
c       Goo3  + sum(be,j) [ (2(be,j|iu) - (be,i|ju)) . t1(be,j) ]
c
c       T2(be,ga,u,v) <<-
c       T26  + sum(a)     [ (be,u|ga,a) . t1(a,v)         ]
c       T27  - sum(a,i)   [ (iu|ga,a) . t1(be,i) . t1(a,v)]
c       T28  - sum(i)     [ (be,u|iv)          . t1(ga,i) ]
c       T29  - sum(a,i)   [ (be,u|ia). t1(a,v) . t1(ga,i) ]
c
c       N.B. T26-T29 implementede in expression:
c       T2(a,be,u,v)  <-
c       T26  + sum(b)     [           (u,a|be,b) . t1(b,v)  ]
c       T27  - sum(b,j)   [ t1(a,j) . (u,j|be,b) . t1(b,v)  ]
c       T28  - sum(j)     [           (u,a|j,v)  . t1(be,j) ]
c       T29  - sum(b,j)   [ t1(b,v) . (u,a|j,b)  . t1(be,j) ]
c        T2(a,be,u,v) is stored as X(a,u,be,v) in Tmp2Name(a,be)
c
c
c        N.B. Process V3(j,v,u,a') <- V1(m,j,v) . L1(T)(m,u,a') run
c        for Na*Nbe times, ie. Nbe . o3mv competee with Nbe . I+O (o3v)
c
c
c       Extra cost - multiple reading:
c                    1) N_be . T;
c                    2) N_ga . Q
c                    3) N_ga . K
c
c
c1      intermediates Q and K will be temporarry stored in files
c       Tmp1Name(be,a)
c       do beGrp
c       do aGrp
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
        integer NvGrp,LunAux

c       help variables
        integer dim1,dim2,dima,dimb,dimbe
        integer aGrp,bGrp,beGrp,adda,addb,addbe
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
c
        integer PossT,maxdim

c
        character*6 LunName
c
c       --- introduction part ---
c
cx      Distribute memory
c       @@ vyladit distmem, mozno netreba vsetko, toto je len prefackany
c        DistMemo3v3
c
        PossT=PossFree
        if (printkey.ge.10) then
        write (6,*) 'PossFree',PossFree,NvGrp,maxdim
        end if
        call DistMemo3v3chol (NvGrp,maxdim,
     c       PossV1,PossV2,PossV3,PossV4,
     c       PossH1,PossH2,PossH3,PossH4,
     c       PossM1,PossM2,PossM3,PossM4,PossM5,
     c       PossK,PossQ,
     c       PossT)
c
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
cmp!           write (6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
       write (6,*) ' Not Enough memory in o3v3chol step!',
     & 'Increase large and/or small segmentation ',
     c                  (1.0d0*PossT)/(1.0d0*wrksize)
           call abend()
        end if
c
c
c@@
c        call Chck_mkJ
c        call Chck_mkK
c@@
c
c       read V1(m,ij) <- L1(m,ij)
        LunName=L0Name
        dim1=nc*no*(no+1)/2
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c        Expand M4(m,i,j) <- V1(m,ij)
        dim1=no*(no+1)/2
        call Exp1 (wrk(PossV1),wrk(PossM4),nc,dim1,no)
c
cG        Vanish Gvv,Goo
        dim1=nv*nv
        call mv0zero (dim1,dim1,wrk(PossGvv))
        dim1=no*no
        call mv0zero (dim1,dim1,wrk(PossGoo))
c
c
c##G     Vanish Xyes
c
        do dim1=1,NvGrp
        do dim2=1,NvGrp
          Xyes(dim1,dim2)=0
        end do
        end do
c
c
        addbe=0
c
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
c
c##     test, if this beGrp is planed to be run on this node
        dim1=0
        do dim2=1,NvGrp
          dim1=dim1+BeAID(myRank,beGrp,dim2)
        end do
        dim1=dim1+BetaID(myRank,beGrp)
        if (dim1.eq.0) goto 11
        if (printkey.gt.1) then
        write (6,*) ' Chol PID, beGrp',myRank,beGrp
cmp
        Call CWTime(TCpu,TWall)
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
c          Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
c
c         vanish M1(m,i,u),M2(m,be',u),M3(m,be',u)
          dim1=nc*no*no
          call mv0zero (dim1,dim1,wrk(PossM1))
          dim1=nc*no*dimbe
          call mv0zero (dim1,dim1,wrk(PossM2))
          call mv0zero (dim1,dim1,wrk(PossM3))
c
cT1G          vanish H4(be',u)
          dim1=dimbe*no
          call mv0zero (dim1,dim1,wrk(PossH4))
c
c
c
           if (BetaID(myRank,beGrp).eq.1) then
c##       contributions Goo2,3  must be taken only once per Beta'
c
cGoo2.1            Extract H2(be',i) <- Fvo(be,i) (vlastnu rutinu netreba)
            call ExtT1 (wrk(PossH2),wrk(PossFvo),dimbe,addbe)
c
cGoo2.2f    Goo(i,u) <<- H2(T)(be',i) . H1(be',u)
            call mc0c1at3b (dimbe,no,dimbe,no,no,no,
     c                    no,dimbe,no,
     c                    wrk(PossH2),wrk(PossH1),wrk(PossGoo))
c
c
cGoo3.1            read V1(be',I,JK) <- I1(be',I,JK)
            dim1=dimbe*no*no*(no+1)/2
            LunName=I1Name(beGrp)
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
cGoo3.2            Make V2(be',j,i,u) <- 2 V1(be',j|iu) - V1(be',i|ju)
            call MkV_Goo3 (wrk(PossV1),wrk(PossV2),dimbe,no)
c
cGoo3.3f    Goo(i,u) <<- V2(T)(be',j,i,u) . H1(be,j)
            dim1=dimbe*no
c@x            dorobit rutinu  mv0v1at3u
c           call mv0v1at3u (dim1,no*no,dim1,no*no,
c    c                      no*no,dim1,1,1,
c    c                      wrk(PossV2),wrk(PossH1),wrk(PossGoo))
c
c            zatial odflaknute takto
             call mc0c1at3b (dim1,no*no,dim1,1,no*no,1,
     c                     no*no,dim1,1,
     c                     wrk(PossV2),wrk(PossH1),wrk(PossGoo))
c    c
c@x
           end if
c
c
          addb=0
c
          do bGrp=1,NvGrp
          dimb=DimGrpv(bGrp)
c
cG            Extract H2(b',I) <- T1(b,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dimb,addb)
c
cQK3.1      read V3(m,i,b') <- L1(m,i,b')
            LunName=L1Name(bGrp)
            dim1=nc*no*dimb
            call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
c
cQK3.3      M1(m,i,u)   <<-  V3(m,i  ,b') . H2(b',u)
            call mc0c1a3b(nc*no,dimb,dimb,no,nc*no,no,
     c                    nc*no,dimb,no,
     c                    wrk(PossV3),wrk(PossH2),wrk(PossM1))
c
cQ3.4       read V1(m,be',b')<- [V2(m,beb') or V2(m,b'be)]<- L2(m,be',b')
c            in cases be.le.b through interstep V2
            if (beGrp.gt.bGrp) then
              LunName=L2Name(beGrp,bGrp)
              dim1=nc*dimb*dimbe
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
            else if (beGrp.eq.bGrp) then
              LunName=L2Name(beGrp,bGrp)
              dim1=nc*dimb*(dimbe+1)/2
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
              dim1=dimb*(dimbe+1)/2
              call Exp1 (wrk(PossV2),wrk(PossV1),nc,dim1,dimb)
            else
              LunName=L2Name(bGrp,beGrp)
              dim1=nc*dimb*dimbe
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
              call Map3_132 (wrk(PossV2),wrk(PossV1),nc,dimb,dimbe)
            end if
c
cQ3.5       M2(m,be',u) <<-  V1(m,be',b') . H2(b',u)
            call mc0c1a3b(nc*dimbe,dimb,dimb,no,nc*dimbe,no,
     c                    nc*dimbe,dimb,no,
     c                    wrk(PossV1),wrk(PossH2),wrk(PossM2))
c
c
c           read V2(be'b',I,J) <- T2(be',b,I,J)
            LunName=T2Name(beGrp,bGrp)
            if (bGrp.eq.beGrp) then
              dim1=dimbe*(dimbe+1)*no*no/2
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c             Expand V2(be',b',I,J) <- V1(be'b',I,J)
              dim1=dimb*(dimb+1)/2
              call ExpT2 (wrk(possV1),wrk(PossV2),dimbe,dim1,no)
            else
              dim1=dimbe*dimb*no*no
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
            end if
c
c           Make Tau(be',b',I,J) in V2(be',b',I,J)
               call MkTau_chcc (wrk(PossV2),wrk(PossH1),wrk(PossH2),
     c                       dimbe,dimb,no,1.0d0,1.0d0)

c           Make V1(i,b',be',u) <- 2 V2(be',b',u,i) - V2(be',b',i,u)
            call MkT_T17 (wrk(PossV1),wrk(PossV2),dimb,dimbe,no)
c
c           M3(m,be',u) <<- V3(m,i,b') . V1(i,b',be',u)
            dim1=dimb*no
            dim2=dimbe*no
            call mc0c1a3b (nc,dim1,dim1,dim2,nc,dim2,
     c                     nc,dim1,dim2,
     c                     wrk(PossV3),wrk(PossV1),wrk(PossM3))
c
c
          addb=addb+dimb
          end do
c
c
          if (BetaID(myRank,beGrp).eq.1) then
c##       contribution T16  must be taken only once per Beta'
c
c           Map V1(be',m,i) <- M2(m,be',i)
            call Map3_213 (wrk(PossM2),wrk(PossV1),nc,dimbe,no)
c
cT162.f     H4(be',u) <<- - V1(be',m,i) . M4(m,i,u)
            call mc0c2a3b (dimbe,nc*no,nc*no,no,dimbe,no,
     c                   dimbe,nc*no,no,
     c                   wrk(PossV1),wrk(PossM4),wrk(PossH4))
c
          end if
c
c
cT27.1    V4(be',v,u,j ) <- M2(T)(m,be',v) . M4(m,u,j)
          dim1=dimbe*no*no*no
          call mv0zero (dim1,dim1,wrk(PossV4))
          call mc0c1at3b (nc,dimbe*no,nc,no*no,dimbe*no,no*no,
     c                    dimbe*no,nc,no*no,
     c                    wrk(PossM2),wrk(PossM4),wrk(PossV4))
c
c
          adda=0
c
          do aGrp=1,NvGrp
          dima=DimGrpv(aGrp)
c
c         Test, if this Be' A' combination is to be run on this node
          if (BeAID(myRank,beGrp,aGrp).eq.0) goto 12
          if (printkey.ge.10) then
          write (6,*) ' Chol PID, be,a',myRank,beGrp,aGrp
          end if
c
c            Extract H2(a',I) <- T1(a,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dima,adda)
c
c           read Q(be',u,i,a'),K(be',u,i,a')
            LunName=Tmp1Name(beGrp,aGrp)
            dim1=dimbe*dima*no*no
            call GetX (wrk(PossQ),dim1,LunAux,LunName,1,0)
            call GetX (wrk(PossK),dim1,LunAux,LunName,0,1)
c
c           read V1(m,be',a')<- [V2(m,bea') or V2(m,a',be')] <- L2(m,be',a')
c            in cases be.le.b through interstep V2
            if (beGrp.gt.aGrp) then
              LunName=L2Name(beGrp,aGrp)
              dim1=nc*dima*dimbe
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
            else if (beGrp.eq.aGrp) then
              LunName=L2Name(beGrp,aGrp)
              dim1=nc*dima*(dimbe+1)/2
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
              dim1=dima*(dimbe+1)/2
              call Exp1 (wrk(PossV2),wrk(PossV1),nc,dim1,dimbe)
            else
              LunName=L2Name(aGrp,beGrp)
              dim1=nc*dima*dimbe
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
              call Map3_132 (wrk(PossV2),wrk(PossV1),nc,dima,dimbe)
            end if
c
c           V3(a',u) <- V1(T)(m,be',a')(L2) . M3(m,be',u)
            call mv0zero (dima*no,dima*no,wrk(PossV3))
            call mc0c1at3b (nc*dimbe,dima,nc*dimbe,no,dima,no,
     c                      dima,nc*dimbe,no,
     c                      wrk(PossV1),wrk(PossM3),wrk(PossV3))
c
cT17.f      Add T1n(a',u) <<- V3(a',u)
                 call AdT_T17 (wrk(PossT1n),wrk(PossV3),
     c                    nv,dima,no,adda,1.0d0)
c
c           V3(be',a',i,u) <-  V1(m,be',a')(L2) . M1(m,i,u)
            dim1=dima*dimbe*no*no
            call mv0zero (dim1,dim1,wrk(PossV3))
            call mc0c1at3b (nc,dimbe*dima,nc,no*no,dimbe*dima,no*no,
     c                      dima*dimbe,nc,no*no,
     c                      wrk(PossV1),wrk(PossM1),wrk(PossV3))
c
c           map V2(be',u,i,a') <-  V3(be',a',i,u)
            call Map4_1432 (wrk(PossV3),wrk(PossV2),dimbe,dima,no,no)
c
c           Q(be',u,i,a')  <<- -V2(be',u,i,a')
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossQ),1,-1.0d0)
c
c           Gvv(be',a') <<- 2 sum(i) [V2(be',i,i,a')
            call AdV_G2 (wrk(PossGvv),wrk(PossV2),
     c                   nv,dimbe,dima,no,addbe,adda,2.0d0)
c
c           K(be',u,i,a')  <<-  V2(be',u,i,a')
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossK),1,1.0d0)
c
c           read M5(m,i,a') <- L1(m,i,a')
            LunName=L1Name(aGrp)
            dim1=dima*nc*no
            call GetX (wrk(PossM5),dim1,LunAux,LunName,1,1)
c
c           V2(be',u,i,a') <-  M2(m,be',u ) . M5(m,i,a')(L1)
            dim1=dima*dimbe*no*no
            call mv0zero (dim1,dim1,wrk(PossV2))
            call mc0c1at3b (nc,dimbe*no,nc,no*dima,dimbe*no,dima*no,
     c                      dimbe*no,nc,dima*no,
     c                      wrk(PossM2),wrk(PossM5),wrk(PossV2))
c
c           Q(be',u,i,a') <<- 2 V2(be',u,i,a')
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossQ),1,2.0d0)
c
cx          write Q and K submatrix to corresponding files
            LunName=Tmp1Name(beGrp,aGrp)
            dim1=dimbe*dima*no*no
            call SaveX (wrk(PossQ),dim1,LunAux,LunName,1,0)
            call SaveX (wrk(PossK),dim1,LunAux,LunName,0,1)
c@@
c        call Chck_K (wrk(PossK),dimbe,addbe,dima,adda)
c        call Chck_Q (wrk(PossQ),dimbe,addbe,dima,adda)
c@@
c
c           Gvv(be',a') <<- - sum(i) [V2(be,i,i,a')
            call AdV_G2 (wrk(PossGvv),wrk(PossV2),
     c                   nv,dimbe,dima,no,addbe,adda,-1.0d0)
c
c
c            Q array (alredy relesed) will be served for X
c
cT26.1f     Map Q(a',U,be',V) <-  2 . V2(be',V,U,a')
c           ie. here a' is a firs index, that means X (Q) is saved
c           in different order as K,Q
            call Map4_3421 (wrk(PossV2),wrk(PossQ),dimbe,no,no,dima)
            dim1=dima*no*no*dimbe
            call mv0sv (dim1,dim1,wrk(PossQ),2.0d0)
c
cT27.2      Map H3(j,a') <- H2(a',j)(T1)
            call Map2_21 (wrk(PossH2),wrk(PossH3),dima,no)
c
cT27.3f            V2(be',v,u,a') <- - V4(be',v,u,j) . H3(j,a')
            dim1=dimbe*no*no*dima
            call mv0zero (dim1,dim1,wrk(PossV2))
            dim1=dimbe*no*no
                 call mc0c2a3b (dim1,no,no,dima,dim1,dima,
     c                     dim1,no,dima,
     c                     wrk(PossV4),wrk(PossH3),wrk(PossV2))
c
cT289.1            V1(m,j,v) <- M4(m,j,v) (L1)
            dim1=nc*no*no
            call mv0u (dim1,wrk(PossM4),1,wrk(PossV1),1)
c
cT289.2     V1(m,j,v) <<- M1(m,j,v)
            dim1=nc*no*no
            call mv0v1u (dim1,wrk(PossM1),1,wrk(PossV1),1,1.0d0)
c
cT289.3     V3(j,v,u,a') <- V1(T) (m,j,v) . M5(m,u,a')(L1)
            dim1=no*no*no*dima
            call mv0zero (dim1,dim1,wrk(PossV3))
               call mc0c1at3b (nc,no*no,nc,no*dima,no*no,no*dima,
     c                      no*no,nc,no*dima,
     c                      wrk(PossV1),wrk(PossM5),wrk(PossV3))
c
cT289.4f    V2(be',v,u,a') <<- - H1(be',j)(T1) . V3(j,v,u,a')
            dim1=no*no*dima
                 call mc0c2a3b (dimbe,no,no,dim1,dimbe,dim1,
     c                     dimbe,no,dim1,
     c                     wrk(PossH1),wrk(PossV3),wrk(PossV2))
c
cT2G        Map V1(a',U,be',V) <- V2(be',V,U,a')
            call Map4_3421 (wrk(PossV2),wrk(PossV1),dimbe,no,no,dima)
c
cT2G        Add Q(a',u,be',v)(X) <<- V1(a',U,be',V) (Faktor 2 koli X)
            dim1=dimbe*dima*no*no
            call mv0v1u (dim1,wrk(PossV1),1,wrk(PossQ),1,2.0d0)
c
cT2G        write X(a',u,be',v) = Q(a',u,be',v)
            LunName=Tmp2Name(aGrp,beGrp)
            dim1=dimbe*dima*no*no
            Xyes(aGrp,beGrp)=1
            call SaveX (wrk(PossQ),dim1,LunAux,LunName,1,1)
c
12        adda=adda+dima
          end do
c
cT1G        Add T1n(be,u) <<- H4(be',u)
        call AdT_T17 (wrk(PossT1n),wrk(PossH4),nv,dimbe,no,addbe,1.0d0)
c
11      addbe=addbe+dimbe
        end do
c
c
#ifdef _MOLCAS_MPP_
c##        Synchronizacny bod:
c        Allreduce Goo,Gvv
        dim1=no*no
        call gadgop (wrk(PossGoo),dim1,'+')
        dim1=nv*nv
        call gadgop (wrk(PossGvv),dim1,'+')
#endif
c
c
cGoo1.1f Add Goo(i,u) <<- Hoo(i,u)
        dim1=no*no
        call mv0v1u (dim1,wrk(PossHoo),1,wrk(PossGoo),1,1.0d0)
c
c
cGvv1.1 Map V1(be,a) <- Hvv(a,be)
        call Map2_21 (wrk(PossHvv),wrk(PossV1),nv,nv)
c
cGvv1.2f Add Gvv(be,a) <<- V1(be,a)
        dim1=nv*nv
        call mv0v1u (dim1,wrk(PossV1),1,wrk(PossGvv),1,1.0d0)
c
c
cGvv2.1 Map V2(i,a) <- Fvo(a,i)
        call Map2_21 (wrk(PossFvo),wrk(PossV2),nv,no)
c
cGvv2.2 Calc V1(be,a) <- t1(be,i) . V2(i,a)
        dim1=nv*nv
        call mv0zero (dim1,dim1,wrk(PossV1))
        call mc0c1a3b (nv,no,no,nv,nv,nv,
     c                 nv,no,nv,
     c                 wrk(PossT1o),wrk(PossV2),wrk(PossV1))
c
cGvv2.3f Add Gvv(be,a) <<- - V1(be,a)
        dim1=nv*nv
        call mv0v1u (dim1,wrk(PossV1),1,wrk(PossGvv),1,-1.0d0)
c
c@@
c        call Chck_Goo (wrk(PossGoo))
c        call Chck_Gvv (wrk(PossGvv))
c@@
c
        return
        end
c
c        --------------------
c
        subroutine Chck_Goo (Goo)
c
c        check Goo (i,u)
c
        implicit none
#include "chcc1.fh"
        real*8 Goo(1:no,1:no)
c
c        help var
        integer i,u,j,a,bad
c       integer b
        real*8 s
c
        bad=0
c
        do i=1,no
        do u=1,no
c
c          s=0.0d0
c          do j=1,no
c          do a=1,nv
c          do b=1,nv
c          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(a,b,u,j)
c          end do
c          end do
c          end do
c
           s=Hooc(i,u)
c
          do j=1,no
          do a=1,nv
           s=s+(2.0d0*Q1(a,j,i,u)-Q1(a,i,j,u))*T1c(a,j)
          end do
          end do
c
          Gooc(i,u)=s
c
          if (abs(Goo(i,u)-s).gt.1.0d-10) then
          bad=bad+1
          end if
c
        end do
        end do
c
        write (6,*) ' Goo Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Gvv (Gvv)
c
c        check Gvv(be,a)
c
        implicit none
#include "chcc1.fh"
        real*8 Gvv(1:nv,1:nv)
c
c        help var
        integer i,a,b,be,bad
c       integer j
        real*8 s
c
        bad=0
c
        do be=1,nv
        do a=1,nv
c
c          s=0.0d0
c          do i=1,no
c          do j=1,no
c          do b=1,nv
c          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(be,b,i,j)
c          end do
c          end do
c          end do
c          s=-s
c
          s=Hvvc(be,a)
c
          do i=1,no
          do b=1,nv
          s=s+(2.0d0*Q3(a,be,b,i)-Q3(b,be,a,i))*T1c(b,i)
          end do
          end do
c
          Gvvc(be,a)=s
c
          if (abs(Gvv(be,a)-s).gt.1.0d-10) then
          bad=bad+1
          end if
c
        end do
        end do
c
        write (6,*) ' Gvv Chck :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_0 (dim,A)
c
c        check zero
c
        implicit none
        integer dim
        real*8 A(1:dim)
c
c        help var
        integer i,bad
c
        bad=0
        do i=1,dim
          if (abs(A(i)).gt.1.0d-10) then
          bad=bad+1
          end if
        end do
c
        write (6,*) ' Nonzero elements ',bad,dim
c
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_K (K,dimbe,addbe,dima,adda)
c
c        check K(be,u,i,a)
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dima,adda
        real*8 K(1:dimbe,1:no,1:no,1:dima)
c
        integer be,u,i,a,bad
        real*8 s
c
        bad=0
        do a=adda+1,adda+dima
        do i=1,no
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          s=Kc(i,be,u,a)
c
          if (abs(K(be-addbe,u,i,a-adda)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck K :',bad
c99        format (a9,1x,i8,1x,4(i3,1x))
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Q (Q,dimbe,addbe,dima,adda)
c
c        check Q(be,u,i,a)
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dima,adda
        real*8 Q(1:dimbe,1:no,1:no,1:dima)
c
        integer be,u,i,a,bad
        real*8 s,sk,sj
c
        bad=0
c
        do a=adda+1,adda+dima
        do i=1,no
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          sj=Jc(be,i,u,a)
          sk=Kc(i,be,u,a)
c
          s=2*sj-sk
c
          if (abs(Q(be-addbe,u,i,a-adda)-s).gt.1.0d-10) then
            bad=bad+1
c           Q(be-addbe,u,i,a-adda)=s
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck Q :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_AA (A)
c
c        check T(a,b,i,j)
c
        implicit none
#include "chcc1.fh"
c        real*8 T(1:nv,1:nv,1:no,1:no)
         real*8 A(1:no*(no+1)/2,no,no)
c
        integer j,i,ij,u,v,bad
        real*8 s
c
        bad=0
        do v=1,no
        do u=1,no
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
c
           s=Ac(i,j,u,v)
c
          if (abs(A(ij,u,v)-s).gt.1.0d-10) then
            bad=bad+1
c    A(ij,u,v)=s
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck AA :',bad
c
        return
        end
c
c        check Y(be,u,ga,v)
c
c        ---------------------
c
        subroutine Chck_Tx (T)
c
c        check T(a,b,i,j)
c
        implicit none
#include "chcc1.fh"
c        real*8 T(1:nv,1:nv,1:no,1:no)
         real*8 T(1:nv,1:no,1:nv,1:no)
c        real*8 T(1:nv,1:nv,1:no)
c
        integer b,j,a,i,bad
        real*8 s
c
        bad=0
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
c
           s=T2c(a,b,i,j)
c
          if (abs(T(b,i,a,j)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck T2 :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Vx (V)
c
c        check V
c
        implicit none
#include "chcc1.fh"
         real*8 V(1:nv,1:no,1:nv,1:no)
c
        integer be,u,i,a,bad
        real*8 s
c
        bad=0
        do a=1,nv
        do i=1,no
        do u=1,no
        do be=1,nv
c
           s=2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a)
           s=Kc(i,be,u,a)
c
          if (abs(V(be,u,a,i)-s).gt.1.0d-10) then
            bad=bad+1
            V(be,u,a,i)=s
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck Vx :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Y (Y,dimbe,addbe,dimga,addga)
c
c        check Y(be,u,ga,v)
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dimga,addga
        real*8 Y(1:dimbe,1:no,1:dimga,1:no)
c
        integer be,u,ga,v,bad
        integer a,i
        real*8 s
c
        bad=0
        do v=1,no
        do ga=addga+1,addga+dimga
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          s=0.0d0
          do i=1,no
          do a=1,nv
           s=s+Kc(i,be,u,a)*T2c(ga,a,i,v)
          end do
          end do
c
          if (abs(Y(be-addbe,u,ga-addga,v)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck Y :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_X (X,dimbe,addbe,dimga,addga)
c
c        check X(be,u,ga,v)
c        do X clenu nieje zahrnuty prispevok od Aex, preto picuje
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dimga,addga
         real*8 X(1:dimbe,1:no,1:dimga,1:no)
c
        integer be,u,ga,v,bad
        integer a,i,j
        real*8 s,s1
c
        bad=0
        do v=1,no
        do ga=addga+1,addga+dimga
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          s=0.0d0
c
c         X2 (T24)   + Gvv(a,be)   . t2(ga,a,v,u)
          s1=0.0d0
          do a=1,nv
           s1=s1+Gvvc(be,a)*T2c(a,ga,u,v)
          end do
          s=s+2.0d0*s1
c
c         X3 (T25)   - Goo(i,u)    . t2(ga,be,v,i)
          s1=0.0d0
          do i=1,no
           s1=s1+Gooc(i,u)*T2c(ga,be,v,i)
          end do
          s=s-2.0d0*s1
c
c
c          T1 cleny
          s1=0.0d0
c
          do a=1,nv
           s1=s1+Q3(a,ga,be,u)*T1c(a,v)
          end do
c
          do i=1,no
           s1=s1-Q1(be,u,i,v)*T1c(ga,i)
          end do
c
          do i=1,no
          do a=1,nv
           s1=s1-(Q22(a,ga,i,u)*T1c(be,i))*T1c(a,v)
           s1=s1-(Q21(a,i,be,u)*T1c(a,v))*T1c(ga,i)
          end do
          end do
c
          s=s+2.0d0*s1
c
c
c         X4 (T22)   + sum(i,j)   [ Ta(be,ga,i,j) . A(i,j,u,v)  ]
          s1=0.0d0
          do i=1,no
           do j=1,no
          s1=s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
          end do
           s=s+s1
c
c
c         X1   <-  Q(be,u,i,a) . (2 t2(ga,a,v,i) - t2(ga,a,i,v))
c         calc as (2J(be,i,u,a)-K(i,be,u,a)*(2t2(a,ga,i,v)-t2(ga,a,i,v)
          s1=0.d0
          do i=1,no
          do a=1,nv
          s1=s1+(2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a))*
     c          (2.0d0*T2c(a,ga,i,v)-T2c(ga,a,i,v))
          end do
          end do
          s=s+s1
c
c         <- V1(be',u,ga',v)
          s1=Q21(be,u,ga,v)
          s=s+s1

           if (abs(X(be-addbe,u,ga-addga,v)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck X :',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Xred (X,dimbe,addbe,dimga,addga)
c
c        check X(be,u,ga,v)
c        do X clenu nieje zahrnuty prispevok od Aex, preto picuje
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dimga,addga
         real*8 X(1:dimbe,1:no,1:dimga,1:no)
c
        integer be,u,ga,v,bad
        integer a,i,j
        real*8 s,s1
c
        bad=0
        do v=1,no
        do ga=addga+1,addga+dimga
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          s=0.0d0
c
c         X2 (T24)   + Gvv(a,be)   . t2(ga,a,v,u)
          s1=0.0d0
          do a=1,nv
           s1=s1+Gvvc(be,a)*T2c(a,ga,u,v)
          end do
          s=s+2.0d0*s1
c
c         X3 (T25)   - Goo(i,u)    . t2(ga,be,v,i)
          s1=0.0d0
          do i=1,no
           s1=s1+Gooc(i,u)*T2c(ga,be,v,i)
          end do
          s=s-2.0d0*s1
c
c
c          T1 cleny
          s1=0.0d0
c
          do a=1,nv
           s1=s1+Q3(a,ga,be,u)*T1c(a,v)
          end do
c
          do i=1,no
           s1=s1-Q1(be,u,i,v)*T1c(ga,i)
          end do
c
          do i=1,no
          do a=1,nv
           s1=s1-(Q22(a,ga,i,u)*T1c(be,i))*T1c(a,v)
           s1=s1-(Q21(a,i,be,u)*T1c(a,v))*T1c(ga,i)
          end do
          end do
c
cred          s=s+2.0d0*s1
c
c
c         X4 (T22)   + sum(i,j)   [ Ta(be,ga,i,j) . A(i,j,u,v)  ]
          s1=0.0d0
          do i=1,no
           do j=1,no
          s1=s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
          end do
           s=s+s1
c
c
c         X1   <-  Q(be,u,i,a) . (2 t2(ga,a,v,i) - t2(ga,a,i,v))
c         calc as (2J(be,i,u,a)-K(i,be,u,a)*(2t2(a,ga,i,v)-t2(ga,a,i,v)
          s1=0.d0
          do i=1,no
          do a=1,nv
          s1=s1+(2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a))*
     c          (2.0d0*T2c(a,ga,i,v)-T2c(ga,a,i,v))
          end do
          end do
          s=s+s1
c
c         <- V1(be',u,ga',v)
          s1=Q21(be,u,ga,v)
          s=s+s1

           if (abs(X(be-addbe,u,ga-addga,v)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck X :',bad
c
        return
        end








c
c        ---------------------
c
        subroutine Chck_mkK
c
c        make K(be,u,i,a)
c
        implicit none
#include "chcc1.fh"
c
        integer be,u,i,a
        integer j,b
        real*8 s
c
        do a=1,nv
        do i=1,no
        do u=1,no
        do be=1,nv
c
          s=0.0d0
c
          s=s+Q22(be,a,i,u)
c
          do j=1,no
          s=s-Q1(a,j,i,u)*T1c(be,j)
          end do
c
          do b=1,nv
          s=s+Q3(a,be,b,i)*T1c(b,u)
          end do
c
          do j=1,no
          do b=1,nv
           s=s-Q21(b,i,a,j)*(T2c(b,be,u,j)/2+T1c(b,u)*T1c(be,j))
          end do
          end do
c
          Kc(i,be,u,a)=s
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' K done '
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_mkJ
c
c        make J(be,u,i,a)
c
        implicit none
#include "chcc1.fh"
c
        integer be,u,i,a
        integer j,b
        real*8 sj
c
        do a=1,nv
        do i=1,no
        do u=1,no
        do be=1,nv
c
          sj=0.0d0
c
          do j=1,no
           sj=sj-Q1(a,i,j,u)*T1c(be,j)
          end do
c
          do b=1,nv
           sj=sj+Q3(b,be,a,i)*T1c(b,u)
          end do
c
          do j=1,no
          do b=1,nv
           sj=sj+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(be,b,u,j)/2
           sj=sj-Q21(a,i,b,j)*(T2c(b,be,u,j)/2+T1c(b,u)*T1c(be,j))
          end do
          end do
c
          sj=sj+Q21(be,u,a,i)
c
          Jc(be,i,u,a)=sj
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' J done'
c
        return
        end
