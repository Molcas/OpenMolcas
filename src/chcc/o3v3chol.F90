!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine o3v3chol (wrk,wrksize,NvGrp,maxdim,LunAux)
!
!       This routine do:
!
!       --- Part II - generation of Q3,K3,Gvv,Goo,T26-9 in X form
!
!       T(be,u)      <<-
!       T162 - sum(b,i)   [ (iu|be,b) . t(b,i) ]
!       T17  + sum(a,b,i) [ 2(a,i|b,be) - (b,i|a,be) ] . Ta(a,b,i,u)
!          = + sum(a,b,i)  (a,i|b,be) . [ 2 Ta(a,b,i,u) - Ta(a,b,u,i) ]
!
!       N.B. T17 implementede in expression:
!       T(a,u)      <<-
!       T17  + sum(be,b,i)  (b,i|be,a) . [ 2 Ta(b,be,i,u) - Ta(b,be,u,i) ]
!
!
!       Q(be,u,i,a) <-
!              Q3     +   [ 2(b,be|a,i) - (b,i|a,be) ] . t(b,u )
!
!       K(be,u,i,a) <-
!              K3     +   [  (b,i |a,be)             ] . t(b ,u)
!
!       Gvv(be,a)    <-
!       Gvv1  +           [ Hvv(a,be) ]
!       Gvv2  - sum(i)    [ Fvo(a,i) . t1(be,i) ]
!       Gvv3  + sum(b,i)  [ -(b,be|a,i) +2(b,i|a,be) ] . t(b,i )
!
!       Goo(i,u) <-
!       Goo1  +  <-       [ Hoo(i,u) ]
!       Goo2  + sum(a)    [ Fvo(a,i) . t1(a,u) ]
!       Goo3  + sum(a,j)  [ (2(aj|iu) - (ai|ju)) . t1(a,j) ]
!
!       terms Goo2,3 implemented as:
!       Goo2  + sum(be)   [ Fvo(be,i) . t1(be,u) ]
!       Goo3  + sum(be,j) [ (2(be,j|iu) - (be,i|ju)) . t1(be,j) ]
!
!       T2(be,ga,u,v) <<-
!       T26  + sum(a)     [ (be,u|ga,a) . t1(a,v)         ]
!       T27  - sum(a,i)   [ (iu|ga,a) . t1(be,i) . t1(a,v)]
!       T28  - sum(i)     [ (be,u|iv)          . t1(ga,i) ]
!       T29  - sum(a,i)   [ (be,u|ia). t1(a,v) . t1(ga,i) ]
!
!       N.B. T26-T29 implementede in expression:
!       T2(a,be,u,v)  <-
!       T26  + sum(b)     [           (u,a|be,b) . t1(b,v)  ]
!       T27  - sum(b,j)   [ t1(a,j) . (u,j|be,b) . t1(b,v)  ]
!       T28  - sum(j)     [           (u,a|j,v)  . t1(be,j) ]
!       T29  - sum(b,j)   [ t1(b,v) . (u,a|j,b)  . t1(be,j) ]
!        T2(a,be,u,v) is stored as X(a,u,be,v) in Tmp2Name(a,be)
!
!
!        N.B. Process V3(j,v,u,a') <- V1(m,j,v) . L1(T)(m,u,a') run
!        for Na*Nbe times, ie. Nbe . o3mv competee with Nbe . I+O (o3v)
!
!
!       Extra cost - multiple reading:
!                    1) N_be . T;
!                    2) N_ga . Q
!                    3) N_ga . K
!
!
!1      intermediates Q and K will be temporarry stored in files
!       Tmp1Name(be,a)
!       do beGrp
!       do aGrp
!         Q(K) (be',u,i,a')
!       end do
!       end do
!
!2      Structure of files, where selected group of (pq|rs) are
!       stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
!
!       (A'I|JK)  I1inxx xx - Group of A'
!
!       (A'I|B'J) I2xxyy xx - Group of A'
!                        yy - Group of B'
!
!       (A'B'|IJ) I3xxyy xx - Group of A'
!                        yy - Group of B'
!
!
!3      Structure of Cholesky vector files
!
!       L1(m,I ,A')  L1vcxx xx - Group of A'
!
!       L2(m,A',B')  L2xxyy xx - Group of A'
!                           yy - Group of B'
!
!4      Structure of Amplitude file
!       t2(A',B',IJ)  T2xxyy xx - Group of A'
!                            yy - Group of B'
!
!
        use Para_Info, only: MyRank
        implicit none
#include "chcc1.fh"
#include "parcc.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
#include "wrk.fh"
#include "chcc_casy.fh"
!
        integer NvGrp,LunAux

!       help variables
        integer dim1,dim2,dima,dimb,dimbe
        integer aGrp,bGrp,beGrp,adda,addb,addbe
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
!
        integer PossT,maxdim

!
        character*6 LunName
!
!       --- introduction part ---
!
!x      Distribute memory
!       @@ vyladit distmem, mozno netreba vsetko, toto je len prefackany
!        DistMemo3v3
!
        PossT=PossFree
        if (printkey.ge.10) then
        write (6,*) 'PossFree',PossFree,NvGrp,maxdim
        end if
        call DistMemo3v3chol (NvGrp,maxdim,                             &
     &       PossV1,PossV2,PossV3,PossV4,                               &
     &       PossH1,PossH2,PossH3,PossH4,                               &
     &       PossM1,PossM2,PossM3,PossM4,PossM5,                        &
     &       PossK,PossQ,                                               &
     &       PossT)
!
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
!mp!           write (6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
       write (6,*) ' Not Enough memory in o3v3chol step!',              &
     & 'Increase large and/or small segmentation ',                     &
     &                  (1.0d0*PossT)/(1.0d0*wrksize)
           call abend()
        end if
!
!
!@@
!        call Chck_mkJ
!        call Chck_mkK
!@@
!
!       read V1(m,ij) <- L1(m,ij)
        LunName=L0Name
        dim1=nc*no*(no+1)/2
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!        Expand M4(m,i,j) <- V1(m,ij)
        dim1=no*(no+1)/2
        call Exp1 (wrk(PossV1),wrk(PossM4),nc,dim1,no)
!
!G        Vanish Gvv,Goo
        dim1=nv*nv
        call mv0zero (dim1,dim1,wrk(PossGvv))
        dim1=no*no
        call mv0zero (dim1,dim1,wrk(PossGoo))
!
!
!##G     Vanish Xyes
!
        do dim1=1,NvGrp
        do dim2=1,NvGrp
          Xyes(dim1,dim2)=0
        end do
        end do
!
!
        addbe=0
!
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
!
!##     test, if this beGrp is planed to be run on this node
        dim1=0
        do dim2=1,NvGrp
          dim1=dim1+BeAID(myRank,beGrp,dim2)
        end do
        dim1=dim1+BetaID(myRank,beGrp)
        if (dim1.eq.0) goto 11
        if (printkey.gt.1) then
        write (6,*) ' Chol PID, beGrp',myRank,beGrp
!mp
        Call CWTime(TCpu,TWall)
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',                  &
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',                  &
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',                      &
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',                      &
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',                      &
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
!mp
!
!          Extract H1(be',I) <- T1(be,I)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
!
!         vanish M1(m,i,u),M2(m,be',u),M3(m,be',u)
          dim1=nc*no*no
          call mv0zero (dim1,dim1,wrk(PossM1))
          dim1=nc*no*dimbe
          call mv0zero (dim1,dim1,wrk(PossM2))
          call mv0zero (dim1,dim1,wrk(PossM3))
!
!T1G          vanish H4(be',u)
          dim1=dimbe*no
          call mv0zero (dim1,dim1,wrk(PossH4))
!
!
!
           if (BetaID(myRank,beGrp).eq.1) then
!##       contributions Goo2,3  must be taken only once per Beta'
!
!Goo2.1            Extract H2(be',i) <- Fvo(be,i) (vlastnu rutinu netreba)
            call ExtT1 (wrk(PossH2),wrk(PossFvo),dimbe,addbe)
!
!Goo2.2f    Goo(i,u) <<- H2(T)(be',i) . H1(be',u)
            call mc0c1at3b (dimbe,no,dimbe,no,no,no,                    &
     &                    no,dimbe,no,                                  &
     &                    wrk(PossH2),wrk(PossH1),wrk(PossGoo))
!
!
!Goo3.1            read V1(be',I,JK) <- I1(be',I,JK)
            dim1=dimbe*no*no*(no+1)/2
            LunName=I1Name(beGrp)
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!Goo3.2            Make V2(be',j,i,u) <- 2 V1(be',j|iu) - V1(be',i|ju)
            call MkV_Goo3 (wrk(PossV1),wrk(PossV2),dimbe,no)
!
!Goo3.3f    Goo(i,u) <<- V2(T)(be',j,i,u) . H1(be,j)
            dim1=dimbe*no
!@x            dorobit rutinu  mv0v1at3u
!           call mv0v1at3u (dim1,no*no,dim1,no*no,
!    c                      no*no,dim1,1,1,
!    c                      wrk(PossV2),wrk(PossH1),wrk(PossGoo))
!
!            zatial odflaknute takto
             call mc0c1at3b (dim1,no*no,dim1,1,no*no,1,                 &
     &                     no*no,dim1,1,                                &
     &                     wrk(PossV2),wrk(PossH1),wrk(PossGoo))
!    c
!@x
           end if
!
!
          addb=0
!
          do bGrp=1,NvGrp
          dimb=DimGrpv(bGrp)
!
!G            Extract H2(b',I) <- T1(b,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dimb,addb)
!
!QK3.1      read V3(m,i,b') <- L1(m,i,b')
            LunName=L1Name(bGrp)
            dim1=nc*no*dimb
            call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
!
!QK3.3      M1(m,i,u)   <<-  V3(m,i  ,b') . H2(b',u)
            call mc0c1a3b(nc*no,dimb,dimb,no,nc*no,no,                  &
     &                    nc*no,dimb,no,                                &
     &                    wrk(PossV3),wrk(PossH2),wrk(PossM1))
!
!Q3.4       read V1(m,be',b')<- [V2(m,beb') or V2(m,b'be)]<- L2(m,be',b')
!            in cases be.le.b through interstep V2
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
!
!Q3.5       M2(m,be',u) <<-  V1(m,be',b') . H2(b',u)
            call mc0c1a3b(nc*dimbe,dimb,dimb,no,nc*dimbe,no,            &
     &                    nc*dimbe,dimb,no,                             &
     &                    wrk(PossV1),wrk(PossH2),wrk(PossM2))
!
!
!           read V2(be'b',I,J) <- T2(be',b,I,J)
            LunName=T2Name(beGrp,bGrp)
            if (bGrp.eq.beGrp) then
              dim1=dimbe*(dimbe+1)*no*no/2
              call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!             Expand V2(be',b',I,J) <- V1(be'b',I,J)
              dim1=dimb*(dimb+1)/2
              call ExpT2 (wrk(possV1),wrk(PossV2),dimbe,dim1,no)
            else
              dim1=dimbe*dimb*no*no
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
            end if
!
!           Make Tau(be',b',I,J) in V2(be',b',I,J)
               call MkTau_chcc (wrk(PossV2),wrk(PossH1),wrk(PossH2),    &
     &                       dimbe,dimb,no,1.0d0,1.0d0)

!           Make V1(i,b',be',u) <- 2 V2(be',b',u,i) - V2(be',b',i,u)
            call MkT_T17 (wrk(PossV1),wrk(PossV2),dimb,dimbe,no)
!
!           M3(m,be',u) <<- V3(m,i,b') . V1(i,b',be',u)
            dim1=dimb*no
            dim2=dimbe*no
            call mc0c1a3b (nc,dim1,dim1,dim2,nc,dim2,                   &
     &                     nc,dim1,dim2,                                &
     &                     wrk(PossV3),wrk(PossV1),wrk(PossM3))
!
!
          addb=addb+dimb
          end do
!
!
          if (BetaID(myRank,beGrp).eq.1) then
!##       contribution T16  must be taken only once per Beta'
!
!           Map V1(be',m,i) <- M2(m,be',i)
            call Map3_213 (wrk(PossM2),wrk(PossV1),nc,dimbe,no)
!
!T162.f     H4(be',u) <<- - V1(be',m,i) . M4(m,i,u)
            call mc0c2a3b (dimbe,nc*no,nc*no,no,dimbe,no,               &
     &                   dimbe,nc*no,no,                                &
     &                   wrk(PossV1),wrk(PossM4),wrk(PossH4))
!
          end if
!
!
!T27.1    V4(be',v,u,j ) <- M2(T)(m,be',v) . M4(m,u,j)
          dim1=dimbe*no*no*no
          call mv0zero (dim1,dim1,wrk(PossV4))
          call mc0c1at3b (nc,dimbe*no,nc,no*no,dimbe*no,no*no,          &
     &                    dimbe*no,nc,no*no,                            &
     &                    wrk(PossM2),wrk(PossM4),wrk(PossV4))
!
!
          adda=0
!
          do aGrp=1,NvGrp
          dima=DimGrpv(aGrp)
!
!         Test, if this Be' A' combination is to be run on this node
          if (BeAID(myRank,beGrp,aGrp).eq.0) goto 12
          if (printkey.ge.10) then
          write (6,*) ' Chol PID, be,a',myRank,beGrp,aGrp
          end if
!
!            Extract H2(a',I) <- T1(a,I)
            call ExtT1 (wrk(PossH2),wrk(PossT1o),dima,adda)
!
!           read Q(be',u,i,a'),K(be',u,i,a')
            LunName=Tmp1Name(beGrp,aGrp)
            dim1=dimbe*dima*no*no
            call GetX (wrk(PossQ),dim1,LunAux,LunName,1,0)
            call GetX (wrk(PossK),dim1,LunAux,LunName,0,1)
!
!           read V1(m,be',a')<- [V2(m,bea') or V2(m,a',be')] <- L2(m,be',a')
!            in cases be.le.b through interstep V2
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
!
!           V3(a',u) <- V1(T)(m,be',a')(L2) . M3(m,be',u)
            call mv0zero (dima*no,dima*no,wrk(PossV3))
            call mc0c1at3b (nc*dimbe,dima,nc*dimbe,no,dima,no,          &
     &                      dima,nc*dimbe,no,                           &
     &                      wrk(PossV1),wrk(PossM3),wrk(PossV3))
!
!T17.f      Add T1n(a',u) <<- V3(a',u)
                 call AdT_T17 (wrk(PossT1n),wrk(PossV3),                &
     &                    nv,dima,no,adda,1.0d0)
!
!           V3(be',a',i,u) <-  V1(m,be',a')(L2) . M1(m,i,u)
            dim1=dima*dimbe*no*no
            call mv0zero (dim1,dim1,wrk(PossV3))
            call mc0c1at3b (nc,dimbe*dima,nc,no*no,dimbe*dima,no*no,    &
     &                      dima*dimbe,nc,no*no,                        &
     &                      wrk(PossV1),wrk(PossM1),wrk(PossV3))
!
!           map V2(be',u,i,a') <-  V3(be',a',i,u)
            call Map4_1432 (wrk(PossV3),wrk(PossV2),dimbe,dima,no,no)
!
!           Q(be',u,i,a')  <<- -V2(be',u,i,a')
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossQ),1,-1.0d0)
!
!           Gvv(be',a') <<- 2 sum(i) [V2(be',i,i,a')
            call AdV_G2 (wrk(PossGvv),wrk(PossV2),                      &
     &                   nv,dimbe,dima,no,addbe,adda,2.0d0)
!
!           K(be',u,i,a')  <<-  V2(be',u,i,a')
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossK),1,1.0d0)
!
!           read M5(m,i,a') <- L1(m,i,a')
            LunName=L1Name(aGrp)
            dim1=dima*nc*no
            call GetX (wrk(PossM5),dim1,LunAux,LunName,1,1)
!
!           V2(be',u,i,a') <-  M2(m,be',u ) . M5(m,i,a')(L1)
            dim1=dima*dimbe*no*no
            call mv0zero (dim1,dim1,wrk(PossV2))
            call mc0c1at3b (nc,dimbe*no,nc,no*dima,dimbe*no,dima*no,    &
     &                      dimbe*no,nc,dima*no,                        &
     &                      wrk(PossM2),wrk(PossM5),wrk(PossV2))
!
!           Q(be',u,i,a') <<- 2 V2(be',u,i,a')
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossQ),1,2.0d0)
!
!x          write Q and K submatrix to corresponding files
            LunName=Tmp1Name(beGrp,aGrp)
            dim1=dimbe*dima*no*no
            call SaveX (wrk(PossQ),dim1,LunAux,LunName,1,0)
            call SaveX (wrk(PossK),dim1,LunAux,LunName,0,1)
!@@
!        call Chck_K (wrk(PossK),dimbe,addbe,dima,adda)
!        call Chck_Q (wrk(PossQ),dimbe,addbe,dima,adda)
!@@
!
!           Gvv(be',a') <<- - sum(i) [V2(be,i,i,a')
            call AdV_G2 (wrk(PossGvv),wrk(PossV2),                      &
     &                   nv,dimbe,dima,no,addbe,adda,-1.0d0)
!
!
!            Q array (alredy relesed) will be served for X
!
!T26.1f     Map Q(a',U,be',V) <-  2 . V2(be',V,U,a')
!           ie. here a' is a firs index, that means X (Q) is saved
!           in different order as K,Q
            call Map4_3421 (wrk(PossV2),wrk(PossQ),dimbe,no,no,dima)
            dim1=dima*no*no*dimbe
            call mv0sv (dim1,dim1,wrk(PossQ),2.0d0)
!
!T27.2      Map H3(j,a') <- H2(a',j)(T1)
            call Map2_21 (wrk(PossH2),wrk(PossH3),dima,no)
!
!T27.3f            V2(be',v,u,a') <- - V4(be',v,u,j) . H3(j,a')
            dim1=dimbe*no*no*dima
            call mv0zero (dim1,dim1,wrk(PossV2))
            dim1=dimbe*no*no
                 call mc0c2a3b (dim1,no,no,dima,dim1,dima,              &
     &                     dim1,no,dima,                                &
     &                     wrk(PossV4),wrk(PossH3),wrk(PossV2))
!
!T289.1            V1(m,j,v) <- M4(m,j,v) (L1)
            dim1=nc*no*no
            call mv0u (dim1,wrk(PossM4),1,wrk(PossV1),1)
!
!T289.2     V1(m,j,v) <<- M1(m,j,v)
            dim1=nc*no*no
            call mv0v1u (dim1,wrk(PossM1),1,wrk(PossV1),1,1.0d0)
!
!T289.3     V3(j,v,u,a') <- V1(T) (m,j,v) . M5(m,u,a')(L1)
            dim1=no*no*no*dima
            call mv0zero (dim1,dim1,wrk(PossV3))
               call mc0c1at3b (nc,no*no,nc,no*dima,no*no,no*dima,       &
     &                      no*no,nc,no*dima,                           &
     &                      wrk(PossV1),wrk(PossM5),wrk(PossV3))
!
!T289.4f    V2(be',v,u,a') <<- - H1(be',j)(T1) . V3(j,v,u,a')
            dim1=no*no*dima
                 call mc0c2a3b (dimbe,no,no,dim1,dimbe,dim1,            &
     &                     dimbe,no,dim1,                               &
     &                     wrk(PossH1),wrk(PossV3),wrk(PossV2))
!
!T2G        Map V1(a',U,be',V) <- V2(be',V,U,a')
            call Map4_3421 (wrk(PossV2),wrk(PossV1),dimbe,no,no,dima)
!
!T2G        Add Q(a',u,be',v)(X) <<- V1(a',U,be',V) (Faktor 2 koli X)
            dim1=dimbe*dima*no*no
            call mv0v1u (dim1,wrk(PossV1),1,wrk(PossQ),1,2.0d0)
!
!T2G        write X(a',u,be',v) = Q(a',u,be',v)
            LunName=Tmp2Name(aGrp,beGrp)
            dim1=dimbe*dima*no*no
            Xyes(aGrp,beGrp)=1
            call SaveX (wrk(PossQ),dim1,LunAux,LunName,1,1)
!
12        adda=adda+dima
          end do
!
!T1G        Add T1n(be,u) <<- H4(be',u)
        call AdT_T17 (wrk(PossT1n),wrk(PossH4),nv,dimbe,no,addbe,1.0d0)
!
11      addbe=addbe+dimbe
        end do
!
!
#ifdef _MOLCAS_MPP_
!##        Synchronizacny bod:
!        Allreduce Goo,Gvv
        dim1=no*no
        call gadgop (wrk(PossGoo),dim1,'+')
        dim1=nv*nv
        call gadgop (wrk(PossGvv),dim1,'+')
#endif
!
!
!Goo1.1f Add Goo(i,u) <<- Hoo(i,u)
        dim1=no*no
        call mv0v1u (dim1,wrk(PossHoo),1,wrk(PossGoo),1,1.0d0)
!
!
!Gvv1.1 Map V1(be,a) <- Hvv(a,be)
        call Map2_21 (wrk(PossHvv),wrk(PossV1),nv,nv)
!
!Gvv1.2f Add Gvv(be,a) <<- V1(be,a)
        dim1=nv*nv
        call mv0v1u (dim1,wrk(PossV1),1,wrk(PossGvv),1,1.0d0)
!
!
!Gvv2.1 Map V2(i,a) <- Fvo(a,i)
        call Map2_21 (wrk(PossFvo),wrk(PossV2),nv,no)
!
!Gvv2.2 Calc V1(be,a) <- t1(be,i) . V2(i,a)
        dim1=nv*nv
        call mv0zero (dim1,dim1,wrk(PossV1))
        call mc0c1a3b (nv,no,no,nv,nv,nv,                               &
     &                 nv,no,nv,                                        &
     &                 wrk(PossT1o),wrk(PossV2),wrk(PossV1))
!
!Gvv2.3f Add Gvv(be,a) <<- - V1(be,a)
        dim1=nv*nv
        call mv0v1u (dim1,wrk(PossV1),1,wrk(PossGvv),1,-1.0d0)
!
!@@
!        call Chck_Goo (wrk(PossGoo))
!        call Chck_Gvv (wrk(PossGvv))
!@@
!
        return
        end
