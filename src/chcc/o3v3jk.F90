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
        subroutine o3v3jk (wrk,wrksize,NvGrp,maxdim,LunAux)
!
!
!       This routine do:
!
!       --- Part I - generation of Q and K intermediates
!
!       Q(be,u,i,a) = 2 J(be,u,i,a)   - K(be,u,i,a)
!       Q(be,u,i,a) <-
!              Q1     + [ 2(u,be|a,i) - (u,i|a,be) ]
!              Q2     - [ 2(u,j |a,i) - (u,i|a,j ) ] . t(be,j)
!              Q4     + [ 2(a,i |b,j) - (a,j|b,i ) ] . D(be,b,j,u)
!
!       K(be,u,i,a) <-
!              K1     + [  (u,i |a,be)             ]
!              K2     - [  (u,i |a,j )             ] . t(be,j)
!              K4     - [  (b,i |a,j )             ] . T(be,b,j,u)
!
!       Hoo(i,u)  <<-
!       Hoo1        +  Foo(i,u)
!       Hoo2        +  sum(a,b,j) [2 (ai|bj) - (aj|bi)] Ta(a,b,u,j)
!        implemented as:
!       Hoo2        +  sum(a,b,j) [2 (bi|aj) - (bj|ai)] Ta(a,b,j,u)
!
!       Hvv(a,be) <<-
!       Hvv1        +  Fvv(a,be)
!       Hvv2        +  sum(b,i,j) [2 (ai|bj) - (aj|bi)] Ta(be,b,i,j)
!
!       Hvo(a,i)  <<-
!       Hvo1        +  Fvo(a,i)
!       Hvo2        +  sum(be,j)   [2 (be,j|ai) - (be,i|aj)] T1(be,j)
!
!        Aex(ij,u,v) <-
!       Aex1  + sum(a,b)  [  (ai|bj) . Ta(a,b,u,v)
!
!       A(ij,u,v) <-
!       A1    +   <-      [  (iu,jv)           ]
!       A2    + sum(a)    [  (aj|iu) . t1(a,v) ]
!       A3    + sum(a)    [  (ai|jv) . t1(a,u) ]
!       A4    +           [   Aex(ij,u,v)      ]
!
!        T1(be,u) <<-
!        T161 + sum(a,i)   [ (2(ai|be,u) . t(a,i) ]
!        T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(a,be,i,j) ]
!
!        T18 is implemented in the form:
!        T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(be,a,j,i) ]
!
!        where:
!       D(a,b,i,j) =  t(a,b,j,i) - T(a,b,i,j)
!       T(a,b,i,j) =  t(a,b,i,j)/2 + t(a,i) . t(b,j)
!        Ta(a,b,i,j)=  t(a,b,i,j)   + t(a,i) . t(b,j)
!
!
!1      intermediates Q and K will be temporarry stored in files
!       QFil, KFil as follows
!       do aGrp
!       do beGrp
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
        integer NvGrp,maxdim,LunAux

!       help variables
        integer dim1,dim2,dima,dimb,dimbe
        integer aGrp,bGrp,beGrp,adda,addb,addbe
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4,PossH5
        integer PossK,PossQ
!
        integer PossT

!
        character*6 LunName
!
!       --- introduction part ---
!
!x      Distribute memory
!        @@ yet distribution is the same for all three o3v3 drivers,
!        v buducnosti nosnost usposobit vsetky 3 vlastne
!
!
!
        PossT=PossFree
        call DistMemo3v3jk (NvGrp,maxdim,                               &
     &       PossV1,PossV2,PossV3,PossV4,                               &
     &       PossH1,PossH2,PossH3,PossH4,PossH5,                        &
     &       PossK,PossQ,                                               &
     &       PossT)
!
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
!mp!           write (6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
       write (6,*) ' Not Enough memory in o3v3jk step!',                &
     & 'Increase large and/or small segmentation ',                     &
     &                  (1.0d0*PossT)/(1.0d0*wrksize)
           call abend()
        end if
!
!       if (generkey.eq.1) then
!x      vytvorenie suboru T2       vektorov (docasne)
!       call UrobT2 (wrk(PossV1),NbeGrp2,NbGrp2,LunAux)
!       write (6,*) 'T2Vc    done'
!       end if
!
!
!
!##        Vanish Hoo,Hvo,Hvv, T1n, A, Aex (because of paralelization)
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
!
!
!*      cycle over be'
!
        addbe=0
!
        do beGrp=1,NvGrp
        dimbe=DimGrpv(beGrp)
!
!##        test, if something for this beGrp is planed to be run on this node
        dim1=0
        do dim2=1,NvGrp
          dim1=dim1+BeAID(myRank,beGrp,dim2)
        end do
!       dim1=dim1+BetaID(myRank,beGrp)
        if (dim1.eq.0) goto 11
        if (printkey.gt.1) then ! toto som si nie isty ...
        write (6,*) ' o3v3 JK - ID,beGrp',myRank,beGrp
        end if
!mp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
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
!
!T1G          vanish H4(be',u)
          dim1=dimbe*no
          call mv0zero (dim1,dim1,wrk(PossH4))
!
!G             Extract H1(be',J) <- T1(be,j)
          call ExtT1 (wrk(PossH1),wrk(PossT1o),dimbe,addbe)
!
!
!
!*        cycle over a'
!
          adda=0
!
          do aGrp=1,NvGrp
          dima=DimGrpv(aGrp)
!
!          Test, if this Be' A' combination is to be run on this node
          if (BeAID(myRank,beGrp,aGrp).eq.0) goto 12
          if (printkey.ge.10) then
          write (6,*) ' o3v3 JK - ID,be,a',myRank,beGrp,aGrp
          end if
!
!G             Extract H2(a',J) <- T1(a,j)
          call ExtT1 (wrk(PossH2),wrk(PossT1o),dima,adda)
!
!Gx         vanish Q(K)(be',u,i,a')

            dim1=no*no*dimbe*dima
            call mv0zero (dim1,dim1,wrk(PossQ))
            call mv0zero (dim1,dim1,wrk(PossK))
!
!QK2.1.1    read V2(a',o_a,JK) <- I2 (a',o_a|JK)
            LunName=I1name(aGrp)
            dim1=no*(no+1)*no*dima/2
            call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!QK2.1.2    Expand V1(a',o_a,J,K) <- V2(a',o_a,JK)
            dim1=no*(no+1)/2
                   call Exp2 (wrk(PossV2),wrk(PossV1),dima,no,dim1,no)
!
!QK2.2      map  V4(j,u,i,a') <- V1(a',j,u,i)
            call Map4_4123 (wrk(PossV1),wrk(PossV4),dima,no,no,no)
!
!
            if (aGrp.eq.beGrp) then
!            term A23 only for a'=be'
!
!A23.1        Calc V1(I,JK,L) <- V2(T)(a',I,JK) . H2(a',L)
              dim1=no*no*no*(no+1)/2
              call mv0zero (dim1,dim1,wrk(PossV1))
              dim1=no*no*(no+1)/2
              call mc0c1at3b (dima,dim1,dima,no,dim1,no,                &
     &                        dim1,dima,no,                             &
     &                        wrk(PossV2),wrk(PossH2),wrk(PossV1))
!
!A23.2f       Add A(ij,u,v) <<- V1(j,iu,v) + V1(i,jv,u)
              dim1=no*(no+1)/2
              call AdV_A23 (wrk(PossV1),wrk(PossA),dim1,no)
!
            end if
!
!
!QK1.1.1    read V2(bea',ui) = (be',a'|IJ)
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
!
!QK1.1.2    Expand V1(be',a',u,i) <- V2(bea',ui)
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
!
!QK1.2      map V2(be'u,i,a') <- V1(be',a',u,i)
            call Map4_1423 (wrk(PossV1),wrk(PossV2),dimbe,dima,no,no)
!
!Q1.3       Q(be',u,i,a) <-  - V2(be',u,i,a')
            dim1=no*no*dimbe*dima
            call mv0v1u (dim1,wrk(PossV2),1,wrk(PossQ),1,-1.0d0)
!
!K1.3f      K(be',u,i,a) <-  V2(be',u,i,a')
            dim1=no*no*dimbe*dima
            call mv0v1u(dim1,wrk(PossV2),1,wrk(PossK),1,1.0d0)
!
!Q1.4       read V1(be',o_be,a',o_a) = (be'I|a'J)
            LunName=I2name(beGrp,aGrp)
            dim1=no*no*dima*dimbe
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!Q1.5       Map V2(be',u,i,a') <- V1(be',u,a',i)
            call Map4_1243 (wrk(PossV1),wrk(PossV2),dimbe,no,dima,no)
!
!Q1.6f      Q(be',u,i,a) <-  2 V2(be',u,i,a')
            call mv0v1u(no*no*dimbe*dima,                               &
     &                  wrk(PossV2),1,wrk(PossQ),1,2.0d0)
!
!T161.1     V2(be',u) <- V1(be',u,a',i) . H2(a',i)
            dim1=dimbe*no
            dim2=dima*no
            call mv0zero (dim1,dim1,wrk(PossV2))
               call mv0v1a3u (dim1,dim2,dim2,dim1,                      &
     &                     dim1,dim2,1,1,                               &
     &                     wrk(PossV1),wrk(PossH2),wrk(PossV2))
!
!T161.2f    Add H4(be',u) <<- 2 V2(be',u)
            dim1=dimbe*no
               call mv0v1u (dim1,wrk(PossV2),1,wrk(PossH4),1,2.0d0)
!
!Hvo2.1            Make V2(a',i,be',j) <- [2 V1(be',j|a'i) - V1(be',i|a'j)]
            call MkV_Hvo2 (wrk(PossV1),wrk(PossV2),dimbe,dima,no)
!
!Hvo2.2     V3(a',i) <- V2(a',i,be',j) . H1(be',j)
            call mv0zero (dima*no,dima*no,wrk(PossV3))
            call mv0v1a3u (dima*no,dimbe*no,dimbe*no,dima*no,           &
     &                     dima*no,dimbe*no,1,1,                        &
     &                     wrk(PossV2),wrk(PossH1),wrk(PossV3))
!
!Hvo2.3f    Add Hvo(a,i) <<- V3(a',i) (da sa spravit aj vlastna rutina)
            call AdT_T17 (wrk(PossHvo),wrk(PossV3),                     &
     &                    nv,dima,no,adda,1.0d0)
!
!
!K2.3       map  V2(j,u,i,a') = -V4(j,u,i,a')
            dim1=no*no*no*dima
            call MkV_K22 (wrk(PossV2),wrk(PossV4),dim1)
!
!Q2.3       make V1(j,u,i,a') = -2V4(i,u,j,a')+V4(j,u,i,a') uz prepermutovanuo
            call MkV_Q22 (wrk(PossV4),wrk(PossV1),dima)
!
!Q2.4f      Q(be',u,i,a') <-  H1(be',j).V1(j,u,i,a')
            call mc0c1a3b (dimbe,no,no,no*no*dima,dimbe,no*no*dima,     &
     &                     dimbe,no,no*no*dima,                         &
     &                     wrk(PossH1),wrk(PossV1),wrk(PossQ))
!
!K2.4f      K(be',u,i,a') <- H1(be',j).V2(j,u,i,a')
            call mc0c1a3b (dimbe,no,no,no*no*dima,dimbe,no*no*dima,     &
     &                     dimbe,no,no*no*dima,                         &
     &                     wrk(PossH1),wrk(PossV2),wrk(PossK))
!
!Hvv2.1            vanish H3(be',a')
            dim1=dimbe*dima
            call mv0zero(dim1,dim1,wrk(PossH3))
!
!*            cycle over b'
!
              addb=0
!
              do bGrp=1,NvGrp
              dimb=DimGrpv(bGrp)
!
!Gx           Extract H5(b',i) <- T1(b,i)
              call ExtT1 (wrk(PossH5),wrk(PossT1o),dimb,addb)
!
!
!QK4.1.12     read V1(be'b',o_be,o_b) <- t2(be',b',I,J)
              LunName=T2name(beGrp,bGrp)
              if (bGrp.eq.beGrp) then
                dim1=dimbe*(dimbe+1)*no*no/2
                call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!QK4.1.2        Expand V1(be',b',p,q) <- V2(be'b',p,q)
                dim1=dimb*(dimb+1)/2
                call ExpT2 (wrk(PossV2),wrk(PossV1),dimbe,dim1,no)
              else
                dim1=dimbe*dimb*no*no
                call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
              end if
!
!
!        T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(be,a,j,i) ]
              if (aGrp.eq.bGrp) then
!              term T18 only in the case a'=b'
!              preserve: V1(be',a',I,J) = T2(be',a',I,J)
!                        V4(J,K,L,a') = (a',J|K,L)
!              distroy : V2,V3
!
!T18.1                Make V2(a',j,i,u) <- - [2(ai|ju)-(aj|iu)] from V4
                call MkV_T18 (wrk(PossV2),wrk(PossV4),dima,no)
!
!T18.2                Set V3(be',a',j,i) <- V1(be',a',j,i)
                dim1=dimbe*dima*no*no
                call mv0u(dim1,wrk(PossV1),1,wrk(PossV3),1)
!
!T18.3                Make Tau in V3
                   call MkTau_chcc (wrk(PossV3),wrk(PossH1),wrk(PossH2),&
     &                           dimbe,dima,no,1.0d0,1.0d0)
!
!T18.4f                Calc H4(be',u) <<- V3(be',a',j,i) . V2(a',j,i,u)
                   dim1=no*no*dima
                     call mc0c1a3b (dimbe,dim1,dim1,no,dimbe,no,        &
     &                         dimbe,dim1,no,                           &
     &                         wrk(PossV3),wrk(PossV2),wrk(PossH4))
!
              end if
!
!
!QK4.2        read V2(b',o_b,a',o_a) = (b'I|a'J)
              LunName=I2Name(bGrp,aGrp)
              dim1=dimb*dima*no*no
              call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
!
              if (aGrp.eq.beGrp) then
!              term A4 only for a'=be' ,If Tau is needed instead
!             of T, A4 term can be joined with Hoo2 subpart
!
!Aex1.1         Extract V3(a',b',ij) <- (ai|bj) from V2(b',J,a',I)
                dim1=no*(no+1)/2
                call MkV_A4 (wrk(PossV3),wrk(PossV2),dimb,dima,no,dim1)
!
!Aex1.2         Make Tau in V1 (@@@ toto sa menilo oproti T1=0 - OK)
                call MkTau_chcc (wrk(PossV1),wrk(PossH2),wrk(PossH5),   &
     &                           dima,dimb,no,1.0d0,1.0d0)
!
                if (intkey.eq.0) then
!                cholesky generation of integrals
!Aex1.3f          Aex(ij,u,v) <<- V3(T)(a',b',ij) . V1(a',b',u,v)
                  dim1=no*(no+1)/2
                  call mc0c1at3b (dimb*dima,dim1,dimb*dima,no*no,       &
     &                            dim1,no*no,                           &
     &                            dim1,dima*dimb,no*no,                 &
     &                            wrk(PossV3),wrk(PossV1),wrk(PossAex))
                else
!                W4 and W3 integrals from disc
!Aex(A4).3f       A(ij,u,v) <<- V3(T)(a',b',ij) . V1(a',b',u,v)
                  dim1=no*(no+1)/2
                  call mc0c1at3b (dimb*dima,dim1,dimb*dima,no*no,       &
     &                            dim1,no*no,                           &
     &                            dim1,dima*dimb,no*no,                 &
     &                            wrk(PossV3),wrk(PossV1),wrk(PossA))
                end if
!
!Aex1.4ff       nesuspended: V pripade Tau rekonstruovat naspat T2 vo V1
                call MkTau_chcc (wrk(PossV1),wrk(PossH2),wrk(PossH5),   &
     &                           dima,dimb,no,1.0d0,-1.0d0)
!
              end if
!
!Hvv2.2       make Tau(be',b',o_be,o_b) in V1
              call MkTau_chcc (wrk(PossV1),wrk(PossH1),wrk(PossH5),     &
     &                         dimbe,dimb,no,1.0d0,1.0d0)
!
!Hvv2.3       Make V3(b',i,j,a') <-2(a'i|b'j)-(a'j|b'i) from V2(b'J|a'I)
              call MkV_Hvv2 (wrk(PossV3),wrk(PossV2),dima,dimb,no)
!
!Hvv2.4       Calc H3(be',a') <<- V1(be',b',i,j) . V3(b',i,j,a')
              dim1=dimb*no*no
              call mc0c1a3b (dimbe,dim1,dim1,dima,dimbe,dima,           &
     &                       dimbe,dim1,dima,                           &
     &                       wrk(PossV1),wrk(PossV3),wrk(PossH3))
!
              if (aGrp.eq.beGrp) then
!              terms Hoo2 only for a'=be'
!              @ musim tu sachovat, lebo nemam dalsie V a tak
!                je tu jeden lacny perm navyse Hoo2.4 :-(, ak bude
!                vytvorene nove V5, tak  to treba prerobit @
!
!Hoo2.1                Map V3(a',b',j,u) <- V1(a',b',u,j)
                call Map3_132 (wrk(PossV1),wrk(PossV3),dima*dimb,no,no)
!
!Hoo2.2                Make V1(i,a',b',j) <- 2 V2(b'j|a'i) - V2(b'i|a'j)
                call MkV_Hoo2 (wrk(PossV1),wrk(PossV2),dima,dimb,no)
!
!Hoo2.3f        Hoo(i,u) <<- + V1(i,a',b',j) . V3(a',b',j,u)
                dim1=dima*dimb*no
                call mc0c1a3b (no,dim1,dim1,no,no,no,                   &
     &                         no,dim1,no,                              &
     &                         wrk(PossV1),wrk(PossV3),wrk(PossHoo))
!Hoo2.4post        Map V1(a',b',u,j) <- V3(a',b',j,u)
                call Map3_132 (wrk(PossV3),wrk(PossV1),dima*dimb,no,no)
!
              end if
!
!QK4.3        make T(be',b',o_be,o_b) in V1 (in V1 is Tau from Hvv2.2)
              call MkTau_chcc (wrk(PossV1),wrk(PossH1),wrk(PossH5),     &
     &                         dimbe,dimb,no,0.5d0,0.5d0)
!
!QK4.4        map V3(be',o_b,b',o_be) <-V1(be',b',o_be,o_b)(now T in V3)
              call Map4_1342 (wrk(PossV1),wrk(PossV3),dimbe,dimb,no,no)
!
!QK4.5        map V1(b',o_a,o_b,a') <- V2(b',o_b,a',o_a) (now I2 in V1)
              call Map4_1342 (wrk(PossV2),wrk(PossV1),dimb,no,dima,no)
!
!K4.6f        K(be',u,i,a') <<- - V3(be',u_b,b',j_be) . V1(b',j_a,i_b,a')
              call mc0c2a3b (dimbe*no,dimb*no,dimb*no,dima*no,          &
     &                       dimbe*no,dima*no,                          &
     &                       dimbe*no,dimb*no,dima*no,                  &
     &                       wrk(PossV3),wrk(PossV1),wrk(PossK))
!
!Q4.6         Make D: V2(be',u_b,b',j_be) (now D in V2)
!Q4.6c        from T - V3(be',u_b,b',j_be) and T1 - H1(be',j),H2(b',u)
!             velice specialna procedurka, ale da sa urobit fok :-)))
              call MkD_Q46                                              &
     &             (wrk(PossV2),wrk(PossV3),wrk(PossH1),wrk(PossH5),    &
     &              dimbe,dimb,no)
!
!Q4.7         Make V3(b',o_a,o_b,a') = 2(a',i |b',j) - (a',j|b',i )
!Q4.7c        from V1(b',o_a,o_b,a')
!             dalsia, o nieco menej korenista rutinka fok :-)))
              call MkI_Q47 (wrk(PossV3),wrk(PossV1),dimb,dima,no)
!
!Q4.8f        Q(be',u,i,a') <<- V2(be',u_b,b',j_be) . V3(b',j_a,i_b,a')
              call mc0c1a3b (dimbe*no,dimb*no,dimb*no,dima*no,          &
     &                       dimbe*no,dima*no,                          &
     &                       dimbe*no,dimb*no,dima*no,                  &
     &                       wrk(PossV2),wrk(PossV3),wrk(PossQ))
!
            addb=addb+dimb
            end do
!
!
!Hvv2.5f    Add Hvv(a,be) <<- - H3(be',a')
              call AdH_Hvv2 (wrk(PossH3),wrk(PossHvv),                  &
     &                     dima,dimbe,adda,addbe,nv)
!
!x          write Q and K submatrix to corresponding files
            LunName=Tmp1Name(beGrp,aGrp)
            dim1=dimbe*dima*no*no
            call SaveX (wrk(PossQ),dim1,LunAux,LunName,1,0)
            call SaveX (wrk(PossK),dim1,LunAux,LunName,0,1)
!
12        adda=adda+dima
          end do
!
!T1G        Add T1n(be,u) <<- H4(be',u)
        call AdT_T17(wrk(PossT1n),wrk(PossH4),nv,dimbe,no,addbe,1.0d0)
!
11      addbe=addbe+dimbe
        end do
!
!
#ifdef _MOLCAS_MPP_
!##        Synchronizacny bod:
!        Allreduce Hoo,Hvv,Hvo
!        Allreduce A,Aex
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
!
!Hoo1.1 Hoo(i,u) <<- Foo(i,u)
        dim1=no*nv
        call mv0v1u (dim1,wrk(PossFoo),1,wrk(PossHoo),1,1.0d0)
!
!Hvv1.1 Hvv(a,be) <<- Fvv(a,be)
        dim1=nv*nv
        call mv0v1u (dim1,wrk(PossFvv),1,wrk(PossHvv),1,1.0d0)
!
!Hvo1.1 Hvo(a,i) <<- Foo(a,i)
        dim1=no*nv
        call mv0v1u (dim1,wrk(PossFvo),1,wrk(PossHvo),1,1.0d0)
!
!A1.1   read V1(IJ,KL) <- I0(ij,kl)
        LunName=I0Name
        dim1=no*(no+1)*no*(no+1)/4
        call  GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
!
!A1.2f  Add A(ij,u,v) <<- (iu|jv) from V1(iu,jv)
        dim1=no*(no+1)/2
        call MkV_A1 (wrk(PossA),wrk(PossV1),dim1,no)
!
        if (intkey.eq.0) then
!A4.1f          Add A(ij,u,v) <<- Aex(ij,u,v)
          dim1=no*no*no*(no+1)/2
           call mv0v1u (dim1,wrk(PossAex),1,wrk(PossA),1,1.0d0)
        end if
!
!
!@@
!        call Chck_Hoo (wrk(PossHoo))
!        call Chck_Hvv (wrk(PossHvv))
!        call Chck_Hvo (wrk(PossHvo))
!        call Chck_A (wrk(PossA))
!@@
        return
        end
