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

subroutine o3v3jk(wrk,wrksize,NvGrp,maxdim,LunAux)
! This routine does:
!
! --- Part I - generation of Q and K intermediates
!
! Q(be,u,i,a) = 2 J(be,u,i,a)   - K(be,u,i,a)
! Q(be,u,i,a) <-
!        Q1     + [ 2(u,be|a,i) - (u,i|a,be) ]
!        Q2     - [ 2(u,j |a,i) - (u,i|a,j ) ] . t(be,j)
!        Q4     + [ 2(a,i |b,j) - (a,j|b,i ) ] . D(be,b,j,u)
!
! K(be,u,i,a) <-
!        K1     + [  (u,i |a,be)             ]
!        K2     - [  (u,i |a,j )             ] . t(be,j)
!        K4     - [  (b,i |a,j )             ] . T(be,b,j,u)
!
! Hoo(i,u)  <<-
! Hoo1        +  Foo(i,u)
! Hoo2        +  sum(a,b,j) [2 (ai|bj) - (aj|bi)] Ta(a,b,u,j)
!  implemented as:
! Hoo2        +  sum(a,b,j) [2 (bi|aj) - (bj|ai)] Ta(a,b,j,u)
!
! Hvv(a,be) <<-
! Hvv1        +  Fvv(a,be)
! Hvv2        +  sum(b,i,j) [2 (ai|bj) - (aj|bi)] Ta(be,b,i,j)
!
! Hvo(a,i)  <<-
! Hvo1        +  Fvo(a,i)
! Hvo2        +  sum(be,j)   [2 (be,j|ai) - (be,i|aj)] T1(be,j)
!
! Aex(ij,u,v) <-
! Aex1  + sum(a,b)  [  (ai|bj) . Ta(a,b,u,v)
!
! A(ij,u,v) <-
! A1    +   <-      [  (iu,jv)           ]
! A2    + sum(a)    [  (aj|iu) . t1(a,v) ]
! A3    + sum(a)    [  (ai|jv) . t1(a,u) ]
! A4    +           [   Aex(ij,u,v)      ]
!
!  T1(be,u) <<-
!  T161 + sum(a,i)   [ (2(ai|be,u) . t(a,i) ]
!  T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(a,be,i,j) ]
!
!  T18 is implemented in the form:
!  T18  - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(be,a,j,i) ]
!
! where:
! D(a,b,i,j) =  t(a,b,j,i) - T(a,b,i,j)
! T(a,b,i,j) =  t(a,b,i,j)/2 + t(a,i) . t(b,j)
! Ta(a,b,i,j)=  t(a,b,i,j)   + t(a,i) . t(b,j)
!
!1 intermediates Q and K will be temporarily stored in files
!  QFil, KFil as follows
!  do aGrp
!    do beGrp
!      Q(K) (be',u,i,a')
!    end do
!  end do
!
!2 Structure of files, where selected group of (pq|rs) are
!  stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
!
!  (A'I|JK)  I1inxx xx - Group of A'
!
!  (A'I|B'J) I2xxyy xx - Group of A'
!                   yy - Group of B'
!
!  (A'B'|IJ) I3xxyy xx - Group of A'
!                   yy - Group of B'
!
!3 Structure of Cholesky vector files
!
!  L1(m,I ,A')  L1vcxx xx - Group of A'
!
!  L2(m,A',B')  L2xxyy xx - Group of A'
!                      yy - Group of B'
!
!4 Structure of Amplitude file
!  t2(A',B',IJ)  T2xxyy xx - Group of A'
!                       yy - Group of B'

use Index_Functions, only: nTri_Elem
use chcc_global, only: BeAID, DimGrpv, I0Name, I1Name, I2Name, I3Name, intkey, no, nv, PosA, PosAex, PosFoo, PosFree, PosFvo, &
                       PosFvv, PosHoo, PosHvo, PosHvv, PosT1n, PosT1o, printkey, T2Name, TCpu, TCpu_l, Tmp1Name, TWall, TWall_l, &
                       TWall0
use Para_Info, only: MyRank
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NvGrp, maxdim, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: adda, addb, addbe, aGrp, beGrp, bGrp, dim_1, dim_2, dima, dimb, dimbe, PosH1, PosH2, PosH3, PosH4, PosH5, &
                     PosK, PosQ, PosT, PosV1, PosV2, PosV3, PosV4
character(len=6) :: LunName

! --- introduction part ---

!x Distribute memory
!  @@ yet distribution is the same for all three o3v3 drivers,
!  v buducnosti nosnost usposobit vsetky 3 vlastne

PosT = PosFree
call DistMemo3v3jk(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosH5,PosK,PosQ,PosT)
if (printkey >= 10) write(u6,*) ' Last Value :',PosT,wrksize
if (PosT > wrksize) then
  !mp write(u6,*) ' Nieje dobre, Dr. Ch.  Kokotopuloss',
  write(u6,*) ' Not Enough memory in o3v3jk step! Increase large and/or small segmentation ', &
              real(PosT,kind=wp)/real(wrksize,kind=wp)
  call abend()
end if

!if (generkey == 1) then
!  !vytvorenie suboru T2       vektorov (docasne)
!  call UrobT2(wrk(PosV1),NbeGrp2,NbGrp2,LunAux)
!  write(u6,*) 'T2Vc    done'
!end if

!## Vanish Hoo,Hvo,Hvv, T1n, A, Aex (because of paralelization)
dim_1 = no*no
wrk(PosHoo:PosHoo+dim_1-1) = Zero
dim_1 = no*nv
wrk(PosHvo:PosHvo+dim_1-1) = Zero
dim_1 = nv*nv
wrk(PosHvv:PosHvv+dim_1-1) = Zero
dim_1 = no*nv
wrk(PosT1n:PosT1n+dim_1-1) = Zero
dim_1 = no*no*nTri_Elem(no)
wrk(PosA:PosA+dim_1-1) = Zero
if (intkey == 0) wrk(PosAex:PosAex+dim_1-1) = Zero

! cycle over be'

addbe = 0

do beGrp=1,NvGrp
  dimbe = DimGrpv(beGrp)

  !## test, if something for this beGrp is planed to be run on this node
  dim_1 = sum(BeAID(myRank,beGrp,1:NvGrp))
  !dim_1 = dim_1+BetaID(myRank,beGrp)
  if (dim_1 /= 0) then
    if (printkey > 1) write(u6,*) ' o3v3 JK - ID,beGrp',myRank,beGrp ! toto som si nie isty ...
    !mp
    call CWTime(TCpu,TWall)
    if (printkey > 1) then
      write(u6,*)
      write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
      write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
      write(u6,*)
      write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
      write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
      write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
      write(u6,*)
    end if
    TCpu_l = TCpu
    TWall_l = TWall
    !mp

    !T1G vanish H4(be',u)
    dim_1 = dimbe*no
    wrk(PosH4:PosH4+dim_1-1) = Zero

    !G Extract H1(be',J) <- T1(be,j)
    call ExtT1(wrk(PosH1),wrk(PosT1o),dimbe,addbe)

    ! cycle over a'

    adda = 0

    do aGrp=1,NvGrp
      dima = DimGrpv(aGrp)

      ! Test, if this Be' A' combination is to be run on this node
      if (BeAID(myRank,beGrp,aGrp) /= 0) then
        if (printkey >= 10) write(u6,*) ' o3v3 JK - ID,be,a',myRank,beGrp,aGrp

        !G Extract H2(a',J) <- T1(a,j)
        call ExtT1(wrk(PosH2),wrk(PosT1o),dima,adda)

        !G vanish Q(K)(be',u,i,a')

        dim_1 = no*no*dimbe*dima
        wrk(PosQ:PosQ+dim_1-1) = Zero
        wrk(PosK:PosK+dim_1-1) = Zero

        !QK2.1.1 read V2(a',o_a,JK) <- I2 (a',o_a|JK)
        LunName = I1name(aGrp)
        dim_1 = nTri_Elem(no)*no*dima
        call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
        !QK2.1.2 Expand V1(a',o_a,J,K) <- V2(a',o_a,JK)
        dim_1 = nTri_Elem(no)
        call Exp2(wrk(PosV2),wrk(PosV1),dima,no,dim_1,no)

        !QK2.2 map  V4(j,u,i,a') <- V1(a',j,u,i)
        call Map4_4123(wrk(PosV1),wrk(PosV4),dima,no,no,no)

        if (aGrp == beGrp) then
          ! term A23 only for a'=be'

          !A23.1 Calc V1(I,JK,L) <- V2(T)(a',I,JK) . H2(a',L)
          dim_1 = no*no*nTri_Elem(no)
          wrk(PosV1:PosV1+dim_1-1) = Zero
          dim_1 = no*nTri_Elem(no)
          call mc0c1at3b(dima,dim_1,dima,no,dim_1,no,dim_1,dima,no,wrk(PosV2),wrk(PosH2),wrk(PosV1))

          !A23.2f Add A(ij,u,v) <<- V1(j,iu,v) + V1(i,jv,u)
          dim_1 = nTri_Elem(no)
          call AdV_A23(wrk(PosV1),wrk(PosA),dim_1,no)

        end if

        !QK1.1.1 read V2(bea',ui) = (be',a'|IJ)
        if (beGrp > aGrp) then
          LunName = I3name(beGrp,aGrp)
          dim_1 = nTri_Elem(no)*dima*dimbe
        else if (beGrp == aGrp) then
          LunName = I3name(beGrp,aGrp)
          dim_1 = nTri_Elem(no)*nTri_Elem(dima)
        else
          LunName = I3name(aGrp,beGrp)
          dim_1 = nTri_Elem(no)*dima*dimbe
        end if
        call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

        !QK1.1.2 Expand V1(be',a',u,i) <- V2(bea',ui)
        if (beGrp > aGrp) then
          dim_1 = nTri_Elem(no)
          call Exp2(wrk(PosV2),wrk(PosV1),dimbe,dima,dim_1,no)
        else if (beGrp == aGrp) then
          dim_1 = nTri_Elem(dima)
          dim_2 = nTri_Elem(no)
          call Exp4(wrk(PosV2),wrk(PosV1),dim_1,dima,dim_2,no)
        else
          dim_1 = nTri_Elem(no)
          call Exp2i(wrk(PosV2),wrk(PosV1),dima,dimbe,dim_1,no)
        end if

        !QK1.2 map V2(be'u,i,a') <- V1(be',a',u,i)
        call Map4_1423(wrk(PosV1),wrk(PosV2),dimbe,dima,no,no)

        !Q1.3 Q(be',u,i,a) <-  - V2(be',u,i,a')
        dim_1 = no*no*dimbe*dima
        wrk(PosQ:PosQ+dim_1-1) = wrk(PosQ:PosQ+dim_1-1)-wrk(PosV2:PosV2+dim_1-1)

        !K1.3f K(be',u,i,a) <-  V2(be',u,i,a')
        dim_1 = no*no*dimbe*dima
        wrk(PosK:PosK+dim_1-1) = wrk(PosK:PosK+dim_1-1)+wrk(PosV2:PosV2+dim_1-1)

        !Q1.4 read V1(be',o_be,a',o_a) = (be'I|a'J)
        LunName = I2name(beGrp,aGrp)
        dim_1 = no*no*dima*dimbe
        call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

        !Q1.5 Map V2(be',u,i,a') <- V1(be',u,a',i)
        call Map4_1243(wrk(PosV1),wrk(PosV2),dimbe,no,dima,no)

        !Q1.6f Q(be',u,i,a) <-  2 V2(be',u,i,a')
        dim_1 = no*no*dimbe*dima
        wrk(PosQ:PosQ+dim_1-1) = wrk(PosQ:PosQ+dim_1-1)+Two*wrk(PosV2:PosV2+dim_1-1)

        !T161.1 V2(be',u) <- V1(be',u,a',i) . H2(a',i)
        dim_1 = dimbe*no
        dim_2 = dima*no
        wrk(PosV2:PosV2+dim_1-1) = Zero
        call mv0v1a3u(dim_1,dim_2,dim_2,dim_1,dim_1,dim_2,1,1,wrk(PosV1),wrk(PosH2),wrk(PosV2))

        !T161.2f Add H4(be',u) <<- 2 V2(be',u)
        dim_1 = dimbe*no
        wrk(PosH4:PosH4+dim_1-1) = wrk(PosH4:PosH4+dim_1-1)+Two*wrk(PosV2:PosV2+dim_1-1)

        !Hvo2.1 Make V2(a',i,be',j) <- [2 V1(be',j|a'i) - V1(be',i|a'j)]
        call MkV_Hvo2(wrk(PosV1),wrk(PosV2),dimbe,dima,no)

        !Hvo2.2 V3(a',i) <- V2(a',i,be',j) . H1(be',j)
        wrk(PosV3:PosV3+dima*no-1) = Zero
        call mv0v1a3u(dima*no,dimbe*no,dimbe*no,dima*no,dima*no,dimbe*no,1,1,wrk(PosV2),wrk(PosH1),wrk(PosV3))

        !Hvo2.3f Add Hvo(a,i) <<- V3(a',i) (da sa spravit aj vlastna rutina)
        call AdT_T17(wrk(PosHvo),wrk(PosV3),nv,dima,no,adda,One)

        !K2.3 map  V2(j,u,i,a') = -V4(j,u,i,a')
        dim_1 = no*no*no*dima
        call MkV_K22(wrk(PosV2),wrk(PosV4),dim_1)

        !Q2.3 make V1(j,u,i,a') = -2V4(i,u,j,a')+V4(j,u,i,a') uz prepermutovanuo
        call MkV_Q22(wrk(PosV4),wrk(PosV1),dima)

        !Q2.4f Q(be',u,i,a') <-  H1(be',j).V1(j,u,i,a')
        call mc0c1a3b(dimbe,no,no,no*no*dima,dimbe,no*no*dima,dimbe,no,no*no*dima,wrk(PosH1),wrk(PosV1),wrk(PosQ))

        !K2.4f K(be',u,i,a') <- H1(be',j).V2(j,u,i,a')
        call mc0c1a3b(dimbe,no,no,no*no*dima,dimbe,no*no*dima,dimbe,no,no*no*dima,wrk(PosH1),wrk(PosV2),wrk(PosK))

        !Hvv2.1 vanish H3(be',a')
        dim_1 = dimbe*dima
        wrk(PosH3:PosH3+dim_1-1) = Zero

        ! cycle over b'

        addb = 0

        do bGrp=1,NvGrp
          dimb = DimGrpv(bGrp)

          !G Extract H5(b',i) <- T1(b,i)
          call ExtT1(wrk(PosH5),wrk(PosT1o),dimb,addb)

          !QK4.1.12 read V1(be'b',o_be,o_b) <- t2(be',b',I,J)
          LunName = T2name(beGrp,bGrp)
          if (bGrp == beGrp) then
            dim_1 = nTri_Elem(dimbe)*no*no
            call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)
            !QK4.1.2 Expand V1(be',b',p,q) <- V2(be'b',p,q)
            dim_1 = nTri_Elem(dimb)
            call ExpT2(wrk(PosV2),wrk(PosV1),dimbe,dim_1,no)
          else
            dim_1 = dimbe*dimb*no*no
            call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)
          end if

          !T18 - sum(a,i,j) [ (2(ai|ju)   - (aj|iu)  ) . Ta(be,a,j,i) ]
          if (aGrp == bGrp) then
            ! term T18 only in the case a'=b'
            ! preserve: V1(be',a',I,J) = T2(be',a',I,J)
            !           V4(J,K,L,a') = (a',J|K,L)
            ! destroy : V2,V3

            !T18.1 Make V2(a',j,i,u) <- - [2(ai|ju)-(aj|iu)] from V4
            call MkV_T18(wrk(PosV2),wrk(PosV4),dima,no)

            !T18.2 Set V3(be',a',j,i) <- V1(be',a',j,i)
            dim_1 = dimbe*dima*no*no
            wrk(PosV3:PosV3+dim_1-1) = wrk(PosV1:PosV1+dim_1-1)

            !T18.3 Make Tau in V3
            call MkTau_chcc(wrk(PosV3),wrk(PosH1),wrk(PosH2),dimbe,dima,no,One,One)

            !T18.4f Calc H4(be',u) <<- V3(be',a',j,i) . V2(a',j,i,u)
            dim_1 = no*no*dima
            call mc0c1a3b(dimbe,dim_1,dim_1,no,dimbe,no,dimbe,dim_1,no,wrk(PosV3),wrk(PosV2),wrk(PosH4))

          end if

          !QK4.2 read V2(b',o_b,a',o_a) = (b'I|a'J)
          LunName = I2Name(bGrp,aGrp)
          dim_1 = dimb*dima*no*no
          call GetX(wrk(PosV2),dim_1,LunAux,LunName,1,1)

          if (aGrp == beGrp) then
            ! term A4 only for a'=be' ,If Tau is needed instead
            ! of T, A4 term can be joined with Hoo2 subpart

            !Aex1.1 Extract V3(a',b',ij) <- (ai|bj) from V2(b',J,a',I)
            dim_1 = nTri_Elem(no)
            call MkV_A4(wrk(PosV3),wrk(PosV2),dimb,dima,no,dim_1)

            !Aex1.2 Make Tau in V1 (@@@ toto sa menilo oproti T1=0 - OK)
            call MkTau_chcc(wrk(PosV1),wrk(PosH2),wrk(PosH5),dima,dimb,no,One,One)

            if (intkey == 0) then
              ! cholesky generation of integrals
              !Aex1.3f Aex(ij,u,v) <<- V3(T)(a',b',ij) . V1(a',b',u,v)
              dim_1 = nTri_Elem(no)
              call mc0c1at3b(dimb*dima,dim_1,dimb*dima,no*no,dim_1,no*no,dim_1,dima*dimb,no*no,wrk(PosV3),wrk(PosV1),wrk(PosAex))
            else
              ! W4 and W3 integrals from disc
              !Aex(A4).3f A(ij,u,v) <<- V3(T)(a',b',ij) . V1(a',b',u,v)
              dim_1 = nTri_Elem(no)
              call mc0c1at3b(dimb*dima,dim_1,dimb*dima,no*no,dim_1,no*no,dim_1,dima*dimb,no*no,wrk(PosV3),wrk(PosV1),wrk(PosA))
            end if

            !Aex1.4ff nesuspended: V pripade Tau rekonstruovat naspat T2 vo V1
            call MkTau_chcc(wrk(PosV1),wrk(PosH2),wrk(PosH5),dima,dimb,no,One,-One)

          end if

          !Hvv2.2 make Tau(be',b',o_be,o_b) in V1
          call MkTau_chcc(wrk(PosV1),wrk(PosH1),wrk(PosH5),dimbe,dimb,no,One,One)

          !Hvv2.3 Make V3(b',i,j,a') <-2(a'i|b'j)-(a'j|b'i) from V2(b'J|a'I)
          call MkV_Hvv2(wrk(PosV3),wrk(PosV2),dima,dimb,no)

          !Hvv2.4 Calc H3(be',a') <<- V1(be',b',i,j) . V3(b',i,j,a')
          dim_1 = dimb*no*no
          call mc0c1a3b(dimbe,dim_1,dim_1,dima,dimbe,dima,dimbe,dim_1,dima,wrk(PosV1),wrk(PosV3),wrk(PosH3))

          if (aGrp == beGrp) then
            ! terms Hoo2 only for a'=be'
            ! @ musim tu sachovat, lebo nemam dalsie V a tak
            !   je tu jeden lacny perm navyse Hoo2.4 :-(, ak bude
            !   vytvorene nove V5, tak  to treba prerobit @

            !Hoo2.1 Map V3(a',b',j,u) <- V1(a',b',u,j)
            call Map3_132(wrk(PosV1),wrk(PosV3),dima*dimb,no,no)

            !Hoo2.2 Make V1(i,a',b',j) <- 2 V2(b'j|a'i) - V2(b'i|a'j)
            call MkV_Hoo2(wrk(PosV1),wrk(PosV2),dima,dimb,no)

            !Hoo2.3f Hoo(i,u) <<- + V1(i,a',b',j) . V3(a',b',j,u)
            dim_1 = dima*dimb*no
            call mc0c1a3b(no,dim_1,dim_1,no,no,no,no,dim_1,no,wrk(PosV1),wrk(PosV3),wrk(PosHoo))
            !Hoo2.4post Map V1(a',b',u,j) <- V3(a',b',j,u)
            call Map3_132(wrk(PosV3),wrk(PosV1),dima*dimb,no,no)

          end if

          !QK4.3 make T(be',b',o_be,o_b) in V1 (in V1 is Tau from Hvv2.2)
          call MkTau_chcc(wrk(PosV1),wrk(PosH1),wrk(PosH5),dimbe,dimb,no,Half,Half)

          !QK4.4 map V3(be',o_b,b',o_be) <-V1(be',b',o_be,o_b)(now T in V3)
          call Map4_1342(wrk(PosV1),wrk(PosV3),dimbe,dimb,no,no)

          !QK4.5 map V1(b',o_a,o_b,a') <- V2(b',o_b,a',o_a) (now I2 in V1)
          call Map4_1342(wrk(PosV2),wrk(PosV1),dimb,no,dima,no)

          !K4.6f K(be',u,i,a') <<- - V3(be',u_b,b',j_be) . V1(b',j_a,i_b,a')
          call mc0c2a3b(dimbe*no,dimb*no,dimb*no,dima*no,dimbe*no,dima*no,dimbe*no,dimb*no,dima*no,wrk(PosV3),wrk(PosV1),wrk(PosK))

          !Q4.6 Make D: V2(be',u_b,b',j_be) (now D in V2)
          !Q4.6c from T - V3(be',u_b,b',j_be) and T1 - H1(be',j),H2(b',u)
          ! velice specialna procedurka, ale da sa urobit fok :-)))
          call MkD_Q46(wrk(PosV2),wrk(PosV3),wrk(PosH1),wrk(PosH5),dimbe,dimb,no)

          !Q4.7 Make V3(b',o_a,o_b,a') = 2(a',i |b',j) - (a',j|b',i )
          !Q4.7c from V1(b',o_a,o_b,a')
          ! dalsia, o nieco menej korenista rutinka fok :-)))
          call MkI_Q47(wrk(PosV3),wrk(PosV1),dimb,dima,no)

          !Q4.8f Q(be',u,i,a') <<- V2(be',u_b,b',j_be) . V3(b',j_a,i_b,a')
          call mc0c1a3b(dimbe*no,dimb*no,dimb*no,dima*no,dimbe*no,dima*no,dimbe*no,dimb*no,dima*no,wrk(PosV2),wrk(PosV3),wrk(PosQ))

          addb = addb+dimb
        end do

        !Hvv2.5f Add Hvv(a,be) <<- - H3(be',a')
        call AdH_Hvv2(wrk(PosH3),wrk(PosHvv),dima,dimbe,adda,addbe,nv)

        ! write Q and K submatrix to corresponding files
        LunName = Tmp1Name(beGrp,aGrp)
        dim_1 = dimbe*dima*no*no
        call SaveX(wrk(PosQ),dim_1,LunAux,LunName,1,0)
        call SaveX(wrk(PosK),dim_1,LunAux,LunName,0,1)

      end if
      adda = adda+dima
    end do

    !T1G Add T1n(be,u) <<- H4(be',u)
    call AdT_T17(wrk(PosT1n),wrk(PosH4),nv,dimbe,no,addbe,One)

  end if
  addbe = addbe+dimbe
end do

#ifdef _MOLCAS_MPP_
!## Synchronizacny bod:
! Allreduce Hoo,Hvv,Hvo
! Allreduce A,Aex
dim_1 = no*no
call gadgop(wrk(PosHoo),dim_1,'+')
dim_1 = nv*nv
call gadgop(wrk(PosHvv),dim_1,'+')
dim_1 = no*nv
call gadgop(wrk(PosHvo),dim_1,'+')
dim_1 = no*no*nTri_Elem(no)
call gadgop(wrk(PosA),dim_1,'+')
if (intkey == 0) call gadgop(wrk(PosAex),dim_1,'+')
#endif

!Hoo1.1 Hoo(i,u) <<- Foo(i,u)
dim_1 = no*nv
wrk(PosHoo:PosHoo+dim_1-1) = wrk(PosHoo:PosHoo+dim_1-1)+wrk(PosFoo:PosFoo+dim_1-1)

!Hvv1.1 Hvv(a,be) <<- Fvv(a,be)
dim_1 = nv*nv
wrk(PosHvv:PosHvv+dim_1-1) = wrk(PosHvv:PosHvv+dim_1-1)+wrk(PosFvv:PosFvv+dim_1-1)

!Hvo1.1 Hvo(a,i) <<- Foo(a,i)
dim_1 = no*nv
wrk(PosHvo:PosHvo+dim_1-1) = wrk(PosHvo:PosHvo+dim_1-1)+wrk(PosFvo:PosFvo+dim_1-1)

!A1.1 read V1(IJ,KL) <- I0(ij,kl)
LunName = I0Name
dim_1 = nTri_Elem(no)**2
call GetX(wrk(PosV1),dim_1,LunAux,LunName,1,1)

!A1.2f Add A(ij,u,v) <<- (iu|jv) from V1(iu,jv)
dim_1 = nTri_Elem(no)
call MkV_A1(wrk(PosA),wrk(PosV1),dim_1,no)

if (intkey == 0) then
  !A4.1f Add A(ij,u,v) <<- Aex(ij,u,v)
  dim_1 = no*no*nTri_Elem(no)
  wrk(PosA:PosA+dim_1-1) = wrk(PosA:PosA+dim_1-1)+wrk(PosAex:PosAex+dim_1-1)
end if

!@@
!call Chck_Hoo(wrk(PosHoo))
!call Chck_Hvv(wrk(PosHvv))
!call Chck_Hvo(wrk(PosHvo))
!call Chck_A(wrk(PosA))
!@@

return

end subroutine o3v3jk
