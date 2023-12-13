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

subroutine o2v4(wrk,wrksize,NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,LunAux)
! this routine does:
!
! t(u,v,bega) <- sum(a,b) [ b(a,b,be,ga) . Tau(u,v,a,b) ]
!
! b(a,b,be,ga) =             (a,be|b,ga)
!              - sum (i)   [ (a,be|b,i) . t(i,ga) ]
!              - sum (i)   [ (a,i|b,ga) . t(i,be) ]
!              + sum (i,j) [ (a,i|b,j) . t(i,be) . t(j,ga) ]
!
! Tau(u,v,ab) = t(u,v,ab) + t(i,a) . t(j,b) : stored for: i,j,a>=b
!
!
!A List of main arrays:
!
!* tau(i,j,(a,b)') - PosTau
!
!* t2n1(ij,(bega)")- PosT2n1
!  t2n2(ij,(bega)")- PosT2n2
!  t2w(ij,(ab)")   - PosT2w - work file for T(+-)
!
!* L21 (m,a',be')  - PosL21
!  L22 (m,a',ga')  - PosL22
!  L23 (m,b',be')  - PosL23
!  L24 (m,b',ga')  - PosL24
!  L2W (m,c',de')  - PosL2W - Used for Tmp file in GetCHV
!  some (or all) of L2s might be identical
!
!* M1(m,a",be")  - PosM1, used also for M(m,b'',be'')
!  M2(m,b",ga")  - PosM2, used also for M(m,a'',ga'')
!  Both Ms might be identical (mozno + M3 ako manipulacne)
!
!* W1(a",be",b",ga") - PosW1 - b(a",be"|b",ga")
!  W2(b",be",a",ga") - PosW2 - b(b",be"|a",ga")
!  Ww((ab)",(bega)") - PosWw - work file for b(+-)
!
!B list of additional arrays
!
!C list of used pre- and suffixes and assignements
!
!  _Grp         - Group
!  _SGrp        - SubGroup
!  Pos_         - Position
!  PsAc_        - Actual Position
!  (...)'       - Group, Block
!  (...)"       - SubGroup, SubBlock
!
!D list of most important variables

use Index_Functions, only: nTri_Elem
use chcc_global, only: ABID, DimGrpa, DimGrpbe, DimSGrpa, DimSGrpbe, GrpaLow, GrpbeLow, GrpaUp, GrpbeUp, intkey, L1Name, maxSGrp, &
                       nc, no, nv, PosFree, PosT1o, printkey, T2o2v4yes, Tmp3Name
use Para_Info, only: MyRank
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NaGrp, NbeGrp, NaSGrp, NbeSgrp, mdGrpa, mdGrpbe, mdSGrpa, mdSGrpbe, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: adda, addb, addbe, addbepp, addga, addgapp, aGrp, aSGrp, beGrp, beSGrp, bGrp, bSGrp, bSGrpUp, choleskykey, &
                     dim_1, dim_2, dim_3, dima, dimb, dimbe, dimga, FirstT2n(maxSGrp,maxSGrp), gaGrp, gaSGrp, gaSGrpUp, i, &
                     L2Status(4,3), lent2n1, lent2n2, NL2, pL21, pL22, pL23, pL24, PosH1, PosH2, PosL11, PosL12, PosL21, PosL22, &
                     PosL23, PosL24, PosL2W, PosM1, PosM2, PosT, PosT2n1, PosT2n2, PosT2w, PosTau, PosW1, PosW2, PosWw, PosWx, &
                     PsAcL21, PsAcL22, PsAcL23, PsAcL24
character(len=6) :: LunName

if (intkey == 1) then
  choleskykey = 0
else
  choleskykey = 1
end if

! distribute memory

PosT = PosFree
call DistMemo2v4(NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,PosTau,PosT2n1,PosT2n2,PosT2w,PosL11,PosL12,PosL21, &
                 PosL22,PosL23,PosL24,PosL2W,PosH1,PosH2,PosM1,PosM2,PosW1,PosW2,PosWw,PosWx,PosT,NL2)
if (printkey >= 10) write(u6,*) ' Last Value :',PosT,wrksize
if (PosT > wrksize) then
  !mp write(u6,*) ' Nieje dobre - o2v4, Dr. Ch. Kokotopuloss',
  write(u6,*) ' Not Enough memory in o2v4 step! Increase large and/or small segmentation ',real(PosT,kind=wp)/real(wrksize,kind=wp)
  call abend()
end if

!@@
! call Calc_Bc
!@@

! initialize L2Status

L2Status(1:NL2,:) = 0

if (NL2 == 1) then
  L2Status(1,3) = posL21
else if (NL2 == 2) then
  if (NbeGrp == 1) then
    L2Status(1,3) = posL21
    L2Status(2,3) = posL23
  else
    L2Status(1,3) = posL21
    L2Status(2,3) = posL22
  end if
else
  L2Status(1,3) = posL21
  L2Status(2,3) = posL22
  L2Status(3,3) = posL23
  L2Status(4,3) = posL24
end if

! initialize FirstT2n and Generate Tmp3Name s
!## initialize T2o2v4yes

do beGrp=1,nbeGrp
  do gaGrp=1,beGrp

    do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
      if (beGrp == gaGrp) then
        gaSGrpUp = beSGrp
      else
        gaSGrpUp = GrpbeUp(gaGrp)
      end if
      do gaSGrp=GrpbeLow(gaGrp),gaSGrpUp

        FirstT2n(beSGrp,gaSGrp) = 1

        T2o2v4yes(beSGrp,gaSGrp) = 0

      end do
    end do

  end do
end do

! cycle over all groups of (a>=b)
adda = 0
do aGrp=1,NaGrp
  dima = DimGrpa(aGrp)

  !## test, if on this node at least one combination with this
  ! aGrp is scheduled. Skip if no
  i = sum(ABID(myRank,aGrp,1:NaGrp))

  if (i /= 0) then

    if (choleskykey == 1) then
      ! read L2W(m,i,a') <- L1(m,i,a')
      LunName = L1Name(aGrp)
      dim_1 = nc*dima*no
      call GetX(wrk(PosL2W),dim_1,LunAux,LunName,1,1)
      ! Map L11(m,a',i) <- L2W(m,i,a')
      call Map3_132(wrk(PosL2W),wrk(PosL11),nc,no,dima)
    end if

    addb = 0
    do bGrp=1,aGrp
      dimb = DimGrpa(bGrp)

      !## test, if this a'b' combination is planed to be run on this node
      if (ABID(myRank,aGrp,bGrp) /= 0) then
        if (printkey >= 10) write(u6,*) ' O2V4 MyRank, aGrp, bGrp',myRank,aGrp,bGrp

        if (choleskykey == 1) then
          ! read L12(m,b',i) <- L1(m,b',i)
          LunName = L1Name(bGrp)
          dim_1 = nc*dimb*no
          if (NaGrp > 1) then
            if (aGrp == bGrp) then
              wrk(PosL12:PosL12+dim_1-1) = wrk(PosL11:PosL11+dim_1-1)
            else
              call GetX(wrk(PosL2W),dim_1,LunAux,LunName,1,1)
              call Map3_132(wrk(PosL2W),wrk(PosL12),nc,no,dimb)
            end if
          end if
        end if

        ! read the block of  amplitudes T2((ab)',ij) for given aGrp,bGrp
        ! and make Tau from them
        call getTau(wrk(PosTau),wrk(PosT1o),aGrp,bGrp,dima,dimb,adda,addb,LunAux)

        ! cycle over all groups of (be>=ga)
        addbe = 0
        do beGrp=1,nbeGrp
          dimbe = DimGrpbe(beGrp)
          addga = 0
          do gaGrp=1,beGrp
            dimga = DimGrpbe(gaGrp)

            if (choleskykey == 1) then
              ! read Cholesky vectors
              !  L21 (m,a',be')
              !  L22 (m,a',ga')
              !  L23 (m,b',be')
              !  L24 (m,b',ga')
              call getChV(wrk,wrksize,aGrp,bGrp,beGrp,gaGrp,NL2,L2Status,pL21,pL22,pL23,pL24,PosL2W,PosL11,PosL12,LunAux)
              PsAcL21 = L2Status(pL21,3)
              PsAcL22 = L2Status(pL22,3)
              PsAcL23 = L2Status(pL23,3)
              PsAcL24 = L2Status(pL24,3)
            else
              PsAcL21 = 0
              PsAcL22 = 0
              PsAcL23 = 0
              PsAcL24 = 0
            end if

            ! cycle over all subgroups of (be>=ga)'
            addbepp = addbe
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)

              ! get H1(i,be") <- t1o(be,i)
              if (intkey == 1) call Mk_T1t(wrk(PosT1o),wrk(PosH1),DimSGrpbe(beSGrp),no,nv,addbepp)

              if (beGrp == gaGrp) then
                gaSGrpUp = beSGrp
              else
                gaSGrpUp = GrpbeUp(gaGrp)
              end if

              addgapp = addga
              do gaSGrp=GrpbeLow(gaGrp),gaSGrpUp

                ! get H2(i,ga") <- t1o(ga,i)
                if (intkey == 1) call Mk_T1t(wrk(PosT1o),wrk(PosH2),DimSGrpbe(gaSGrp),no,nv,addgapp)

                ! vanish arrays for new (final) amplitudes, if this is a
                !   first use, or read from Tmp3Name(be",ga") actual stage
                ! T2n1((bega)",ij)
                ! T2n2((bega)",ij)
                if (FirstT2n(beSGrp,gaSgrp) == 1) then
                  call VanishT2n(wrk(PosT2n1),wrk(PosT2n2),beSGrp,gaSGrp)
                  FirstT2n(beSGrp,gaSgrp) = 0
                else
                  call GetT2n(wrk(PosT2n1),wrk(PosT2n2),beSGrp,gaSGrp,LunAux)
                end if

                ! cycle over all subgroups of (a>=b)'
                do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
                  if (aGrp == bGrp) then
                    bSGrpUp = aSGrp
                  else
                    bSGrpUp = GrpaUp(bGrp)
                  end if
                  do bSGrp=GrpaLow(bGrp),bSGrpUp
                    if (printkey >= 10) write(u6,99) aGrp,bGrp,beGrp,gaGrp,aSGrp,bSGrp,beSGrp,gaSGrp

                    if (choleskykey == 1) then

                      ! cholesky generation of (VV|VV) integrals

                      ! Extract M1(m,a",be") from L21(m,a',be')
                      call ExtractM(wrk(PosM1),wrk(PsAcL21),aGrp,beGrp,aSGrp,beSGrp)
                      ! Extract M2(m,b",ga") from L24(m,b',ga')
                      if ((aSGrp == bSGrp) .and. (beSGrp == gaSGrp)) then
                        if (PosM1 /= PosM2) then
                          dim_1 = nc*DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)
                          wrk(PosM2:PosM2+dim_1-1) = wrk(PosM1:PosM1+dim_1-1)
                        end if
                      else
                        call ExtractM(wrk(PosM2),wrk(PsAcL24),bGrp,gaGrp,bSGrp,gaSGrp)
                      end if
                      ! Calc W1(a",be",b",ga")=M1(m,a",be")(T).M2(m,b",ga")
                      dim_1 = DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)
                      dim_2 = DimSGrpa(bSGrp)*DimSGrpbe(gaSGrp)
                      wrk(PosW1:PosW1+dim_1*dim_2-1) = Zero
                      call mc0c1at3b(nc,dim_1,nc,dim_2,dim_1,dim_2,dim_1,nc,dim_2,wrk(posM1),wrk(posM2),wrk(posW1))

                      if (aSGrp /= bSGrp) then
                        ! Extract M1(m,b",be") from L23(m,b',be')
                        call ExtractM(wrk(PosM1),wrk(PsAcL23),bGrp,beGrp,bSGrp,beSGrp)
                        ! Extract M2(m,a",ga") from L22(m,a',ga')
                        call ExtractM(wrk(PosM2),wrk(PsAcL22),aGrp,gaGrp,aSGrp,gaSGrp)
                        ! Calc W2(b",be",a",ga")=M1(m,b",be")(T).M2(m,a",ga")
                        dim_1 = DimSGrpa(bSGrp)*DimSGrpbe(beSGrp)
                        dim_2 = DimSGrpa(aSGrp)*DimSGrpbe(gaSGrp)
                        wrk(posW2:posW2+dim_1*dim_2-1) = Zero
                        call mc0c1at3b(nc,dim_1,nc,dim_2,dim_1,dim_2,dim_1,nc,dim_2,wrk(posM1),wrk(posM2),wrk(posW2))
                      end if

                    else

                      ! clasical reading or generating of (VV|VV) integrals

                      ! term W1(a",be",b",ga") <- -(b"ga"|a"i) . T1(i,be")
                      ! Def:W1 destroy:Ww
                      ! Get Ww(b",ga",a",i) (Wx used as Aux)
                      call ReaW3(wrk(PosWw),wrk(PosWx),bSGrp,gaSGrp,aSGrp,LunAux)
                      ! Map V1(a",ga",b",i) <- Ww(b",ga",a",i)
                      call Map4_3214(wrk(PosWw),wrk(PosW1),DimSGrpa(bSGrp),DimSGrpbe(gaSGrp),DimSGrpa(aSGrp),no)
                      ! Calc Ww(a",ga",b",be") <- - W1(a",ga",b",i) . H1(i,be")
                      dim_1 = DimSGrpa(bSGrp)*DimSGrpa(aSGrp)*DimSGrpbe(gaSGrp)
                      dim_2 = DimSGrpbe(beSGrp)
                      wrk(PosWw:PosWw+dim_1*dim_2-1) = Zero
                      call mc0c2a3b(dim_1,no,no,dim_2,dim_1,dim_2,dim_1,no,dim_2,wrk(PosW1),wrk(PosH1),wrk(PosWw))
                      ! Map W1(a",be",b",ga") <<- Ww(a",ga",b",be")
                      call Map4_1432(wrk(PosWw),wrk(PosW1),DimSGrpa(aSGrp),DimSGrpbe(gaSGrp),DimSGrpa(bSGrp),DimSGrpbe(beSGrp))

                      ! term W1(a",be",b",ga") <<- -(a"be"|b"i) . T1(i,ga")
                      ! Upgrade: W1, destroy: Ww
                      ! Get Ww(a",be",b",i) (Wx used as Aux)
                      call ReaW3(wrk(PosWw),wrk(PosWx),aSGrp,beSGrp,bSGrp,LunAux)
                      ! Add W1(a",be",b",ga") <<- - Ww(a",be",b",i) . H2(i,ga")
                      dim_1 = DimSGrpa(aSGrp)*DimSGrpa(bSGrp)*DimSGrpbe(beSGrp)
                      dim_2 = DimSGrpbe(gaSGrp)
                      call mc0c2a3b(dim_1,no,no,dim_2,dim_1,dim_2,dim_1,no,dim_2,wrk(PosWw),wrk(PosH2),wrk(PosW1))

                      ! Upgrade W1(a",be",b",ga") <<- (a",be"|b",ga")
                      ! Upgrade:W1 destroy:Ww
                      call ReaW4(wrk(PosW1),wrk(PosWw),aSGrp,beSGrp,bSGrp,gaSGrp,LunAux)

                      if (aSGrp /= bSGrp) then

                        ! term W2(b",be",a",ga") <- -(a"ga"|b"i) . T1(i,ga")
                        ! Def:W2 destroy: Ww
                        ! Get Ww(a",ga",b",i) (Wx used as Aux)
                        call ReaW3(wrk(PosWw),wrk(PosWx),aSGrp,gaSGrp,bSGrp,LunAux)
                        ! Map W2(b",ga",a",i) <- Ww(a",ga",b",i) @@Nepreskusany
                        call Map4_3214(wrk(PosWw),wrk(PosW2),DimSGrpa(aSGrp),DimSGrpbe(gaSGrp),DimSGrpa(bSGrp),no)
                        ! Calc Ww(b",ga",a",be") <- -W2(b",ga",a",i) . H1(i,be")
                        dim_1 = DimSGrpa(aSGrp)*DimSGrpbe(gaSGrp)*DimSGrpa(bSGrp)
                        dim_2 = DimSGrpbe(beSGrp)
                        wrk(PosWw:PosWw+dim_1*dim_2-1) = Zero
                        call mc0c2a3b(dim_1,no,no,dim_2,dim_1,dim_2,dim_1,no,dim_2,wrk(PosW2),wrk(PosH1),wrk(PosWw))
                        ! Map W2 (b",be",a",ga") <- Ww(b",ga",a",be")
                        call Map4_1432(wrk(PosWw),wrk(PosW2),DimSGrpa(bSGrp),DimSGrpbe(gaSGrp),DimSGrpa(aSGrp),DimSGrpbe(beSGrp))

                        ! term W2(b",be",a",ga") <<- -(b"be"|a"i) . T1(i,ga")
                        ! Upgrade:W2 destroy: 0
                        ! Get Ww(b",be",a",i) (Wx used as Aux)
                        call ReaW3(wrk(PosWw),wrk(PosWx),bSGrp,beSGrp,aSGrp,LunAux)
                        ! Add W2(b",be",a",ga") <<- - Ww(b",be",a",i) .H2(i,ga")
                        dim_1 = DimSGrpa(aSGrp)*DimSGrpa(bSGrp)*DimSGrpbe(beSGrp)
                        dim_2 = DimSGrpbe(gaSGrp)
                        call mc0c2a3b(dim_1,no,no,dim_2,dim_1,dim_2,dim_1,no,dim_2,wrk(PosWw),wrk(PosH2),wrk(PosW2))

                        ! Add W2(b",be",a",ga") <<- (b",be"Ia",ga")
                        ! Upgrade:W2, destroy:Ww
                        call ReaW4(wrk(PosW2),wrk(PosWw),bSGrp,beSGrp,aSGrp,gaSGrp,LunAux)
                      end if

                    end if

                    ! Make t2w(+)((ab)",ij) from Tau((ab)',ij)
                    call MakeT2p(wrk(PosT2w),wrk(PosTau),aGrp,bGrp,aSGrp,bSGrp,1)
                    ! Make Ww(+)((ab)",(bega)") from W1(a",be",b",ga")
                    !                            and W2(b",be",a",ga")
                    !      N.B. for a"=b" we have only W1 defined (=W2)
                    !           but never mind (solved within MakeWw)
                    call MakeWw(wrk(PosWw),wrk(PosW1),wrk(PosW2),aSGrp,bSGrp,beSGrp,gaSGrp,1)
                    ! Add T2n1((bega)",ij) <<-
                    ! Ww(+)(T)((ab)",(bega)").t2w(+)((ab)",ij)
                    dim_1 = nTri_Elem(no)
                    if (aSGrp == bSGrp) then
                      dim_2 = DimSGrpa(aSGrp)*(DimSGrpa(bSGrp)-1)/2
                    else
                      dim_2 = DimSGrpa(aSGrp)*DimSGrpa(bSGrp)
                    end if
                    if (beSGrp == gaSGrp) then
                      dim_3 = DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
                    else
                      dim_3 = DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                    end if
                    call mc0c1at3b(dim_2,dim_3,dim_2,dim_1,dim_3,dim_1,dim_3,dim_2,dim_1,wrk(posWw),wrk(posT2w),wrk(posT2n1))

                    if (aSGrp == bSGrp) then
                      ! Make t2w(+)((aa)",ij) from Tau((ab)',ij)
                      call MakeT2pd(wrk(PosT2w),wrk(PosTau),aGrp,aSGrp)
                      ! Make  Ww(+)((aa)",(bega)") from W1(a",be",a",ga")
                      call MakeWwd(wrk(PosWw),wrk(PosW1),aSGrp,beSGrp,gaSGrp)
                      ! Add T2n1((bega)",ij) <<-
                      ! Ww(+)(T)((aa)",(bega)").t2w(+)((aa)",ij)
                      dim_1 = nTri_Elem(no)
                      dim_2 = DimSGrpa(aSGrp)
                      if (beSGrp == gaSGrp) then
                        dim_3 = DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
                      else
                        dim_3 = DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                      end if
                      call mc0c1at3b(dim_2,dim_3,dim_2,dim_1,dim_3,dim_1,dim_3,dim_2,dim_1,wrk(posWw),wrk(posT2w),wrk(posT2n1))
                    end if

                    ! Make t2w(-)((ab)",ij) from Tau((ab)',ij)
                    call MakeT2m(wrk(PosT2w),wrk(PosTau),aGrp,bGrp,aSGrp,bSGrp,1)
                    ! Make Ww(-)((ab)",(bega)") from W1(a",be",b",ga")
                    !                            and W2(b",be",a",ga")
                    !           but never mind (solved within MakeWw)
                    call MakeWw(wrk(PosWw),wrk(PosW1),wrk(PosW2),aSGrp,bSGrp,beSGrp,gaSGrp,2)
                    !write(u6,*) 'V'
                    ! Add T2n2((bega)",ij) <<-
                    ! Ww(-)(T)((ab)",(bega)").t2w(-)((ab)",ij)
                    dim_1 = nTri_Elem(no-1)
                    if (aSGrp == bSGrp) then
                      dim_2 = DimSGrpa(aSGrp)*(DimSGrpa(bSGrp)-1)/2
                    else
                      dim_2 = DimSGrpa(aSGrp)*DimSGrpa(bSGrp)
                    end if
                    if (beSGrp == gaSGrp) then
                      dim_3 = DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/2
                    else
                      dim_3 = DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                    end if
                    call mc0c1at3b(dim_2,dim_3,dim_2,dim_1,dim_3,dim_1,dim_3,dim_2,dim_1,wrk(posWw),wrk(posT2w),wrk(posT2n2))

                    ! end cycle over all subgroups of (a>=b)'
                  end do
                end do

                ! Save T2n1(bega)",ij) and T2n2((bega)",ij)
                !      in recreated form, corresponding to standard
                !      storage in T2 amplitudes
                LunName = Tmp3Name(beSGrp,gaSGrp)
                if (beSGrp == gaSGrp) then
                  lent2n1 = nTri_Elem(no)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
                  lent2n2 = nTri_Elem(no-1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/2
                else
                  lent2n1 = nTri_Elem(no)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                  lent2n2 = nTri_Elem(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                end if
                T2o2v4yes(beSGrp,gaSGrp) = 1
                call SaveX(wrk(PosT2n1),lent2n1,LunAux,LunName,1,0)
                call SaveX(wrk(PosT2n2),lent2n2,LunAux,LunName,0,1)

                ! end cycle over all subgroups of (be>=ga)'
                addgapp = addgapp+DimSGrpbe(gaSGrp)
              end do
              addbepp = addbepp+DimSGrpbe(beSGrp)
            end do

            ! end cycle over all groups of (be>=ga)
            addga = addga+dimga
          end do
          addbe = addbe+dimbe
        end do

      end if
      ! end cycle over all groups of (a>=b)
      addb = addb+dimb
    end do
  end if
  adda = adda+dima
end do

return

99 format(8(i3,1x))

end subroutine o2v4
