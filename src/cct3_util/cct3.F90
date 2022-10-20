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

subroutine CCT3(ireturn)
! this program calculates noniterative T3 contributions to CCSD energy

use CCT3_global, only: fullprint, ijsegkey, imax, imin, jmax, jmin, mapddp1, mapddp2, mapdfk1, mapdfk2, mapdfk3, mapdfk4, mapdfk5, &
                       mapdfk6, mapdh1, mapdh2, mapdh3, mapdl1, mapdl2, mapdm1, mapdm2, mapdm3, mapdn, mapdp, mapdr1, mapdr2, &
                       mapdr3, mapdt11, mapdt12, mapdt21, mapdt22, mapdt23, mapdv, mapdw, mapdw11, mapdw12, mapdw13, mapdw14, &
                       mapdw21, mapdw22, mapdw23, mapidp1, mapidp2, mapifk1, mapifk2, mapifk3, mapifk4, mapifk5, mapifk6, mapih1, &
                       mapih2, mapih3, mapil1, mapil2, mapim1, mapim2, mapim3, mapin, mapip, mapir1, mapir2, mapir3, mapit11, &
                       mapit12, mapit21, mapit22, mapit23, mapiv, mapiw, mapiw11, mapiw12, mapiw13, mapiw14, mapiw21, mapiw22, &
                       mapiw23, maxspace, mmul, noa, nob, nsym, posh10, posh20, posh30, posl20, posl10, posm10, posm20, posm30, &
                       posr10, posr20, posr30, posv0, posw0, symimin, symimax, symjmax, symjmin, typt3
use Para_Info, only: MyRank, nProcs
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ireturn
#include "t3int.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: counter, i, i3, id, ilow, iOff, iPrintLevel, istart, istop, j, j3, jstart, jstop, jup, k, k3, keyyes, nsg, &
                     post, rc1, ssl1, ssl2, ssm1, ssm2, symi, symi3, symij, symijk, symj, symj3, symjstart, symjstop, symk, &
                     symk3, wrksize
real(kind=wp) :: eaaa(1), eaab(1), eabb(1), ebbb(1), ec, eccsd

!par
!stare call SetTim()
!stare call MPI_COMM_RANK(MPI_COMM_WORLD,myRank,rc)
!stare call MPI_COMM_SIZE(MPI_COMM_WORLD,nProcs,rc)

!I **************  start   section **************

!I.1 read input data form INPDAT (reorg) and input file
fullprint = 0
if (iPrintLevel(-1) <= 0) fullprint = -1
if (fullprint >= 0) then
  write(u6,*)
  write(u6,*) ' **********************************'
  write(u6,*) '  Triples Contribution Calculation '
  write(u6,*) ' **********************************'
  write(u6,*)
end if
call t3reainput()

!I.2 write head to output file
if (fullprint >= 0) call t3wrhead()

!I.3.1 calc. work space requirements to fix and help mediates
call t3initfiles(wrksize)
if (fullprint >= 0) write(u6,*) ' Work space requirements :',wrksize

!I.3.2 allocate work space

call GetMem('CCT3','Max','Real',maxspace,maxspace)
maxspace = maxspace-4
if (maxspace < wrksize) then
  write(u6,*) ' Allocation of work space failed',maxspace
  write(u6,*) ' Increase the value of the variable MOLCAS_MEM'
  call Abend()
end if
call GetMem('CCT3','Allo','Real',iOff,wrksize)

!I.3.3 set wrk = 0
call cct3_mv0zero(wrksize,wrksize,Work(iOff))
if (fullprint >= 0) write(u6,*) ' Allocation of work space : Done'

!I.4 read static integrals from INTSTA (reorg) file
call t3reaintsta(Work(iOff),wrksize)

!I.5 divide fok to faa,fai,fii and dp
call cct3_divfok(Work(iOff),wrksize,mapdn,mapin,mapdp,mapip,mapdfk1,mapifk1,mapdfk2,mapifk2,mapdfk3,mapifk3,mapdfk4,mapifk4, &
                 mapdfk5,mapifk5,mapdfk6,mapifk6,mapddp1,mapidp1,mapddp2,mapidp2,rc1)

!I.6 read CCSD energy and amplitudes
call t3reaccsd(Work(iOff),wrksize,eccsd)

!I.7 get address vector T3IntPos
!    they are located at the beginning of the t3nam file
call GetIntPos()

!I.* set energies=0
eaaa = Zero
eaab = Zero
eabb = Zero
ebbb = Zero

!I.par initialize parallel counter
counter = 0

! ***** work section *****

if (ijsegkey == 0) then
  !No Segmentation
  symimin = 1
  symimax = nsym
end if

!Segmented
!noseg do symi=1,nsym
do symi=symimin,symimax
  if (noa(symi) == 0) cycle

  ! def symjstart,symjstop

  if (ijsegkey == 0) then
    !No Segmentation
    symjstart = 1
    symjstop = symi
  else

    !Segmented
    if (symimin == symjmax) then
      ! case symimin=symi=symimax
      symjstart = symjmin
      symjstop = symjmax
    else
      ! case symimin<symjmax

      if (symi == symimin) then
        symjstart = symjmin
        symjstop = symi
      else if (symi == symimax) then
        symjstart = 1
        symjstop = symjmax
      else
        ! sub case symimin<symi<symimax
        symjstart = 1
        symjstop = symi
      end if

    end if
  end if

  !noseg do symj=1,symi
  do symj=symjstart,symjstop
    if ((symi == symj) .and. (noa(symi) <= 1)) then
      cycle
    else if (noa(symj) == 0) then
      cycle
    end if

    ! define sym(ij)
    symij = mmul(symi,symj)

    ! write symmetry status of I and J
    if (fullprint > 0) write(u6,*) ' SYMI',symi,'SYMJ',symj

    if (symi == symj) then
      ilow = 2
    else
      ilow = 1
    end if

    ! loop over index I

    ! define istart,istop

    if (ijsegkey == 0) then
      !No Segmentation
      istart = ilow
      istop = noa(symi)
    else

      !Segmented
      if ((symi == symimin) .and. (symj == symjmin)) then
        istart = imin
      else
        istart = ilow
      end if

      if ((symi == symimax) .and. (symj == symjmax)) then
        istop = imax
      else
        istop = noa(symi)
      end if

    end if

    !noseg do i=ilow,noa(symi)
    do i=istart,istop

      ! get integrals <ab|ic> for given i into R1(a,bc)
      call cct3_getint(Work(iOff),wrksize,i,symi,posr10,mapdr1,mapir1,rc1)

      ! def upper limit for index j
      if (symi == symj) then
        jup = i-1
      else
        jup = noa(symj)
      end if

      ! loop over index J

      ! def jstart,jstop

      if (ijsegkey == 0) then
        !No Segmentation
        jstart = 1
        jstop = jup
      else

        !Segmented
        if ((symimin == symimax) .and. (symjmin == symjmax)) then
          ! we are in the only symmetry block taken into consideration
          if (imin == imax) then
            jstart = jmin
            jstop = jmax
          else if (i == imin) then
            jstart = jmin
            jstop = jup
          else if (i == imax) then
            jstart = 1
            jstop = jmax
          else
            jstart = 1
            jstop = jup
          end if
        else if ((symi == symimin) .and. (symj == symjmin)) then
          ! we are in initial symmetry block
          if (i == imin) then
            jstart = jmin
            jstop = jup
          else
            jstart = 1
            jstop = jup
          end if
        else if ((symi == symimax) .and. (symj == symjmax)) then
          ! we are in terminal symmetry block
          if (i == imax) then
            jstart = 1
            jstop = jmax
          else
            jstart = 1
            jstop = jup
          end if
        else
          ! we are in intermediate symmetry block
          jstart = 1
          jstop = jup
        end if

      end if

      !noseg do j=1,jup
      do j=jstart,jstop

        ! get integrals <ab|jc> for given j into R2(a,bc)
        call cct3_getint(Work(iOff),wrksize,j,symj,posr20,mapdr2,mapir2,rc1)

        do symk=1,nsym
          if (noa(symk) == 0) cycle

          if (fullprint > 1) write(u6,*) ' SYMI',symi,'SYMJ',symj,'SYMK',symk

          ! define sym(ijk)
          symijk = mmul(symij,symk)

          ! loop over index K

          do k=1,noa(symk)
            if (fullprint >= 2) then
              write(u6,999) i,j,k
999           format(' I,J,K ',3(i3,1x))
            end if

            !par update parallel counter, choose proper id for this 'portion'
            !    and skip if this portion is not for myRank
            counter = counter+1
            id = mod(counter,nProcs)
            if (myRank /= id) cycle

            ! get integrals <ab|kc> for given k into R3(a,bc)
            call cct3_getint(Work(iOff),wrksize,k,symk,posr30,mapdr3,mapir3,rc1)

            !1 ***** aaa spin combination *****

            !1* def keyyes
            if (symj > symk) then
              keyyes = 1
            else if (symj == symk) then
              if (j > k) then
                keyyes = 1
              else
                keyyes = 0
              end if
            else
              keyyes = 0
            end if

            if (keyyes == 1) then

              !1.* define maps of W(abc)
              call cct3_t3grc0(3,5,3,3,3,0,symijk,posw0,post,mapdw,mapiw)

              !1.* vanish W(abc)
              call stz(Work(iOff),wrksize,mapdw)

              !1.1 permutations (ijk) P(a,bc) (general sign+)
              nsg = 1

              !1.1.1 V graph

              !1.1.1* def L1(bc,d) <- R3(b,cd) for given k
              call defv(Work(iOff),wrksize,1,posl10,mapdl1,mapil1,ssl1,mapdr3,mapir3,symk,rc1)

              !1.1.1* ext M1(da) <- T2aaaa(da,ij) for given i,j
              call ext(Work(iOff),wrksize,4,7,i,j,symi,symj,0,mapdt21,mapit21,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !1.1.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !1.1.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !1.1.1* pack W(abc) <-  P(a,bc) [L2(bc,a)] (minus due to using Tikda instead of T2ikad)
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !1.1.2 O graph

              !1.1.2* ext L1(bc,l) <- T2aaaa(bc,kl) for given k
              call ext(Work(iOff),wrksize,4,3,k,0,symk,0,0,mapdt21,mapit21,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !1.1.2* ext M1(l,a) <- W11(l,a,ij)=<la||ij>aaaa for given ij
              call ext(Work(iOff),wrksize,4,7,i,j,symi,symj,0,mapdw11,mapiw11,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !1.1.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !1.1.2* pack W(abc) <- P(a,bc) [L2(bc,a)]
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !1.2 permutations (ikj) P(a,bc)
              nsg = -1

              !1.2.1 V graph

              !1.2.1* def L1(bc,d) <- R2(b,cd) for given j
              call defv(Work(iOff),wrksize,1,posl10,mapdl1,mapil1,ssl1,mapdr2,mapir2,symj,rc1)

              !1.2.1* ext M1(da) <- T2aaaa(da,ik) for given i,k
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdt21,mapit21,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !1.2.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !1.2.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !1.2.1* pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to using Tikda instead of T2ikad)
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !1.2.2 O graph

              !1.2.2* ext L1(bc,l) <- T2aaaa(bc,jl) for given j
              call ext(Work(iOff),wrksize,4,3,j,0,symj,0,0,mapdt21,mapit21,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !1.2.2* ext M1(l,a) <- W11(l,a,ik)=<la||ik>aaaa for given ik
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdw11,mapiw11,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !1.2.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !1.2.2* pack W(abc) <- P(a,bc) [L2(bc,a)]
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !1.3 permutations (jki) P(a,bc)
              nsg = 1

              !1.3.1 V graph

              !1.3.1* def L1(bc,d) <- R1(b,cd) for given i
              call defv(Work(iOff),wrksize,1,posl10,mapdl1,mapil1,ssl1,mapdr1,mapir1,symi,rc1)

              !1.3.1* ext M1(da) <- T2aaaa(da,jk) for given j,k
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdt21,mapit21,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !1.3.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !1.3.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !1.3.1* pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to using Tjkda instead of T2jkad)
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !1.3.2 O graph

              !1.3.2* ext L1(bc,l) <- T2aaaa(bc,kl) for given i
              call ext(Work(iOff),wrksize,4,3,i,0,symi,0,0,mapdt21,mapit21,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !1.3.2* ext M1(l,a) <- W11(l,a,jk)=<la||jk>aaaa for given jk
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdw11,mapiw11,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !1.3.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !1.3.2* pack W(abc) <-  P(a,bc) [L2(bc,a)]
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !1.4 add singles

              !1.4.0 mov V <- -W
              call cct3_t3grc0(3,5,3,3,3,0,symijk,posv0,post,mapdv,mapiv)
              call minusa(Work(iOff),wrksize,mapdw,-One)
              call setb(Work(iOff),wrksize,mapdw,mapdv,One)

              if (typt3 > 1) then
                !1.4.1 add part W2 . T1
                call t3sgl(Work(iOff),wrksize,mapdv,mapdt11,mapit11,mapdt12,mapit12,mapdw21,mapiw21,mapdw22,mapiw22,1,i,j,k,symi, &
                           symj,symk,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20, &
                           mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              if (typt3 == 3) then
                !1.4.2 add part T2 . U
                call t3sgl(Work(iOff),wrksize,mapdv,mapdfk3,mapifk3,mapdfk4,mapifk4,mapdt21,mapit21,mapdt22,mapit22,1,i,j,k,symi, &
                           symj,symk,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20, &
                           mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              !1.5  divide by denominators and calc energy contribution

              !1.5.1 divide by den.
              call t3div(Work(iOff),wrksize,mapdw,mapdv,mapddp1,mapidp1,mapddp2,mapidp2,1,i,j,k,symi,symj,symk,ec,rc1)

              !1.5.2 add energy contribution
              eaaa = eaaa+ec

            end if

            !2 ***** aab spin combination *****

            !2* def keyyes
            if (k <= nob(symk)) then
              keyyes = 1
            else
              keyyes = 0
            end if

            if (keyyes == 1) then

              !2.* define maps of W(abc)
              call cct3_t3grc0(3,1,3,3,4,0,symijk,posw0,post,mapdw,mapiw)

              !2.* vanish W(abc)
              call stz(Work(iOff),wrksize,mapdw)

              !2.1 permutations (ijk) P(a,b) (c)
              nsg = 1

              !2.1.1 V graph

              !2.1.1* def L1(b,c,d)aba <- R3(b,cd) for given k
              call defv(Work(iOff),wrksize,4,posl10,mapdl1,mapil1,ssl1,mapdr3,mapir3,symk,rc1)

              !2.1.1* ext M1(da) <- T2aaaa(da,ij) for given i,j
              call ext(Work(iOff),wrksize,4,7,i,j,symi,symj,0,mapdt21,mapit21,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.1.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !2.1.1* mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.1.1* pack W(ab,c) <- - P(a,b) [L2(b,c,a)] (minus is due to using Tijda instead of T2ijad)
              call t3addpck(Work(iOff),wrksize,3,2,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !2.1.2 O graph

              !2.1.2* ext L1(b,c,l) <- T2abab(b,c,l,k) for given k
              call ext(Work(iOff),wrksize,4,4,k,0,symk,0,0,mapdt23,mapit23,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !2.1.2* ext M1(l,a) <- W12(l,a,ij)=<la||ij>aaaa for given ij
              call ext(Work(iOff),wrksize,4,7,i,j,symi,symj,0,mapdw11,mapiw11,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.1.2* mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.1.2* pack W(ab,c) <- - P(a,b) [L2(b,c,a)] (-, premutation in V)
              call t3addpck(Work(iOff),wrksize,3,2,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !2.2 permutations (ijk) (cab) do not nontribute

              !2.3 permutations (ikj) P(a,b) (c)
              nsg = -1

              !2.3.1 V graph

              !2.3.1* def L1(b,c,d)abb <- R2(b,cd) for given j
              call defv(Work(iOff),wrksize,3,posl10,mapdl1,mapil1,ssl1,mapdr2,mapir2,symj,rc1)

              !2.3.1* ext M1(a,d) <- T2abab(a,d,i,k) for given i,k
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.3.1* map M2(d,a) <- M1(a,d)
              call cct3_map(Work(iOff),wrksize,2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,posm20,post,rc1)
              ssm2 = ssm1

              !2.3.1* mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.3.1* pack W(ab,c) <-  P(a,b) [L2(b,c,a)]
              call t3addpck(Work(iOff),wrksize,3,2,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !2.3.2 O graph

              !2.3.2* ext L1(b,c,l) <- T2abab(b,c,j,l) for given j
              call ext(Work(iOff),wrksize,4,3,j,0,symj,0,0,mapdt23,mapit23,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !2.3.2* ext M1(l,a) <- W14(l,a,ik)=<la||ik>baab for given ik
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdw14,mapiw14,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.3.2* mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.3.2* pack W(ab,c) <- P(a,b) [L2(b,c,a)]
              call t3addpck(Work(iOff),wrksize,3,2,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !2.4 permutations (ikj) (cab)
              nsg = -1

              !2.4.1 V graph

              !2.4.1* def L1(ab,d)aaa <- R2(a,bc) for given j
              call defv(Work(iOff),wrksize,1,posl10,mapdl1,mapil1,ssl1,mapdr2,mapir2,symj,rc1)

              !2.4.1* ext M1(d,c) <- T2abab(d,c,i,k) for given i,k
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.4.1* mult L2(ab,c) <- L1(ab,d) . M1(d,c)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.4.1* add (pack) W(ab,c) <-  [L2(b,c,a)] (- due to permuted T)
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(-nsg,kind=wp),mapdl2,symijk,mapdw,mapiw,symijk,rc1)

              !2.4.2 O graph

              !2.4.2* ext L1(ab,l) <- T2aaaa(ab,jl) for given j
              call ext(Work(iOff),wrksize,4,3,j,0,symj,0,0,mapdt21,mapit21,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !2.4.2* ext M1(l,c) <- W13(l,c,ik)=<lc||ik>abab for given ik
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdw13,mapiw13,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.4.2* mult L2(ab,c) <- L1(ab,l) . M1(l,c)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.4.2* add (pack) W(ab,c) <- [L2(ab,c)]
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(nsg,kind=wp),mapdl2,symijk,mapdw,mapiw,symijk,rc1)

              !2.5 permutations (jki) P(a,b) (c)
              nsg = 1

              !2.5.1 V graph

              !2.5.1* def L1(b,c,d)abb <- R1(b,cd) for given i
              call defv(Work(iOff),wrksize,3,posl10,mapdl1,mapil1,ssl1,mapdr1,mapir1,symi,rc1)

              !2.5.1* ext M1(a,d) <- T2abab(a,d,j,k) for given j,k
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.5.1* map M2(d,a) <- M1(a,d)
              call cct3_map(Work(iOff),wrksize,2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,posm20,post,rc1)
              ssm2 = ssm1

              !2.5.1* mult L2(b,c,a) <- L1(b,c,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.5.1* pack W(ab,c) <-  P(a,b) [L2(b,c,a)]
              call t3addpck(Work(iOff),wrksize,3,2,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !2.5.2 O graph

              !2.5.2* ext L1(b,c,l) <- T2abab(b,c,i,l) for given i
              call ext(Work(iOff),wrksize,4,3,i,0,symi,0,0,mapdt23,mapit23,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !2.5.2* ext M1(l,a) <- W14(l,a,jk)=<la||jk>baab for given jk
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdw14,mapiw14,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.5.2* mult L2(b,c,a) <- L1(b,c,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.5.2* pack W(ab,c) <- P(a,b) [L2(b,c,a)]
              call t3addpck(Work(iOff),wrksize,3,2,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !2.6 permutations (ikj) (cab)
              nsg = 1

              !2.6.1 V graph

              !2.6.1* def L1(ab,d)aaa <- R1(a,bc) for given i
              call defv(Work(iOff),wrksize,1,posl10,mapdl1,mapil1,ssl1,mapdr1,mapir1,symi,rc1)

              !2.6.1* ext M1(d,c) <- T2abab(d,c,j,k) for given j,k
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.6.1* mult L2(ab,c) <- L1(ab,d) . M1(d,c)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.6.1* add (pack) W(ab,c) <-  [L2(b,c,a)] (- due to permuted T)
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(-nsg,kind=wp),mapdl2,symijk,mapdw,mapiw,symijk,rc1)

              !2.6.2 O graph

              !2.6.2* ext L1(ab,l) <- T2aaaa(ab,il) for given i
              call ext(Work(iOff),wrksize,4,3,i,0,symi,0,0,mapdt21,mapit21,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !2.6.2* ext M1(l,c) <- W13(l,c,jk)=<lc||jk>abab for given jk
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdw13,mapiw13,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !2.6.2* mult L2(ab,c) <- L1(ab,l) . M1(l,c)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !2.6.2* add (pack) W(ab,c) <- [L2(ab,c)]
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(nsg,kind=wp),mapdl2,symijk,mapdw,mapiw,symijk,rc1)

              !2.7 add singles

              !2.7.0 mov V <- W
              call cct3_t3grc0(3,1,3,3,4,0,symijk,posv0,post,mapdv,mapiv)
              call minusa(Work(iOff),wrksize,mapdw,-One)
              call setb(Work(iOff),wrksize,mapdw,mapdv,One)

              if (typt3 > 1) then
                !2.7.1 add part W2 . T1
                call t3sgl(Work(iOff),wrksize,mapdv,mapdt11,mapit11,mapdt12,mapit12,mapdw21,mapiw21,mapdw23,mapiw23,2,i,j,k,symi, &
                           symj,symk,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20, &
                           mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              if (typt3 == 3) then
                !2.7.2 add part T2 . U
                call t3sgl(Work(iOff),wrksize,mapdv,mapdfk3,mapifk3,mapdfk4,mapifk4,mapdt21,mapit21,mapdt23,mapit23,2,i,j,k,symi, &
                           symj,symk,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20, &
                           mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              !2.8 divide by denominators and calc energy contribution

              !2.8.1 divide by den.
              call t3div(Work(iOff),wrksize,mapdw,mapdv,mapddp1,mapidp1,mapddp2,mapidp2,2,i,j,k,symi,symj,symk,ec,rc1)

              !2.8.2 add energy contribution
              eaab = eaab+ec

            end if

            !3 ***** abb spin combination *****

            ! Note: in spin combination 3-abb indices are changed as follows

            i3 = k
            symi3 = symk
            j3 = i
            symj3 = symi
            k3 = j
            symk3 = symj

            ! therefore, also R files are mixed:
            ! R13 = R3, R23=R1 and R33=R2, symijk3=symijk

            !3.* def keyyes
            if (symj3 > symk3) then
              if ((j3 <= nob(symj3)) .and. (k3 <= nob(symk3))) then
                keyyes = 1
              else
                keyyes = 0
              end if
            else if (symj3 == symk3) then
              if ((j3 > k3) .and. (j3 <= nob(symj3))) then
                keyyes = 1
              else
                keyyes = 0
              end if
            else
              keyyes = 0
            end if

            if (keyyes == 1) then

              !3.* define maps of W(a,bc)
              call cct3_t3grc0(3,2,3,4,4,0,symijk,posw0,post,mapdw,mapiw)

              !3.* vanish W(abc)
              call stz(Work(iOff),wrksize,mapdw)

              !3.1 permutation (ijk) (abc)
              nsg = 1

              !3.1.1 V graph

              !3.1.1* def L1(bc,d) <- R33(b,cd) for given k3
              call defv(Work(iOff),wrksize,2,posl10,mapdl1,mapil1,ssl1,mapdr2,mapir2,symk3,rc1)

              !3.1.1* ext M1(a,d) <- T2abab(a,d,i3,j3) for given i3,j3
              call ext(Work(iOff),wrksize,4,7,i3,j3,symi3,symj3,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.1.1* map M2(d,a) <- M1(a,d)
              call cct3_map(Work(iOff),wrksize,2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,posm20,post,rc1)
              ssm2 = ssm1

              !3.1.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.1.1 map L1(a,bc) <- L2(bc,a)
              call cct3_map(Work(iOff),wrksize,3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,posl10,post,rc1)
              ssl1 = ssl2
              !3.1.1* add (pack) W(a,bc) <- [L1(a,bc)]
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(nsg,kind=wp),mapdl1,symijk,mapdw,mapiw,symijk,rc1)

              !3.1.2 O graph

              !3.1.2* ext L1(bc,l) <- T2bbbb(bc,k3l) for given k3
              call ext(Work(iOff),wrksize,4,3,k3,0,symk3,0,0,mapdt22,mapit22,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !3.1.2* ext M1(l,a) <- W14(l,a,i3j3)=<la||i3j3>baab for given i3,j3
              call ext(Work(iOff),wrksize,4,7,i3,j3,symi3,symj3,0,mapdw14,mapiw14,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.1.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.1.2 map L1(a,bc) <- L2(bc,a)
              call cct3_map(Work(iOff),wrksize,3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,posl10,post,rc1)
              ssl1 = ssl2

              !3.1.2* add (pack) W(a,bc) <- [L1(a,bc)]
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(nsg,kind=wp),mapdl1,symijk,mapdw,mapiw,symijk,rc1)

              !3.2 permutations (ijk)(bac),(ijk)(cab)
              nsg = 1

              !3.2.1 V graph

              !3.2.1* def L1(a,c,d)aba <- R33(a,cd) for given k3
              call defv(Work(iOff),wrksize,4,posl10,mapdl1,mapil1,ssl1,mapdr2,mapir2,symk3,rc1)

              !3.2.1* ext M1(d,b) <- T2abab(d,b,i3,j3) for given i3,j3
              call ext(Work(iOff),wrksize,4,7,i3,j3,symi3,symj3,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.2.1* mult L2(ac,b) <- L1(a,c,d) . M1(d,b)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.2.1* add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
              call t3addpck(Work(iOff),wrksize,3,3,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !3.2.2 O graph

              !3.2.2* ext L1(a,c,l) <- T2abab(a,c,l,k3) for given k3
              call ext(Work(iOff),wrksize,4,4,k3,0,symk3,0,0,mapdt23,mapit23,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !3.2.2* ext M1(l,b) <- W13(l,b,i3j3)=<la||i3j3>abab for given i3,j3
              call ext(Work(iOff),wrksize,4,7,i3,j3,symi3,symj3,0,mapdw13,mapiw13,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.2.2* mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.2.2* add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm V)
              call t3addpck(Work(iOff),wrksize,3,3,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !3.3 permutation (ikj) (abc)
              nsg = -1

              !3.3.1 V graph

              !3.3.1* def L1(bc,d) <- R23(b,cd) for given j3
              call defv(Work(iOff),wrksize,2,posl10,mapdl1,mapil1,ssl1,mapdr1,mapir1,symj3,rc1)

              !3.3.1* ext M1(a,d) <- T2abab(a,d,i3,k3) for given i3,k3
              call ext(Work(iOff),wrksize,4,7,i3,k3,symi3,symk3,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.3.1* map M2(d,a) <- M1(a,d)
              call cct3_map(Work(iOff),wrksize,2,2,1,0,0,mapdm1,mapim1,ssm1,mapdm2,mapim2,posm20,post,rc1)
              ssm2 = ssm1

              !3.3.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.3.1 map L1(a,bc) <- L2(bc,a)
              call cct3_map(Work(iOff),wrksize,3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,posl10,post,rc1)
              ssl1 = ssl2

              !3.3.1* add (pack) W(a,bc) <- [L1(a,bc)]
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(nsg,kind=wp),mapdl1,symijk,mapdw,mapiw,symijk,rc1)

              !3.3.2 O graph

              !3.3.2* ext L1(bc,l) <- T2bbbb(bc,j3l) for given j3
              call ext(Work(iOff),wrksize,4,3,j3,0,symj3,0,0,mapdt22,mapit22,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !3.3.2* ext M1(l,a) <- W14(l,a,i3k3)=<la||i3k3>baab for given i3,k3
              call ext(Work(iOff),wrksize,4,7,i3,k3,symi3,symk3,0,mapdw14,mapiw14,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.3.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.3.2map L1(a,bc) <- L2(bc,a)
              call cct3_map(Work(iOff),wrksize,3,2,3,1,0,mapdl2,mapil2,ssl2,mapdl1,mapil1,posl10,post,rc1)
              ssl1 = ssl2

              !3.3.2* add (pack) W(a,bc) <- [L1(a,bc)]
              call cct3_add(Work(iOff),wrksize,3,3,0,0,0,0,1,1,real(nsg,kind=wp),mapdl1,symijk,mapdw,mapiw,symijk,rc1)

              !3.4 permutations (ikj)(bac),(ikj)(cab)
              nsg = -1

              !3.4.1 V graph

              !3.4.1* def L1(a,c,d)aba <- R23(a,cd) for given j3
              call defv(Work(iOff),wrksize,4,posl10,mapdl1,mapil1,ssl1,mapdr1,mapir1,symj3,rc1)

              !3.4.1* ext M1(d,b) <- T2abab(d,b,i3,k3) for given i3,k3
              call ext(Work(iOff),wrksize,4,7,i3,k3,symi3,symk3,0,mapdt23,mapit23,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.4.1* mult L2(ac,b) <- L1(a,c,d) . M1(d,b)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.4.1* add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
              call t3addpck(Work(iOff),wrksize,3,3,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !3.4.2 O graph

              !3.4.2* ext L1(a,c,l) <- T2abab(a,c,l,j3) for given j3
              call ext(Work(iOff),wrksize,4,4,j3,0,symj3,0,0,mapdt23,mapit23,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !3.4.2* ext M1(l,b) <- W13(l,b,i3k3)=<la||i3k3>abab for given i3,k3
              call ext(Work(iOff),wrksize,4,7,i3,k3,symi3,symk3,0,mapdw13,mapiw13,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.4.2* mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.4.2* add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm V)
              call t3addpck(Work(iOff),wrksize,3,3,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !3.5 permutations (jki)(abc) do not contribute

              !3.6 permutations (jki)(bac),(jki)(cab)
              nsg = 1

              !3.6.1 V graph

              !3.6.1* def L1(a,c,d)abb <- R13(a,cd) for given i3
              call defv(Work(iOff),wrksize,3,posl10,mapdl1,mapil1,ssl1,mapdr3,mapir3,symi3,rc1)

              !3.6.1* ext M1(db) <- T2bbbb(db,j3k3) for given j3,k3
              call ext(Work(iOff),wrksize,4,7,j3,k3,symj3,symk3,0,mapdt22,mapit22,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.6.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !3.6.1* mult L2(ac,b) <- L1(a,c,d) . M2(d,b)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.6.1* add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)] (- do to perm T)
              call t3addpck(Work(iOff),wrksize,3,3,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !3.6.2 O graph

              !3.6.2* ext L1(a,c,l) <- T2abab(a,c,i3,l) for given i3
              call ext(Work(iOff),wrksize,4,3,i3,0,symi3,0,0,mapdt23,mapit23,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !3.6.2* ext M1(l,b) <- W12(l,b,j3k3)=<la||j3k3>bbbb for given j3,k3
              call ext(Work(iOff),wrksize,4,7,j3,k3,symj3,symk3,0,mapdw12,mapiw12,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !3.6.2* mult L2(ac,b) <- L1(a,c,l) . M1(l,b)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !3.6.2* add (pack) W(a,bc) <- - P(a,c) [L2(ac,b)]
              call t3addpck(Work(iOff),wrksize,3,3,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !3.7 add singles

              !3.7.0 mov V <- W
              call cct3_t3grc0(3,2,3,4,4,0,symijk,posv0,post,mapdv,mapiv)
              call minusa(Work(iOff),wrksize,mapdw,-One)
              call setb(Work(iOff),wrksize,mapdw,mapdv,One)

              if (typt3 > 1) then
                !3.7.1 add part W2 . T1
                call t3sgl(Work(iOff),wrksize,mapdv,mapdt11,mapit11,mapdt12,mapit12,mapdw23,mapiw23,mapdw22,mapiw22,3,i3,j3,k3, &
                           symi3,symj3,symk3,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2, &
                           posh20,mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              if (typt3 == 3) then
                !3.7.2 add part T2 . U
                call t3sgl(Work(iOff),wrksize,mapdv,mapdfk3,mapifk3,mapdfk4,mapifk4,mapdt23,mapit23,mapdt22,mapit22,3,i3,j3,k3, &
                           symi3,symj3,symk3,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2, &
                           posh20,mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              !3.8 divide by denominators and calc energy contribution

              !3.8.1 divide by den.
              call t3div(Work(iOff),wrksize,mapdw,mapdv,mapddp1,mapidp1,mapddp2,mapidp2,3,i3,j3,k3,symi3,symj3,symk3,ec,rc1)

              !3.8.2 add energy contribution
              eabb = eabb+ec

            end if

            !4 ***** bbb spin combination *****

            !4.* def keyyes
            if (symj > symk) then
              keyyes = 1
            else if (symj == symk) then
              if (j > k) then
                keyyes = 1
              else
                keyyes = 0
              end if
            else
              keyyes = 0
            end if

            if ((i > nob(symi)) .or. (j > nob(symj)) .or. (k > nob(symk))) then
              keyyes = 0
            end if

            if (keyyes == 1) then

              !4.* define maps of W(abc)
              call cct3_t3grc0(3,5,4,4,4,0,symijk,posw0,post,mapdw,mapiw)

              !4.* vanish W(abc)
              call stz(Work(iOff),wrksize,mapdw)

              !4.1 permutations (ijk) P(a,bc) (general sign+)
              nsg = 1

              !4.1.1 V graph

              !4.1.1* def L1(bc,d) <- R3(b,cd) for given k
              call defv(Work(iOff),wrksize,2,posl10,mapdl1,mapil1,ssl1,mapdr3,mapir3,symk,rc1)

              !4.1.1* ext M1(da) <- T2bbbb(da,ij) for given i,j
              call ext(Work(iOff),wrksize,4,7,i,j,symi,symj,0,mapdt22,mapit22,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !4.1.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !4.1.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !4.1.1* pack W(abc) <-  P(a,bc) [L2(bc,a)] (minus due to using Tikda instead of T2ikad)
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !4.1.2 O graph

              !4.1.2* ext L1(bc,l) <- T2bbbb(bc,kl) for given k
              call ext(Work(iOff),wrksize,4,3,k,0,symk,0,0,mapdt22,mapit22,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !4.1.2* ext M1(l,a) <- W12(l,a,ij)=<la||ij>bbbb for given ij
              call ext(Work(iOff),wrksize,4,7,i,j,symi,symj,0,mapdw12,mapiw12,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !4.1.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !4.1.2* pack W(abc) <- P(a,bc) [L2(bc,a)]
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !4.2 permutations (ikj) P(a,bc)
              nsg = -1

              !4.2.1 V graph

              !4.2.1* def L1(bc,d) <- R2(b,cd) for given j
              call defv(Work(iOff),wrksize,2,posl10,mapdl1,mapil1,ssl1,mapdr2,mapir2,symj,rc1)

              !4.2.1* ext M1(da) <- T2bbbb(da,ik) for given i,k
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdt22,mapit22,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !4.2.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !4.2.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !4.2.1* pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to using Tikda instead of T2ikad)
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !4.2.2 O graph

              !4.2.2* ext L1(bc,l) <- T2bbbb(bc,jl) for given j
              call ext(Work(iOff),wrksize,4,3,j,0,symj,0,0,mapdt22,mapit22,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !4.2.2* ext M1(l,a) <- W12(l,a,ik)=<la||ik>bbbb for given ik
              call ext(Work(iOff),wrksize,4,7,i,k,symi,symk,0,mapdw12,mapiw12,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !4.2.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !4.2.2* pack W(abc) <- P(a,bc) [L2(bc,a)]
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !4.3 permutations (jki) P(a,bc)
              nsg = 1

              !4.3.1 V graph

              !4.3.1* def L1(bc,d) <- R1(b,cd) for given i
              call defv(Work(iOff),wrksize,2,posl10,mapdl1,mapil1,ssl1,mapdr1,mapir1,symi,rc1)

              !4.3.1* ext M1(da) <- T2bbbb(da,jk) for given j,k
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdt22,mapit22,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !4.3.1* exp M2(d,a) <- M1(da)
              call cct3_expand(Work(iOff),wrksize,2,1,mapdm1,ssm1,posm20,mapdm2,mapim2,rc1)
              ssm2 = ssm1

              !4.3.1* mult L2(bc,a) <- L1(bc,d) . M2(d,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm2,mapim2,ssm2,mapdl2,mapil2,ssl2,posl20,rc1)

              !4.3.1* pack W(abc) <- - P(a,bc) [L2(bc,a)] (minus is due to using Tjkda instead of T2jkad)
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,-nsg,0,rc1)

              !4.3.2 O graph

              !4.3.2* ext L1(bc,l) <- T2bbbb(bc,kl) for given i
              call ext(Work(iOff),wrksize,4,3,i,0,symi,0,0,mapdt22,mapit22,1,posl10,mapdl1,mapil1,ssl1,rc1)

              !4.3.2* ext M1(l,a) <- W12(l,a,jk)=<la||jk>bbbb for given jk
              call ext(Work(iOff),wrksize,4,7,j,k,symj,symk,0,mapdw12,mapiw12,1,posm10,mapdm1,mapim1,ssm1,rc1)

              !4.3.2* mult L2(bc,a) <- L1(bc,l) . M1(l,a)
              call cct3_mult(Work(iOff),wrksize,3,2,3,1,mapdl1,mapil1,ssl1,mapdm1,mapim1,ssm1,mapdl2,mapil2,ssl2,posl20,rc1)

              !4.3.2* pack W(abc) <-  P(a,bc) [L2(bc,a)]
              call t3addpck(Work(iOff),wrksize,3,1,mapdl2,mapil2,mapdw,nsg,0,rc1)

              !4.4 add singles

              !4.4.0 mov V <- W
              call cct3_t3grc0(3,5,4,4,4,0,symijk,posv0,post,mapdv,mapiv)
              call minusa(Work(iOff),wrksize,mapdw,-One)
              call setb(Work(iOff),wrksize,mapdw,mapdv,One)

              if (typt3 > 1) then
                !4.4.1 add part W2 . T1
                call t3sgl(Work(iOff),wrksize,mapdv,mapdt12,mapit12,mapdt11,mapit11,mapdw22,mapiw22,mapdw21,mapiw21,1,i,j,k,symi, &
                           symj,symk,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20, &
                           mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              if (typt3 == 3) then
                !4.4.2 add part T2 . U
                call t3sgl(Work(iOff),wrksize,mapdv,mapdfk4,mapifk4,mapdfk3,mapifk3,mapdt22,mapit22,mapdt21,mapit21,1,i,j,k,symi, &
                           symj,symk,rc1,mapdm1,mapim1,posm10,mapdh1,mapih1,posh10,mapdm2,mapim2,posm20,mapdh2,mapih2,posh20, &
                           mapdm3,mapim3,posm30,mapdh3,mapih3,posh30)
              end if

              !4.5 divide by denominators and calc energy contribution

              !4.5.1 divide by den.
              call t3div(Work(iOff),wrksize,mapdw,mapdv,mapddp2,mapidp2,mapddp1,mapidp1,1,i,j,k,symi,symj,symk,ec,rc1)

              !4.5.2 add energy contribution
              ebbb = ebbb+ec

            end if

          end do
        end do

        !par Separate printing of partial energies e... are
        !    useful only in serial run. For parallel run also
        !    cycle over k is segmented (via parallelization),
        !    so e... are not complete contributions (only sum
        !    over all nodes have some sense). Thus, these values
        !    in parallel run are too dangerous to use separately
        !    so their printout is suppressed.

        if (nProcs == 1) then
          call t3wresult(symi,symj,i,j,eaaa(1),eaab(1),eabb(1),ebbb(1))
          if (fullprint > 1) then
            write(u6,*) ' Eaaa =',eaaa
            write(u6,*) ' Eaab =',eaab
            write(u6,*) ' Eabb =',eabb
            write(u6,*) ' Ebbb =',ebbb
          end if
        end if
        !endpar

      end do
    end do
  end do
end do

!o ***** final section *****

!o.* allreduced energy components
call gadgop(eaaa,1,'+')
call gadgop(eaab,1,'+')
call gadgop(eabb,1,'+')
call gadgop(ebbb,1,'+')
!stare call MPI_ALLREDUCE(ebbb,ec,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
!      ebbb = ec

!o.* type results

if (fullprint >= 0) then
  write(u6,'(6X,A,F24.13)') ' CCSD     =',eccsd
  write(u6,'(6X,A,F24.13)') ' T3 corr. =',eaaa+eaab+eabb+ebbb
  write(u6,'(6X,A,F24.13)') ' CCSD + T3=',eccsd+eaaa+eaab+eabb+ebbb
  write(u6,*) ' T3 energy decomposition into spin parts'
  write(u6,*) ' Eaaa =',eaaa
  write(u6,*) ' Eaab =',eaab
  write(u6,*) ' Eabb =',eabb
  write(u6,*) ' Ebbb =',ebbb
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') 'Happy Landing!'
  write(u6,*)
end if
call add_Info('E_CCSD_T',eccsd+eaaa+eaab+eabb+ebbb,1,8)
! Export a method and energy to the MOLCAS runfile
call Put_cArray('Relax Method','CCSDT   ',8)
call Store_Energies(1,eccsd+eaaa+eaab+eabb+ebbb,1)
! Releasing the memory
call GetMem('CCT3','Free','Real',iOff,wrksize)

ireturn = 0

return

end subroutine CCT3
