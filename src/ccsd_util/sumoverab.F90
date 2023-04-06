!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Pavel Neogrady                                   *
!***********************************************************************

subroutine sumoverab(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3,nabstack,posabstack)
! this routine realizes summation over ab
! and calculates following contributions:
!
! T25
! T2n(ab,ij)aaaa   <- sum(e>f-aa) [ <ab||ef>aaaa . Tau(ef,ij)aaaa ]
! T2n(ab,ij)bbbb   <- sum(e>f-bb) [ <ab||ef>bbbb . Tau(ef,ij)bbbb ]
! T2n(a,b,i,j)abab <- sum(e,f-ab) [ <ab||ef>abab . Tau(ef,ij)abab ]
!
! T28
!  Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
! T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
!  Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
! T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
! T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
!    <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
!
! General Status:
!
! 1) Tau ampitudes are stored and properly mapped in V1-V3
! V1(ij,ef)   = Tau(ef,ij)aaaa
! V2(ij,ef)   = Tau(ef,ij)bbbb
! V3(i,j,e,f) = Tau(e,f,i,j)abab
!
! 2) Integrals are stored in file nab
!
! syma=1,nsym
!   symb=1,syma
!     ! one record n - in all cases
!     a = 1,nvb(syma)
!       if (syma == symb) then
!         limb = a
!       else
!         limb = nvb(symb)
!       end if
!       b=1,limb
!         ! one record N _a_b(p,q) - if any
!       end b
!     end a
!   end sumb
! end syma
!
! N.B. nab musi byt dakde deklarovane (daky common)
!
! N.B. II. in parallel.fh there are:
! nprocab - number of nodes chosen for 'sumoverab' part
! idab    - array of id's of these nodes

use ccsd_global, only: fullprint, h1, h2, h3, h4, idab, m1, m2, m3, m4, maxproc, mmul, n, nprocab, nsym, nva, nvb, r1, r2, r3, r4, &
                       r5, r6, t11, t12, t21, t22, t23, v1, v2, v3, v4
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, nabstack, posabstack
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunt2o1, lunt2o2, lunt2o3
integer(kind=iwp) :: a, aeqb, b, bb, bstart, bstop, iab, idtot(maxproc), key, left, limb, lunab, nab, nabnow, nadd, nlength, nsa, &
                     nsb, posab, post, rc, ssh1, ssh3, ssn, syma, symb, todo, yes

todo = 0
posab = posabstack

!I.parallel

!I.par.1 - escape, if this node is not reserved for sumoverab
yes = 0
do a=1,nprocab
  !if ((idab(a) == myRank) .and. (ideffab(a) > Zero)) yes = a
  if (idab(a) == myRank) yes = a
end do

if (yes == 0) return

!I.par.2 - cal overal number of ab records
nab = 0
do syma=1,nsym
  do symb=1,syma
    if (syma == symb) then
      nab = nab+(nvb(syma)*(nvb(syma)+1))/2
    else
      nab = nab+nvb(syma)*nvb(symb)
    end if
  end do
end do

!I.par.3 - distribute and calc add parameter
call sumabdistt(nab,idtot)
nab = idtot(yes)
! escape, if there is nothing to do on this node
if (nab == 0) return
nadd = sum(idtot(1:yes-1))

!I.1 open nab file

lunab = 10
call filemanager(4,lunab,rc)

!I.2  prepair
! V1(ij,ef)   = Tau(ef,ij)aaaa
! V2(ij,ef)   = Tau(ef,ij)bbbb
! V3(i,j,e,f) = Tau(e,f,i,j)abab

call filemanager(2,lunt2o1,rc)
call getmediate(wrk,wrksize,lunt2o1,v4,rc)
call mktau(wrk,wrksize,v4,t11,t12,One,rc)
call map(wrk,wrksize,4,3,4,1,2,v4,1,v1,post,rc)

call filemanager(2,lunt2o2,rc)
call getmediate(wrk,wrksize,lunt2o2,v4,rc)
call mktau(wrk,wrksize,v4,t11,t12,One,rc)
call map(wrk,wrksize,4,3,4,1,2,v4,1,v2,post,rc)

call filemanager(2,lunt2o3,rc)
call getmediate(wrk,wrksize,lunt2o3,v4,rc)
call mktau(wrk,wrksize,v4,t11,t12,One,rc)
call map(wrk,wrksize,4,3,4,1,2,v4,1,v3,post,rc)

!II. sum over ab
outer: do syma=1,nsym
  do symb=1,syma
    if (fullprint >= 2) write(u6,*) ' SymA, SymB ',syma,symb

    !II.1 read n
    call getmap(lunab,nlength,n,rc)

    !II.3 skip sum over a,b if N is empty
    if (nlength == 0) cycle

    !II.4 skip sum over a,b if # of ab is zero
    if ((nvb(syma)*nvb(symb)) == 0) cycle

    !II.5 def symmetry of N
    ssn = mmul(syma,symb)

    !II.5.1 def %d and %i for R1 _a,_b(j,e)aaaa = <ab||je>  (pos
    r1%pos0 = m1%pos0
    call grc0(2,0,1,3,0,0,ssn,post,r1)
    !II.5.2 def %d and %i for R2 _a,_b(j,e)bbbb = <ab||je>  (pos
    r2%pos0 = m2%pos0
    call grc0(2,0,2,4,0,0,ssn,post,r2)
    !II.5.3 def %d and %i for R3 _a,_b(j,e)abba = <ab||je>  (pos
    r3%pos0 = m3%pos0
    call grc0(2,0,2,3,0,0,ssn,post,r3)
    !II.5.4 def %d and %i for R4 _b,_a(j,e)abba = <ba||je>  (pos
    r4%pos0 = m4%pos0
    call grc0(2,0,2,3,0,0,ssn,post,r4)
    r5%pos0 = h1%pos0
    !II.5.5 def %d and %i for R5 _a,_b(j,e)abab = <ab||je>  (pos
    call grc0(2,0,1,4,0,0,ssn,post,r5)
    r6%pos0 = h2%pos0
    !II.5.6 def %d and %i for R6 _b,_a(j,e)abab = <ba||je>  (pos
    call grc0(2,0,1,4,0,0,ssn,post,r6)

    !II.6.1 def %d and %i for M1 _a,_b(ef)aaaa = <ab||ef>
    call grc0(2,1,3,3,0,0,ssn,post,m1)
    !II.6.2 def %d and %i for M2 _a,_b(ef)bbbb = <ab||ef>
    call grc0(2,1,4,4,0,0,ssn,post,m2)
    !II.6.3 def %d and %i for M3 _a,_b(e,f)abab = <ab||ef>
    call grc0(2,0,3,4,0,0,ssn,post,m3)
    !II.6.4 def %d and %i for M4 _b,_a(e,f)abab = <ba||fe>
    call grc0(2,0,3,4,0,0,ssn,post,m4)

    !II.7 def # of S orbitals in syma and symb
    nsa = nvb(syma)-nva(syma)
    nsb = nvb(symb)-nva(symb)

    do a=1,nvb(syma)

      if (fullprint >= 3) write(u6,*) ' A',a

      !II.7 def limitations for b
      if (syma == symb) then
        limb = a
      else
        limb = nvb(symb)
      end if

      !@@
      !par for nodes, selected for sumab process:
      !    1) only blank read (instead of _a_b(e,f) form lunab file for records,
      !    pre-standing to the portion, that is going to be done on this node
      !    2) prepair setup (def bstart,limb), if needed

      if (nadd >= limb) then
        ! all cycle over b need to be skipped
        do b=1,limb
          ! blank read and nothing more
          call reajalovy(lunab,nlength,wrk(n%pos0))
          nadd = nadd-1
        end do
        ! no setup needed
        cycle

      else if (nadd > 0) then
        ! setup first
        bstart = nadd
        if (nab < (limb-bstart)) limb = bstart+nab
        ! part of the cycle over b need to be skipped
        do b=1,nadd
          ! blank read and nothing more
          call reajalovy(lunab,nlength,wrk(n%pos0))
        end do
        nadd = 0

      else
        ! nabadd musi byt 0, inak nieje dobre
        ! no part of the cycle over b need to be skipped
        ! only setup
        bstart = 0
        if (nab < limb) limb = nab
      end if

      !* define conditions for reading to stack
      !  bstart - lower limit of b-summation
      !  bstop  - upper limit of b-summation
      !  nabnow - actual size af stack (if=1, no stacking)

      !*.1 initialization
      !Stare bstart=0

      !*.1 set unstacking conditions
      do
        bstart = bstart+1
        bstop = bstart
        nabnow = 1

        !*.2 without stacking if: 1)a=b; 2)a or b \in S; 3) nabnow <= 3
        !    (N.B. it was found inconvenient to execute stacking if
        !          size of stack is less then 4, or even more)

        ! ad 1)
        if ((bstart /= a) .or. (syma /= symb)) then
          ! ad 2)
          if ((bstart > nsb) .and. (a > nsa)) then
            ! ad 3)
            left = limb-bstart+1
            if (syma == symb) left = left-1
            if (left > 3) then

              !*.3 set stacking conditions for real stacking
              if (left <= nabstack) then
                nabnow = left
              else
                nabnow = nabstack
              end if
              bstop = bstart+nabnow-1
              todo = 1
              posab = posabstack
            end if
          end if
        end if

        !*. reading to N or to stack
        if (nabnow == 1) then
          call rea(lunab,nlength,wrk(n%pos0))
        else
          do iab=1,nabnow
            call rea(lunab,nlength,wrk(posabstack+(iab-1)*nlength))
            !if ((posabstack+iab*nlength) >= m1%pos0) then
            !end if
          end do
          !qq = wrk(posabstack)
        end if

        do b=bstart,bstop

          !III.3 read N
          !Stare call rea(lunab,nlength,wrk(n%pos0))

          !III.4 def aeqb
          if ((syma == symb) .and. (a == b)) then
            aeqb = 1
          else
            aeqb = 0
          end if

          !III.5 def key
          if ((a > nsa) .and. (b > nsb)) then
            key = 4
          else if ((a > nsa) .and. (b <= nsb)) then
            key = 2
          else if ((a <= nsa) .and. (b > nsb)) then
            key = 3
          else
            key = 1
          end if

          ! def yes for T25.1
          if (syma == symb) then
            if (a == b) then
              yes = 0
            else if (key == 4) then
              yes = 1
            else
              yes = 0
            end if
          else
            if (key == 4) then
              yes = 1
            else
              yes = 0
            end if
          end if

          ! T25
          ! T2n(ab,ij)aaaa   <- sum(e>f-aa) [ <ab||ef>aaaa . Tau(ef,ij)aaaa ]
          ! T2n(ab,ij)bbbb   <- sum(e>f-bb) [ <ab||ef>bbbb . Tau(ef,ij)bbbb ]
          ! T2n(a,b,i,j)abab <- sum(e,f-ab) [ <ab||ef>abab . Tau(ef,ij)abab ]

          if (nabnow == 1) then
            ! AB stacking is not active

            !III.6 get
            ! -> M1 _a,_b(ef)aaaa  = <ab||ef>
            ! -> M2 _a,_b(ef)bbbb  = <ab||ef>
            ! -> M3 _a,_b(e,f)abab = <ab||ef>
            ! -> M4 _b,_a(e,f)abab = <ba||fe>
            ! free: H1-H4 (V4 reserved for stacking)
            call unpackab1(wrk,wrksize,n,m1,m2,m3,m4,ssn,key,aeqb)

            if (yes == 1) then
              ! T25.1.1 H1 (ij) = V1(ij,ef) . M1 (ef)
              call ccmult(wrk,wrksize,4,2,2,2,v1,1,m1,ssn,h1,ssh1,rc)
              ! T25.1.2 T2n(ab,ij)aaaa <- 1.0 . H1(ij)
              call add(wrk,wrksize,2,4,2,5,a-nsa,b-nsb,syma,symb,One,h1,ssn,t21,1,rc)
            end if

            if (aeqb == 0) then
              ! T25.2.1 H1 (ij) = V2(ij,ef) . M2 (ef)
              call ccmult(wrk,wrksize,4,2,2,2,v2,1,m2,ssn,h1,ssh1,rc)
              ! T25.2.2 T2n(ab,ij)bbbb <- 1.0 . H1(ij)
              call add(wrk,wrksize,2,4,2,5,a,b,syma,symb,One,h1,ssn,t22,1,rc)
            end if

            if ((key == 2) .or. (key == 4)) then
              ! T25.3.1 H1 (i,j) = V3(i,j,e,f) . M3 (e,f)
              call ccmult(wrk,wrksize,4,2,2,2,v3,1,m3,ssn,h1,ssh1,rc)
              ! T25.3.2 T2n(a,b,i,j)abab <- 1.0 . H1(i,j)
              call add(wrk,wrksize,2,4,2,5,a-nsa,b,syma,symb,One,h1,ssn,t23,1,rc)
            end if

            if (((key == 3) .or. (key == 4)) .and. (aeqb == 0)) then
              ! T25.3.3 H1 (i,j) = V3(i,j,e,f) . M4 (e,f)
              call ccmult(wrk,wrksize,4,2,2,2,v3,1,m4,ssn,h1,ssh1,rc)

              ! T25.3.4 T2n(a,b,i,j)abab <- 1.0 . H1(i,j)
              call add(wrk,wrksize,2,4,2,5,b-nsb,a,symb,syma,One,h1,ssn,t23,1,rc)
            end if

          else
            ! AB stacking is active

            if (todo == 1) then
              ! once it need to be done
              todo = 0

              !III.6 get
              ! -> M1 _a,_b(ef)aaaa  = <ab||ef>
              ! -> M2 _a,_b(ef)bbbb  = <ab||ef>
              ! -> M3 _a,_b(e,f)abab = <ab||ef>
              ! -> M4 _b,_a(e,f)abab = <ba||fe>
              ! free: H1-H4 (V4 reserved for stacking)

              ! def %d and %i for V4 _a,(ef,Bp)
              call grc0stack(nabnow,1,3,3,3,0,ssn,post,v4)
              ! -> V4 _a,_b(ef)aaaa  = <ab||ef>
              call unpackab3(wrk,wrksize,n,v4,ssn,nabnow,posabstack,nlength,1)

              !T25.1.1 H2 (ij,Bp) = V1(ij,ef) . V4 (ef,Bp)
              call multstack(wrk,wrksize,v1,v4,h2,1,ssn,nabnow)

              ! def %d and %i for H1 _a,_b(ij)aa
              call grc0(2,1,1,1,0,0,ssn,post,h1)
              do bb=bstart,bstop
                !T25.1.2 ext H1(ij) <- H2(ij,_Bb) (only data transfer)
                call extstack(wrk,wrksize,h1,h2,(bb-bstart+1),nabnow)
                !T25.1.3 T2n(ab,ij)aaaa <- 1.0 . H1(ij)
                call add(wrk,wrksize,2,4,2,5,a-nsa,bb-nsb,syma,symb,One,h1,ssn,t21,1,rc)
              end do

              ! def %d and %i for V4 _a,(ef,Bp)
              call grc0stack(nabnow,1,4,4,4,0,ssn,post,v4)
              ! -> V4 _a,_b(ef)bbbb  = <ab||ef>
              call unpackab3(wrk,wrksize,n,v4,ssn,nabnow,posabstack,nlength,2)

              !T25.2.1 H2 (ij,Bp) = V2(ij,ef) . V4 (ef,Bp)
              call multstack(wrk,wrksize,v2,v4,h2,1,ssn,nabnow)

              ! def %d and %i for H1 _a,_b(ij)bb
              call grc0(2,1,2,2,0,0,ssn,post,h1)
              do bb=bstart,bstop
                !T25.2.2 ext H1(ij) <- H2(ij,Bb)
                call extstack(wrk,wrksize,h1,h2,(bb-bstart+1),nabnow)
                !T25.2 T2n(ab,ij)bbbb <- 1.0 . H1(ij)
                call add(wrk,wrksize,2,4,2,5,a,bb,syma,symb,One,h1,ssn,t22,1,rc)
              end do

              ! def %d and %i for V4 _a,(ef,Bp)
              call grc0stack(nabnow,0,3,4,4,0,ssn,post,v4)
              ! -> V4 _a,_b(e,f)abab = <ab||ef>
              call unpackab3(wrk,wrksize,n,v4,ssn,nabnow,posabstack,nlength,3)

              !T25.3.1 H2 (i,j,Bp) = V3(i,j,e,f) . V4 (e,f,Bp)
              call multstack(wrk,wrksize,v3,v4,h2,1,ssn,nabnow)

              ! def %d and %i for H1 _a,_b(ij)ab
              call grc0(2,0,1,2,0,0,ssn,post,h1)
              do bb=bstart,bstop
                !T25.3.2 ext H1(i,j) <- H2(ij,Bb)
                call extstack(wrk,wrksize,h1,h2,(bb-bstart+1),nabnow)
                !T25.3.3 T2n(a,b,i,j)abab <- 1.0 . H1(i,j)
                call add(wrk,wrksize,2,4,2,5,a-nsa,bb,syma,symb,One,h1,ssn,t23,1,rc)
              end do

              ! def %d and %i for V4 _a,(ef,Bp)
              call grc0stack(nabnow,0,3,4,3,0,ssn,post,v4)
              ! -> V4 _b,_a(e,f)abab = <ba||fe>
              call unpackab3(wrk,wrksize,n,v4,ssn,nabnow,posabstack,nlength,4)

              !T25.3.4 H2 (i,j,Bp) = V3(i,j,e,f) . V4 (e,f,Bp)
              call multstack(wrk,wrksize,v3,v4,h2,1,ssn,nabnow)

              ! def %d and %i for H1 _a,_b(ij)ab
              call grc0(2,0,1,2,0,0,ssn,post,h1)
              do bb=bstart,bstop
                !T25.3.5 ext H1(i,j) <- H2(i,j,Bb)
                call extstack(wrk,wrksize,h1,h2,(bb-bstart+1),nabnow)
                !T25.3.6 T2n(a,b,i,j)abab <- 1.0 . H1(i,j)
                call add(wrk,wrksize,2,4,2,5,bb-nsb,a,symb,syma,One,h1,ssn,t23,1,rc)
              end do

            end if
          end if

          ! T28
          !  Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
          ! T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
          !  Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
          ! T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
          ! T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
          !    <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]

          ! put appropriate data from AB_stack to N mediate, if needed
          if (nabnow > 1) then
            ! transfer data
            wrk(n%pos0:n%pos0+nlength-1) = wrk(posab:posab+nlength-1)
            ! upgrade address
            posab = posab+nlength
          end if

          !III.7 get
          ! -> M1 _a,_b(j,e)aaaa = <ab||je> (map is R1)
          ! -> M2 _a,_b(j,e)bbbb = <ab||je> (map is R2)
          ! -> M3 _a,_b(j,e)abba = <ab||je> (map is R3)
          ! -> M4 _b,_a(j,e)abba = <ba||je> (map is R4)
          ! -> H1 _a,_b(j,e)abab = <ab||je> (map is R5)
          ! -> H2 _b,_a(j,e)abab = <ba||je> (map is R6)
          ! free: H3,H4 (V4 reserved for stacking)

          call unpackab2(wrk,wrksize,n,r1,r2,r3,r4,r5,r6,ssn,key,aeqb)

          ! T28
          !  Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
          ! T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
          !  Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
          ! T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
          ! T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
          !    <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]

          ! def yes for T28.1
          if (syma == symb) then
            if (a == b) then
              yes = 0
            else if (key == 4) then
              yes = 1
            else
              yes = 0
            end if
          else
            if (key == 4) then
              yes = 1
            else
              yes = 0
            end if
          end if

          if (yes == 1) then
            !T28.1.1 H3(k,l) = M1(k,e) . T1o(e,l)aa
            call ccmult(wrk,wrksize,2,2,2,1,r1,ssn,t11,1,h3,ssh3,rc)
            !T28.1.2 H4(kl)  = H3(k,l)-H3(l,k)
            call fack(wrk,wrksize,2,1,h3,ssh3,h4,rc)
            !T28.1.3 T2n(ab,ij)aaaa <- 1.0 H4(ij)
            call add(wrk,wrksize,2,4,2,5,a-nsa,b-nsb,syma,symb,One,h4,ssn,t21,1,rc)
          end if

          if (aeqb == 0) then
            !T28.2.1 H3(k,l) = M2(k,e) . T1o(e,l)bb
            call ccmult(wrk,wrksize,2,2,2,1,r2,ssn,t12,1,h3,ssh3,rc)
            !T28.2.2 H4(kl)  = H3(k,l)-H3(l,k)
            call fack(wrk,wrksize,2,1,h3,ssh3,h4,rc)
            !T28.2.3 T2n(ab,ij)bbbb <- 1.0 H4(ij)
            call add(wrk,wrksize,2,4,2,5,a,b,syma,symb,One,h4,ssn,t22,1,rc)
          end if

          if ((key == 2) .or. (key == 4)) then
            !T28.3.1 H3(j,i)  = M3(j,e) . T1o(e,i)aa
            call ccmult(wrk,wrksize,2,2,2,1,r3,ssn,t11,1,h3,ssh3,rc)
            !T28.3.2 H4(i,j)  = H3(j,i)
            call map(wrk,wrksize,2,2,1,0,0,h3,ssn,h4,post,rc)
            !T28.3.3 T2n(a,b,i,j)abab <- -1.0 H4(i,j)
            call add(wrk,wrksize,2,4,2,5,a-nsa,b,syma,symb,-One,h4,ssn,t23,1,rc)
          end if

          if (((key == 3) .or. (key == 4)) .and. (aeqb == 0)) then
            !T28.3.4 H3(j,i)  = M4(j,e) . T1o(e,i)aa
            call ccmult(wrk,wrksize,2,2,2,1,r4,ssn,t11,1,h3,ssh3,rc)
            !T28.3.5 H4(i,j)  = H3(j,i)
            call map(wrk,wrksize,2,2,1,0,0,h3,ssn,h4,post,rc)
            !T28.3.6 T2n(a,b,i,j)abab <- -1.0 H4(i,j)
            call add(wrk,wrksize,2,4,2,5,b-nsb,a,symb,syma,-One,h4,ssn,t23,1,rc)
          end if

          if ((key == 2) .or. (key == 4)) then
            !T28.3.7 H3(i,j)  = H1(i,e) . T1o(e,j)bb
            call ccmult(wrk,wrksize,2,2,2,1,r5,ssn,t12,1,h3,ssh3,rc)
            !T28.3.8 T2n(a,b,i,j)abab <- 1.0 H3(i,j)
            call add(wrk,wrksize,2,4,2,5,a-nsa,b,syma,symb,One,h3,ssn,t23,1,rc)
          end if

          if (((key == 3) .or. (key == 4)) .and. (aeqb == 0)) then
            !T28.3.9 H3(j,i)  = H2(i,e) . T1o(e,j)bb
            call ccmult(wrk,wrksize,2,2,2,1,r6,ssn,t12,1,h3,ssh3,rc)
            !T28.3.10 T2n(a,b,i,j)abab <- 1.0 H3(i,j)
            call add(wrk,wrksize,2,4,2,5,b-nsb,a,symb,syma,One,h3,ssn,t23,1,rc)
          end if

        end do

        !par
        nab = nab-(bstop-bstart+1)
        if (nab == 0) exit outer

        if (bstop >= limb) exit
        bstart = bstop
      end do
    end do

  end do
end do outer

!IV. close lunab
call filemanager(3,lunab,rc)

return

end subroutine sumoverab
