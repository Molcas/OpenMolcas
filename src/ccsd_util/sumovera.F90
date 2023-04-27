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

subroutine sumovera(wrk,wrksize,lunt2o1,lunt2o2,lunt2o3,lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,lunw3aabb)
! This routine does:
! I) Prepair phase:
!
! II) sum over A - alpha
! FI3, WIII1, WIII2, T15, T16, T27, T2Ex
!
! FI3
!F13.1 FI(a,e)aa <- -sum(m,f-aa) [ <ma||ef>aaaa . T1o(f,m)aa ] *
!F13.2 FI(a,e)aa <- -sum(m,f-bb) [ <ma||ef>baab . T1o(f,m)bb ] *
!F13.3 FI(a,e)bb <-  sum(m,f-bb) [ <ma||fe>bbbb . T1o(f,m)bb ] +
!f13.4 FI(a,e)bb <-  sum(m,f-aa) [ <ma||fe>abab . T1o(f,m)aa ] +
!
! WIII1
!W31.1 WIII(m,e,b,j)aaaa <= <mb||ej>aaaa *
!W31.2 WIII(m,e,b,j)bbbb <= <mb||ej>bbbb +
!W31.3 WIII(m,e,b,j)aabb <= <mb||ej>abab +
!W31.4 WIII(m,e,b,j)abba <= <mb||ej>abba +
!W31.5 WIII(m,e,b,j)baab <= <mb||ej>baab *
!W31.6 WIII(m,e,b,j)bbaa <= <mb||ej>baba *
!
! WIII2
!W32.1 WIII(m,e,b,j)aaaa <- + sum(f-a) [  <mb||ef>aaaa . T1o(f,j)aa ] *
!W32.2 WIII(m,e,b,j)bbbb <- _ sum(f-b) [  <mb||ef>bbbb . T1o(f,j)bb ] +
!W32.3 WIII(m,e,b,j)aabb <- + sum(f-b) [  <mb||ef>abab . T1o(f,j)bb ] +
!W32.4 WIII(m,e,b,j)abba <- - sum(f-a) [  <mb||fe>abab . T1o(f,j)aa ] +
!W32.5 WIII(m,e,b,j)baab <-  sum(f-b) [  <mb||ef>baab . T1o(f,j)bb ] *
!W32.6 WIII(m,e,b,j)bbaa <- - sum(f-a) [  <mb||fe>baab . T1o(f,j)aa ] *
!
! T15
!T15.1 T1n(a,i)aa <- sum(n,f-aa) [ <na||fi>aaaa . T1o(f,n)aa ] *
!T15.2 T1n(a,i)aa <- sum(n,f-bb) [ <na||fi>baba . T1o(f,n)bb ] *
!T15.3 T1n(a,i)bb <- sum(n,f-bb) [ <na||fi>bbbb . T1o(f,n)bb ] +
!T15.4 T1n(a,i)bb <- sum(n,f-aa) [ <na||fi>abab . T1o(f,n)aa ] +
!
! T16
!T16.1 T1n(a,i)aa <- - sum(m,e>f-aaa) [ T2o(ef,i,m)aaaa  . <ma||ef>aaaa ] *
!T16.2 T1n(a,i)aa <- - sum(m,e,f-bab) [ T2o(e,f,i,m)abab . <ma||ef>baab ] *
!T16.3 T1n(a,i)bb <- - sum(m,e>f-bbb) [ T2o(ef,i,m)bbbb  . <ma||ef>bbbb ] +
!T16.4 T1n(a,i)bb <- + sum(m,e,f-aab) [ T2o(e,f,m,i)abab . <ma||ef>abab ] +
!
! T2Ex
!T2E.12 R(a,m,ij)aaaa   <= - sum(e>f-aa) [ <ma||ef>aaaa . Tau(ef,ij)aaaa ] *
!T2E.1 T2n(ab,ij)aaaa   <- - sum(m-a)    [ T1o(b,m)aa   .  R(a,m,ij)aaaa ] *
!T2E.2 T2n(ab,ij)aaaa         <-   sum(m-a)    [ T1o(a,m)aa   .  R(b,m,ij)aaaa ] *
!T2E.34 R(a,m,ij)bbbb   <= - sum(e>f-bb) [ <ma||ef>bbbb . Tau(ef,ij)bbbb ] +
!T2E.3 T2n(ab,ij)bbbb   <- - sum(m-b)    [ T1o(b,m)bb   .  R(a,m,ij)bbbb ] +
!T2E.4 T2n(ab,ij)bbbb         <-   sum(m-b)    [ T1o(a,m)bb   .  R(b,m,ij)bbbb ] +
!T2E.5 R1(a,m,i,j)abab <= - sum(e,f-ab) [ <ma||ef>baab . Tau(e,f,i,j)abab ] *
!T2E.5 T2n(a,b,i,j)abab <- - sum(m-b)    [ T1o(m,b)bb . R1(a,m,i,j)baab ] *
!T2E.6 R2(a,m,i,j)baab <= - sum(e,f-ab) [ <ma||ef>abab . Tau(e,f,i,j)abab ] +
!T2E.6 T2n(a,b,i,j)abab <-   sum(m-a)    [ T1o(m,a)aa . R2(a,m,i,j)abab ] +
!
! T27
!T27.1 G(b,m,i,j)aaaa  <= - sum(e-a)  [ <mb||ej>aaaa . T1o(e,i)aa ]  *
!T27.1 Q(b,m,ij)aaaa   <= G(b,m,i,j)aaaa - G(b,m,j,i)aaaa            *
!T27.1 R(c,d,ij)aaaa   <= sum(m-a)    [ T1o(c,m)aa   . Q(d,m,ij)aa ] *
!T27.1 T2n(ab,ij)aaaa   <- R(a,b,ij)aaaa - R(b,a,ij)aaaa              *
!T27.2 G(b,m,i,j)bbbb  <= - sum(e-b)  [ <mb||ej>bbbb . T1o(e,i)bb ]  +
!T27.2 Q(b,m,ij)bbbb   <= G(b,m,i,j)bbbb - G(b,m,j,i)bbbb            +
!T27.2 R(c,d,ij)bbbb   <= sum(m-b)    [ T1o(c,m)bb   . Q(d,m,ij)bb ] +
!T27.2 T2n(ab,ij)bbbb   <- R(a,b,ij)bbbb - R(b,a,ij)bbbb              +
!T27.3 G1(b,m,i,j)abab <= + sum(e-a)  [ <mb||ej>baab . T1o(e,i)aa ]  *
!T27.3 <- - sum(e-b)  [ <mb||ei>baba . T1o(e,j)bb ]  *
!T27.3 T2n(a,b,i,j)abab <- sum(m-b)    [ T1o(b,m)bb   . G1(a,m,i,j)abab ] *
!T27.4 G2(b,m,i,j)baab <= - sum(e-a)  [ <mb||ej>abab . T1o(e,i)aa ] +
!T27.4 <- + sum(e-b)  [ <mb||ei>abba . T1o(e,j)bb ] +
!T27.4 T2n(a,b,i,j)abab <- sum(m-a)    [ T1o(a,m)aa   . G2(b,m,i,j)baab ] +

! General status
!
! Temporary files lunt2o1,lunt2o2,lunt2o3
! with T2o4, T2o4b and T2o2b must be opened
!
! Integrals are stored as follows:
!
! 1) <ma||ef> for a - alpha
! lun = n1aalpha
! integrals : <ma||ef>aaaa, <ma||ef>baab
!
! do syma=1,nsym
!   ! one record with %d,%i of <ma||ef>aaaa - always
!   ! one record with %d,%i of <ma||ef>baab - always
!   do a=1,nva(syma)
!     ! one record with <ma||ef>aaaa if any
!     ! one record with <ma||ef>baab if any
!   end do
! end do
!
! 2) <ma||ef> for a - beta
! lun = n1abeta
! integrals : <ma||ef>bbbb, <ma||ef>abab
!
! do syma=1,nsym
!   ! one record with %d,%i of <ma||ef>bbbb - always
!   ! one record with %d,%i of <ma||ef>abab - always
!   do a=1,nvb(syma)
!     ! one record with <ma||ef>bbbb if any
!     ! one record with <ma||ef>abab if any
!   end do
! end do
!
! 3) <ma||ej> for a - alpha
! lun = n2aalpha
! integrals : <ma||ej>aaaa, <ma||ej>baab, <ma||ej>baba
!
! do syma=1,nsym
!   ! one record with %d,%i of <ma||ej>aaaa - always
!   ! one record with %d,%i of <ma||ej>baab - always
!   ! one record with %d,%i of <ma||ej>baba - always
!   do a=1,nva(syma)
!     ! one record with <ma||ej>aaaa if any
!     ! one record with <ma||ej>baab if any
!     ! one record with <ma||ej>baba if any
!   end do
! end do
!
! 4) <ma||ej> for a - beta
! lun = n2abeta
! integrals : <ma||ej>bbbb, <ma||ej>abab, <ma||ej>abba
!
! do syma=1,nsym
!   ! one record with %d,%i of <ma||ej>bbbb - always
!   ! one record with %d,%i of <ma||ej>abba - always
!   ! one record with %d,%i of <ma||ej>abab - always
!   do a=1,nva(syma)
!     ! one record with <ma||ej>bbbb if any
!     ! one record with <ma||ej>abba if any
!     ! one record with <ma||ej>abab if any
!   end do
! end do
!
! N.B. use and destroy : all
! N.B. # of reads      : 4
!
! Parallel status
!
! Alternatives:
!
! A) 1 proc - idaaaa,idbbbb,idaabb,idabba,idbaab,idbbaa = 1
!
! B) 2 proc - idaaaa,idbaab,idbbaa = 1
!             idbbbb,idaabb,idabba = 2
!
! C) 4 proc - idaaaa        = 1
!             idbbaa,idbaab = 2
!             idbbbb        = 3
!             idabba,idaabb = 4
!
! D) 6 proc - idaaaa = 1
!             idbaab = 2
!             idbbaa = 3
!             idbbbb = 4
!             idaabb = 5
!             idabba = 6
!
! Pilot nodes (in the case of multiplicity contribution to Tn
! is evaluated only on pivot node)
!               : idbaab for alpha, idaabb for beta
!          i.e. : A-1; B-1,2; C-2,4; D-2,5
!
! In case D) contributions F13 are calculated separately. In such
! a case on idaaaa, and idbbbb nodes contr. F13.1 and F13.3 resp.
! are set (nod add) directly to corresp. F1. On pilot nodes,
! cont. F13.2 and F13.4 are standardly added to F1's

use ccsd_global, only: f11, f12, fullprint, h1, h2, h3, h4, idaaaa, idaabb, idabba, idbaab, idbbaa, idbbbb, m1, m2, m3, m4, nsym, &
                       nva, nvb, t11, t12, t13, t14, t21, t22, t23, v1, v2, v3, v4
use Para_Info, only: MyRank
use Constants, only: One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunt2o1, lunt2o2, lunt2o3
integer(kind=iwp), intent(inout) :: lunw3aaaa, lunw3baab, lunw3bbaa, lunw3bbbb, lunw3abba, lunw3aabb
integer(kind=iwp) :: a, h1length, h2length, h3length, m1length, m2length, n1aalpha, n1abeta, n2aalpha, n2abeta, post, rc, ssh4, &
                     ssm3, ssm4, syma, yesa, yesb

!A0   parallel

!     I.par.1  - escape, if this node is not reserved for sumoverb_a

if ((myRank == idaaaa) .or. (myRank == idbaab) .or. (myRank == idbbaa)) then
  yesa = 1
else
  yesa = 0
end if

if ((myRank == idbbbb) .or. (myRank == idaabb) .or. (myRank == idabba)) then
  yesb = 1
else
  yesb = 0
end if

if ((yesa == 0) .and. (yesb == 0)) return

! ****************   A - prepair phase part 1 -------------------

!par
if (yesa == 1) then

  !A.1 read T2o4a from opened lunt2o1 file
  !    V1(ef,ij)    <- T2o(ef,ij)aaaa
  call filemanager(2,lunt2o1,rc)
  call getmediate(wrk,wrksize,lunt2o1,v1,rc)

  !A.2 V2(ij,ef) <- V1(ef,ij)
  call map(wrk,wrksize,4,3,4,1,2,v1,1,v2,post,rc)

  !A.3 mktau V1(ef,ij) <- V1(ef,ij)
  call mktau(wrk,wrksize,v1,t11,t12,One,rc)

  !A.4 expand V3(i,j,ef) <- V2(ij,ef)
  call expand(wrk,wrksize,4,5,v2,1,v3,rc)

  !A.5
  if ((idaaaa /= idbaab) .and. (myRank == idaaaa)) call set0(wrk,wrksize,f11)
end if
!parend

!par
if ((yesa == 1) .or. (yesb == 1)) then

  !A.6 read T2o2b from opened lunt2o3 file
  !    V2(e,f,i,j)  <- T2o(e,f,i,j)abab
  call filemanager(2,lunt2o3,rc)
  call getmediate(wrk,wrksize,lunt2o3,v2,rc)

  !A.6 V4(i,j,e,f) <- V2(e,f,i,j)
  call map(wrk,wrksize,4,3,4,1,2,v2,1,v4,post,rc)

  !A.7 mktau V2(e,f,i,j) <- V2(e,f,i,j)
  call mktau(wrk,wrksize,v2,t11,t12,One,rc)
end if
!parend

! now there is:
! V1(ef,ij)   = Tau(ef,ij)aaaa
! V2(e,f,i,j) = Tau(e,f,i,j)abab
! V3(i,j,ef)  = T2o(ef,ij)aaaa
! V4(i,j,e,f) = T2o(e,f,i,j)abab
! free: M1-M4, H1-H4

! ****************   B - sum over a alpha ---------------------

!par
if (yesa == 1) then

  !B.1 open n1aalpha,n2aalpha
  n1aalpha = 11
  n2aalpha = 13
  call filemanager(4,n1aalpha,rc)
  call filemanager(4,n2aalpha,rc)

  !B.2 open temp files lunw3aaaa, lunw3baab, lunw3bbaa
  !par
  if (myRank == idaaaa) call filemanager(1,lunw3aaaa,rc)

  if (myRank == idbaab) call filemanager(1,lunw3baab,rc)

  if (myRank == idbbaa) call filemanager(1,lunw3bbaa,rc)

  do syma=1,nsym

    if (fullprint >= 3) write(6,*) ' SymA alpha ',syma

    ! storing of mediates:
    !
    ! V1(ef,ij)   = Tau(ef,ij)aaaa
    ! V2(e,f,i,j) = Tau(e,f,i,j)abab
    ! V3(i,j,ef)  = T2o(ef,ij)aaaa
    ! V4(i,j,e,f) = T2o(e,f,i,j)abab
    ! M1(m,ef)  - <ma||ef>aaaa
    ! M2(m,e,f) - <ma||ef>baab
    ! H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
    ! H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
    ! H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
    ! free: M3,M4,H4

    !B.3 get %d and %i for M1,M2
    call getmap(n1aalpha,m1length,m1,rc)
    call getmap(n1aalpha,m2length,m2,rc)

    !B.4 get %d and %i for H1 - H3
    call getmap(n2aalpha,h1length,h1,rc)
    call getmap(n2aalpha,h2length,h2,rc)
    call getmap(n2aalpha,h3length,h3,rc)

    !B.5 write %d and %i of W3aaaa(H1), W3baab(H2) and W3bbaa(H3)
    !par
    if (myRank == idaaaa) call wrtmap(lunw3aaaa,h1,rc)

    if (myRank == idbaab) call wrtmap(lunw3baab,h2,rc)

    if (myRank == idbbaa) call wrtmap(lunw3bbaa,h3,rc)

    !B.6 skip cycle over a if length of all files is zero
    if ((m1length+m2length+h1length+h2length+h3length) == 0) cycle

    do a=1,nva(syma)
      if (fullprint >= 3) write(u6,*) ' A alpha ',a

      if (h1length > 0) then
        !W31.1.1 read H1(m,e,j) = <ma||ej>aaaa
        !par
        if (myRank == idaaaa) then
          call rea(n2aalpha,h1length,wrk(h1%pos0))
        else
          call reajalovy(n2aalpha,h1length,wrk(h1%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)aaaa
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(i,j,ef)  = T2o(ef,ij)aaaa
        ! V4(i,j,e,f) = T2o(e,f,i,j)abab
        ! M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
        ! H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
        ! H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab (free, but reserved map's)
        ! H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa (free, but reserved map's)
        ! free: M3,M4,H4

        !par
        if (myRank == idaaaa) then

          ! ------- cont to T15
          !T15.1 T1n(a,i)aa <- sum(n,f-aa) [ <na||fi>aaaa . T1o(f,n)aa ]

          !T15.1.1 H4(i,f,n) <- H1(n,f,i)
          call map(wrk,wrksize,3,3,2,1,0,h1,syma,h4,post,rc)

          !T15.1.2 M3(i) <- H4(i,f,n) . T1o(f,n)aa
          call ccmult(wrk,wrksize,3,2,1,2,h4,syma,t11,1,m3,ssm3,rc)

          !T15.1.3 T1n(a,i)aa <- 1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,m3,syma,t13,1,rc)

          ! ------- cont to T27
          !T27.1 G(b,m,i,j)aaaa  <= - sum(e-a)  [ <mb||ej>aaaa . T1o(e,i)aa ]
          !T27.1 Q(b,m,ij)aaaa   <= G(b,m,i,j)aaaa - G(b,m,j,i)aaaa
          !T27.1 R(c,d,ij)aaaa   <= sum(m-a)    [ T1o(c,m)aa   . Q(d,m,ij)aa ]
          !T27.1 T2n(ab,ij)aaaa   <- R(a,b,ij)aaaa - R(b,a,ij)aaaa

          !T27.1.1 M3(m,j,e) <- H1(m,e,j)
          call map(wrk,wrksize,3,1,3,2,0,h1,syma,m3,post,rc)

          !T27.1.2 H4(m,j,i) <- M3(m,j,e) . T1o(e,i)aa
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t11,1,h4,ssh4,rc)

          !T27.1.3 M3(m,ji) <- H4(m,j,i) - H4(m,i,j)
          call fack(wrk,wrksize,3,2,h4,syma,m3,rc)

          !T27.1.4 M4(c,ji) <- T1o(c,m)aa . M3(m,ji)
          call ccmult(wrk,wrksize,2,3,3,1,t11,1,m3,syma,m4,ssm4,rc)

          !T27.1.5 T2n(ab,ij) <- - (M4(a,ji)-M4(b,ji) = M4(b,ij)-M4(a,ij)
          ! since -1 was skipped in first step (T27.1)
          call add(wrk,wrksize,3,4,1,2,a,0,syma,1,One,m4,syma,t21,1,rc)

        end if
        !parend

      end if

      if (h2length > 0) then
        !W31.5 read H2(m,e,j) = <ma||ej>baab
        !par
        if (myRank == idbaab) then
          call rea(n2aalpha,h2length,wrk(h2%pos0))
        else
          call reajalovy(n2aalpha,h2length,wrk(h2%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)aaaa
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(i,j,ef)  = T2o(ef,ij)aaaa
        ! V4(i,j,e,f) = T2o(e,f,i,j)abab
        ! M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
        ! H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
        ! H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
        ! H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa (free, but reserved map's)
        ! free: M3,M4,H4

        !par
        if (myRank == idbaab) then

          ! ------- cont to T27 (part)
          !T27.3p G1(b,m,i,j)abab <= + sum(e-a)  [ <mb||ej>baab . T1o(e,i)aa ]

          !T27.3.1 M4(m,j,e) <- H2(m,e,j)
          call map(wrk,wrksize,3,1,3,2,0,h2,syma,m4,post,rc)

          !T27.3.2 M3(m,j,i) <- M4(m,j,e) . T1o(e,i)
          call ccmult(wrk,wrksize,3,2,3,1,m4,syma,t11,1,m3,ssm3,rc)

          !T27.3.3 M4(m,i,j) <- M3(m,j,i)
          call map(wrk,wrksize,3,1,3,2,0,m3,syma,m4,post,rc)

          !par
          if (idbaab /= idbbaa) then
            ! if idbaab and idbbaa are different, we need to add contribution
            ! from this part of G1 to T2n, since they will be not present on
            ! idbbaa

            !T27.3 T2n(a,b,i,j)abab <- sum(m-b)    [ T1o(b,m)bb   . G1(a,m,i,j)abab ]

            !T27.3.7 M3(b,i,j) <- T1o(b,m)bb . M4(m,i,j)
            call ccmult(wrk,wrksize,2,3,3,1,t12,1,m4,syma,m3,ssm3,rc)

            !T27.3.8 T2n(a,b,i,j) <- 1.0 . M3(b,i,j)
            call add(wrk,wrksize,3,4,1,1,a,0,syma,1,One,m3,syma,t23,1,rc)
          end if
          !parend

        end if
        !parend

      end if

      if (h3length > 0) then
        !W31.6 read H3(m,e,j) = <ma||ej>baba
        !par
        if (myRank == idbbaa) then
          call rea(n2aalpha,h3length,wrk(h3%pos0))
        else
          call reajalovy(n2aalpha,h3length,wrk(h3%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)aaaa
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(i,j,ef)  = T2o(ef,ij)aaaa
        ! V4(i,j,e,f) = T2o(e,f,i,j)abab
        ! M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
        ! M4(m,j,i) - part of G1 intermediat from previos step (only if idbaab=idbbaa)
        ! H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
        ! H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
        ! H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
        ! free: M3,H4

        !par
        if (myRank == idbbaa) then

          ! ------- cont to T27 (continue)
          !T27.3 G1(b,m,i,j)abab <- - sum(e-b)  [ <mb||ei>baba . T1o(e,j)bb ]
          !T27.3 T2n(a,b,i,j)abab <- sum(m-b)    [ T1o(b,m)bb   . G1(a,m,i,j)abab ]
          ! part of G1 is in M4(m,i,j) from part T27.3.3

          !T27.3.4 M3(m,i,e) <- H3(m,e,i)
          call map(wrk,wrksize,3,1,3,2,0,h3,syma,m3,post,rc)

          !T27.3.5 H4(m,i,j) <- M3(m,i,e) . T1o(e,j)bb
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t12,1,h4,ssh4,rc)

          !par
          if (idbaab == idbbaa) then
            ! if idbaab == idbbaa add H4 to M4 from prev step and realize
            ! contribution from both parts of G1 in one step

            !T27.3.6 (G1) M4(m,i,j) <- -H4(m,i,j)
            call add(wrk,wrksize,3,3,0,0,0,0,1,1,-One,h4,syma,m4,syma,rc)

            !T27.3.7 M3(b,i,j) <- T1o(b,m)bb . M4(m,i,j)
            call ccmult(wrk,wrksize,2,3,3,1,t12,1,m4,syma,m3,ssm3,rc)

            !T27.3.8 T2n(a,b,i,j) <- 1.0 . M3(b,i,j)
            call add(wrk,wrksize,3,4,1,1,a,0,syma,1,One,m3,syma,t23,1,rc)

          else
            ! if idbaab /= idbbaa we need to calc only contribution from
            ! 2nd part of G1

            !T27.3.7 M3(b,i,j) <- T1o(b,m)bb . H4(m,i,j)
            call ccmult(wrk,wrksize,2,3,3,1,t12,1,h4,syma,m3,ssm3,rc)

            !T27.3.8 T2n(a,b,i,j) <- -1.0 . M3(b,i,j)
            !        (minus sign, since we skip sign in previous step)
            call add(wrk,wrksize,3,4,1,1,a,0,syma,1,-One,m3,syma,t23,1,rc)
          end if
          !parend

          ! ------- cont to T15
          !T15.2 T1n(a,i)aa <- sum(n,f-bb) [ <na||fi>baba . T1o(f,n)bb ]

          !T15.2.1 H4(i,f,n) <- H3(n,f,i)
          call map(wrk,wrksize,3,3,2,1,0,h3,syma,h4,post,rc)

          !T15.2.2 M3(i) <- H4(i,f,n) . T1o(f,n)bb
          call ccmult(wrk,wrksize,3,2,1,2,h4,syma,t12,1,m3,ssm3,rc)

          !T15.2.3 T1n(a,i)aa <- 1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,m3,syma,t13,1,rc)

        end if
        !parend

      end if

      if (m1length > 0) then

        !B.7 read M1(m,ef) = <ma||ef>aaaa
        !par
        if (myRank == idaaaa) then
          call rea(n1aalpha,m1length,wrk(m1%pos0))
        else
          call reajalovy(n1aalpha,m1length,wrk(m1%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)aaaa
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(i,j,ef)  = T2o(ef,ij)aaaa
        ! V4(i,j,e,f) = T2o(e,f,i,j)abab
        ! M1(m,ef)  - <ma||ef>aaaa
        ! M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
        ! H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
        ! H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
        ! H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
        ! free: M3,M4,H4

        !par
        if (myRank == idaaaa) then

          ! ------- cont to T16
          !T16.1 T1n(a,i)aa <- - sum(m,e>f-aaa) [ T2o(ef,i,m)aaaa  . <ma||ef>aaaa ]

          !T16.1.1 M3(i) = V3(i,m,ef). M1(m,ef)
          call ccmult(wrk,wrksize,4,3,1,3,v3,1,m1,syma,m3,ssm3,rc)

          !T16.1.2 t1n(a,i)aa <- -1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,-One,m3,syma,t13,1,rc)

          ! ------- cont to T2Ex
          !T2E.12 R(a,m,ij)aaaa   <= - sum(e>f-aa) [ <ma||ef>aaaa . Tau(ef,ij)aaaa ]
          !T2E.1 T2n(ab,ij)aaaa   <- - sum(m-a)    [ T1o(b,m)aa   .  R(a,m,ij)aaaa ]
          !T2E.2 T2n(ab,ij)aaaa   <-   sum(m-a)    [ T1o(a,m)aa   .  R(b,m,ij)aaaa ]

          !T2E.12.1 M3(m,ij) = M1(m,ef) . V1(ef,ij)
          call ccmult(wrk,wrksize,3,4,3,2,m1,syma,v1,1,m3,ssm3,rc)

          !T2E.1,2.1 M4(b,ij) = T1o(b,m)aa . M3(m,ij)
          call ccmult(wrk,wrksize,2,3,3,1,t11,1,m3,ssm3,m4,ssm4,rc)

          !T2E.1,2.2 T2n(ab,ij)aaaa <- 1.0 . M4(b,ij) (@@ pozor na factor @@)
          call add(wrk,wrksize,3,4,1,1,a,0,syma,1,One,m4,syma,t21,1,rc)

          ! ------- cont to W32
          !W32.1 WIII(m,e,b,j)aaaa <- + sum(f-a) [  <mb||ef>aaaa . T1o(f,j)aa ]

          !W32.1.1 M3(m,e,f) <- M1(m,ef)
          call expand(wrk,wrksize,3,2,m1,syma,m3,rc)

          !W32.1.2 M4(m,e,j) <- M3(m,e,f) . T1o(f,j)aa
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t11,1,m4,ssm4,rc)

          !W32.1.3 H1(W3aaaa)(m,e,j) <- +1.0 . M4(m,e,j)
          call add(wrk,wrksize,3,3,0,0,0,0,1,1,One,m4,syma,h1,syma,rc)

          !W3.1 write W3aaaa(m,e,j) to  lunw3aaaa file if size is not zero
          if (h1length > 0) call wri(lunw3aaaa,h1length,wrk(h1%pos0))

          ! ------- cont to FI3
          !F13.1 FI(a,e)aa <- -sum(m,f-aa) [ <ma||ef>aaaa . T1o(f,m)aa ]

          !F13.1.1 M4(e,f,m) <- M3(m,e,f)  (M3 from previous part is not destroyed)
          call map(wrk,wrksize,3,3,1,2,0,m3,syma,m4,post,rc)

          !F13.1.2 H4(e) <- M4(e,f,m) . T1o(f,m)aa
          call ccmult(wrk,wrksize,3,2,1,2,m4,syma,t11,1,h4,ssh4,rc)

          !F13.1.3 F1(a,e)aa <- -1.0 . H4(e)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,-One,h4,syma,f11,1,rc)

        end if
        !parend

      end if

      if (m2length > 0) then
        !B.8 read M2(m,e,f) = <ma||e,f>baab
        !par
        if ((myRank == idbaab) .or. (myRank == idbbaa)) then
          call rea(n1aalpha,m2length,wrk(m2%pos0))
        else
          call reajalovy(n1aalpha,m2length,wrk(m2%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)aaaa
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(i,j,ef)  = T2o(ef,ij)aaaa
        ! V4(i,j,e,f) = T2o(e,f,i,j)abab
        ! M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>baab
        ! H1(m,e,j) - <ma||ej>aaaa  (free, but reserved map's)
        ! H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
        ! H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
        ! free: M3,M4,H4

        !par
        if (myRank == idbaab) then

          ! ------- cont to T16
          !T16.2 T1n(a,i)aa <- - sum(m,e,f-bab) [ T2o(e,f,i,m)abab . <ma||ef>baab ]

          !T16.2.1 M3(i) = V4(i,m,e,f) . M2(m,e,f)
          call ccmult(wrk,wrksize,4,3,1,3,v4,1,m2,syma,m3,ssm3,rc)

          !T16.2.2 t1n(a,i)aa <- -1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,-One,m3,syma,t13,1,rc)

          ! ------- cont to T2Ex
          !T2E.5 R1(a,m,i,j)abab <= - sum(e,f-ab) [ <ma||ef>baab . Tau(e,f,i,j)abab ]
          !T2E.5 T2n(a,b,i,j)abab <- - sum(m-b)    [ T1o(m,b)bb . R1(a,m,i,j)baab ]

          !T2E.5.1 M4(m,i,j) <- M2(m,e,f) . V2(e,f,i,j)
          call ccmult(wrk,wrksize,3,4,3,2,m2,syma,v2,1,m4,ssm4,rc)

          !T2E.5.2 M3(b,i,j) <- T1o(b,m)bb . M4(m,i,j)
          call ccmult(wrk,wrksize,2,3,3,1,t12,1,m4,syma,m3,ssm3,rc)

          !T2E.5.3 T2n(a,b,i,j)abab <- 1.0 M3(b,i,j)
          call add(wrk,wrksize,3,4,1,1,a,0,syma,1,One,m3,ssm3,t23,1,rc)
        end if
        !parend

        ! ------- cont to W32
        !W32.5 WIII(m,e,b,j)baab <-  sum(f-b) [  <mb||ef>baab . T1o(f,j)bb ]
        !W32.6 WIII(m,e,b,j)bbaa <- - sum(f-a) [  <mb||fe>baab . T1o(f,j)aa ]

        !par
        if (myRank == idbaab) then
          !W32.5.1 M4(m,e,j) = M2(m,e,f) . T1o(f,j)bb
          call ccmult(wrk,wrksize,3,2,3,1,m2,syma,t12,1,m4,ssm4,rc)

          !W32.5.2 H2(W3baab)(m,e,j) <- 1.0 M4(m,e,j)
          call add(wrk,wrksize,3,3,0,0,0,0,1,1,One,m4,syma,h2,syma,rc)

          !W3.5 write W3baab(m,e,j) to  lunw3baab file if size is not zero
          if (h2length > 0) call wri(lunw3baab,h2length,wrk(h2%pos0))
        end if
        !parend

        !par
        if (myRank == idbbaa) then
          !W32.6.1 M4(m,e,f) <- M2(m,f,e)
          call map(wrk,wrksize,3,1,3,2,0,m2,syma,m4,post,rc)

          !W32.6.2 M3(m,e,j) <- M4(m,e,f) . T1o(f,j)aa
          call ccmult(wrk,wrksize,3,2,3,1,m4,syma,t11,1,m3,ssm3,rc)

          !W32.6.3 H3(W3bbaa)(m,e,j) <- -1.0 M3(m,e,j)
          call add(wrk,wrksize,3,3,0,0,0,0,1,1,-One,m3,syma,h3,syma,rc)

          !W3.6 write W3bbaa(m,e,j) to  lunw3bbaa file if size is not zero
          if (h3length > 0) call wri(lunw3bbaa,h3length,wrk(h3%pos0))
        end if
        !parend

        ! ------- cont to FI3
        !par
        if (myRank == idbaab) then
          !F13.2 FI(a,e)aa <- -sum(m,f-bb) [ <ma||ef>baab . T1o(f,m)bb ]

          !F13.2.1 M3(e,f,m) <- M2(m,e,f)
          call map(wrk,wrksize,3,3,1,2,0,m2,syma,m3,post,rc)

          !F13.2.2 M4(e)     <- M3(e,f,m) . T1o(f,m)bb
          call ccmult(wrk,wrksize,3,2,1,2,m3,syma,t12,1,m4,ssm4,rc)

          !F13.2.3 F1(a,e)aa <- -1.0 M4(e)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,-One,m4,syma,f11,1,rc)
        end if
      end if
      !parend

    end do

  end do

  !B.9 close n1aalpha,n2aalpha
  call filemanager(3,n1aalpha,rc)
  call filemanager(3,n2aalpha,rc)
end if
!parend

! ****************   A - prepair phase part 2 -------------------

!par
if (yesb == 1) then
  !A.10 V3(j,i,e,f) <- V4(i,j,e,f)
  call map(wrk,wrksize,4,2,1,3,4,v4,1,v3,post,rc)

  !A.11 read T2o4b from opened lunt2o1 file
  !     V4(ef,ij)    <- T2o(ef,ij)aaaa
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v4,rc)

  !A.12 V1(ij,ef) <- V4(ef,ij)
  call map(wrk,wrksize,4,3,4,1,2,v4,1,v1,post,rc)

  !A.13 expand V4(i,j,ef) <- V1(ij,ef)
  call expand(wrk,wrksize,4,5,v1,1,v4,rc)

  !A.14 read T2o4b from opened lunt2o1 file
  !     V1(ef,ij)    <- T2o(ef,ij)bbbb
  call filemanager(2,lunt2o2,rc)
  call getmediate(wrk,wrksize,lunt2o2,v1,rc)

  !A.15 mktau V1(ef,ij) <- V1(ef,ij)
  call mktau(wrk,wrksize,v1,t11,t12,One,rc)
  !A.16
  if ((idbbbb /= idaabb) .and. (myRank == idbbbb)) call set0(wrk,wrksize,f12)

  ! now there is:
  ! V1(ef,ij)   = Tau(ef,ij)bbbb
  ! V2(e,f,i,j) = Tau(e,f,i,j)abab
  ! V3(j,i,e,f) = T2o(e,f,i,j)abab
  ! V4(i,j,ef)  = T2o(ef,ij)bbbb
  ! free: M1-M4, H1-H4

  ! ****************   C - sum over a beta ---------------------

  !C.1 open n1abeta,n2abeta
  n1abeta = 12
  n2abeta = 14
  call filemanager(4,n1abeta,rc)
  call filemanager(4,n2abeta,rc)

  !C.2 open temp files lunw3bbbb, lunw3abba, lunw3aabb
  !par
  if (myRank == idbbbb) call filemanager(1,lunw3bbbb,rc)

  if (myRank == idabba) call filemanager(1,lunw3abba,rc)

  if (myRank == idaabb) call filemanager(1,lunw3aabb,rc)

  do syma=1,nsym

    if (fullprint >= 3) write(u6,*) ' SymA beta ',syma

    ! storing of mediates:
    !
    ! V1(ef,ij)   = Tau(ef,ij)bbbb
    ! V2(e,f,i,j) = Tau(e,f,i,j)abab
    ! V3(j,i,e,f) = T2o(e,f,i,j)abab
    ! V4(i,j,ef)  = T2o(ef,ij)bbbb
    ! M1(m,ef)  - <ma||ef>bbbb
    ! M2(m,e,f) - <ma||ef>abab
    ! H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
    ! H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
    ! H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
    ! free: M3,M4,H4

    !C.3 get %d and %i for M1,M2
    call getmap(n1abeta,m1length,m1,rc)
    call getmap(n1abeta,m2length,m2,rc)

    !C.4 get %d and %i for H1 - H3
    call getmap(n2abeta,h1length,h1,rc)
    call getmap(n2abeta,h2length,h2,rc)
    call getmap(n2abeta,h3length,h3,rc)

    !C.5 write %d and %i of W3bbbb(H1), W3abba(H2) and W3aabb(H3)
    !par
    if (myRank == idbbbb) call wrtmap(lunw3bbbb,h1,rc)

    if (myRank == idabba) call wrtmap(lunw3abba,h2,rc)

    if (myRank == idaabb) call wrtmap(lunw3aabb,h3,rc)
    !parend

    !C.6 skip cycle over a if length of all files is zero
    if ((m1length+m2length+h1length+h2length+h3length) == 0) cycle

    do a=1,nvb(syma)
      if (fullprint >= 3) write(u6,*) ' A beta ',a

      if (h1length > 0) then
        !W31.2.1 read H1(m,e,j) = <ma||ej>bbbb
        !par
        if (myRank == idbbbb) then
          call rea(n2abeta,h1length,wrk(h1%pos0))
        else
          call reajalovy(n2abeta,h1length,wrk(h1%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)bbbb
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(j,i,e,f) = T2o(e,f,i,j)abab
        ! V4(i,j,ef)  = T2o(ef,ij)bbbb
        ! M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
        ! H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
        ! H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba (free, but reserved map's)
        ! H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb (free, but reserved map's)
        ! free: M3,M4,H4

        !par
        if (myRank == idbbbb) then

          ! ------- cont to T15
          !T15.3 T1n(a,i)bb <- sum(n,f-bb) [ <na||fi>bbbb . T1o(f,n)bb ]

          !T15.3.1 H4(i,f,n) <- H1(n,f,i)
          call map(wrk,wrksize,3,3,2,1,0,h1,syma,h4,post,rc)

          !T15.3.2 M3(i) <- H4(i,f,n) . T1o(f,n)bb
          call ccmult(wrk,wrksize,3,2,1,2,h4,syma,t12,1,m3,ssm3,rc)

          !T15.3.3 T1n(a,i)aa <- 1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,m3,syma,t14,1,rc)

          ! ------- cont to T27
          !T27.2 G(b,m,i,j)bbbb  <= - sum(e-b)  [ <mb||ej>bbbb . T1o(e,i)bb ]
          !T27.2 Q(b,m,ij)bbbb   <= G(b,m,i,j)bbbb - G(b,m,j,i)bbbb
          !T27.2 R(c,d,ij)bbbb   <= sum(m-b)    [ T1o(c,m)bb   . Q(d,m,ij)bb ]
          !T27.2 T2n(ab,ij)bbbb   <- R(a,b,ij)bbbb - R(b,a,ij)bbbb

          !T27.2.1 M3(m,j,e) <- H1(m,e,j)
          call map(wrk,wrksize,3,1,3,2,0,h1,syma,m3,post,rc)

          !T27.2.2 H4(m,j,i) <- M3(m,j,e) . T1o(e,i)bb
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t12,1,h4,ssh4,rc)

          !T27.2.3 M3(m,ji) <- H4(m,j,i) - H4(m,i,j)
          call fack(wrk,wrksize,3,2,h4,syma,m3,rc)

          !T27.2.4 M4(c,ji) <- T1o(c,m)bb . M3(m,ji)
          call ccmult(wrk,wrksize,2,3,3,1,t12,1,m3,syma,m4,ssm4,rc)

          !T27.2.5 T2n(ab,ij) <- - (M4(a,ji)-M4(b,ji) = M4(b,ij)-M4(a,ij)
          ! since -1 was skipped in first step (T27.1)
          call add(wrk,wrksize,3,4,1,2,a,0,syma,1,One,m4,syma,t22,1,rc)
        end if
        !parend

      end if

      if (h2length > 0) then
        !W31.4 read H2(m,e,j) = <ma||ej>abba
        !par
        if (myRank == idabba) then
          call rea(n2abeta,h2length,wrk(h2%pos0))
        else
          call reajalovy(n2abeta,h2length,wrk(h2%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)bbbb
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(j,i,e,f) = T2o(e,f,i,j)abab
        ! V4(i,j,ef)  = T2o(ef,ij)bbbb
        ! M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
        ! H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
        ! H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
        ! H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb (free, but reserved map's)
        ! free: M3,M4,H4

        !par
        if (myRank == idabba) then

          ! ------- cont to T27 (part)
          !T27.4p G2(b,m,i,j)baab <- + sum(e-b)  [ <mb||ei>abba . T1o(e,j)bb ] +

          !T27.4.1 M3(m,i,e) <- H2(m,e,i)
          call map(wrk,wrksize,3,1,3,2,0,h2,syma,m3,post,rc)

          !T27.4.2 M4(m,i,j) <- M3(m,i,e) . T1o(e,j)bb
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t12,1,m4,ssm4,rc)

          !par
          if (idaabb /= idabba) then
            ! if idaabb and idabba are different, we need to add contribution
            ! from this part of G2 to T2n, since they will be not present on
            ! idaabb

            !T27.4 T2n(a,b,i,j)abab <- sum(m-a)    [ T1o(a,m)aa   . G2(b,m,i,j)baab ] +

            !T27.4.7 H4(a,i,j) <- T1o(a,m)aa . M4(m,i,j)
            call ccmult(wrk,wrksize,2,3,3,1,t11,1,m4,syma,h4,ssh4,rc)

            !T27.4.8 T2n(c,a,i,j) <- 1.0 . H4(c,i,j)
            call add(wrk,wrksize,3,4,1,2,a,0,syma,1,One,h4,syma,t23,1,rc)
          end if
          !parend

        end if
        !parend

      end if

      if (h3length > 0) then
        !W31.3 read H3(m,e,j) = <ma||ej>abab
        !par
        if (myRank == idaabb) then
          call rea(n2abeta,h3length,wrk(h3%pos0))
        else
          call reajalovy(n2abeta,h3length,wrk(h3%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)bbbb
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(j,i,e,f) = T2o(e,f,i,j)abab
        ! V4(i,j,ef)  = T2o(ef,ij)bbbb
        ! M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
        ! M4(m,j,i) - part of G2 intermediat from previos step
        ! H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
        ! H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
        ! H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
        ! free: M3,H4

        !par
        if (myRank == idaabb) then

          ! ------- cont to T27 (continue)
          !T27.4 G2(b,m,i,j)baab <= - sum(e-a)  [ <mb||ej>abab . T1o(e,i)aa ] +
          !T27.4 T2n(a,b,i,j)abab <- sum(m-a)    [ T1o(a,m)aa   . G2(b,m,i,j)baab ] +
          ! N.B. part of G2(m,i,j) is in M4 from part T27.4.2

          !T27.4.3 M3(m,j,e) <- H3(m,e,j)
          call map(wrk,wrksize,3,1,3,2,0,h3,syma,m3,post,rc)

          !T27.4.4 H4(m,j,i) <- M3(m,j,e) . T1o(e,i)aa
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t11,1,h4,ssh4,rc)

          !T27.4.5 M3(m,i,j) <- H4(m,j,i)
          call map(wrk,wrksize,3,1,3,2,0,h4,syma,m3,post,rc)

          !par
          if (idaabb == idabba) then
            ! if idaabb == idabba add H4 to M4 from prev step and realize
            ! contribution from both parts of G2 in one step

            !T27.4.6 (G2) M4(m,i,j) <- - M3(m,i,j)
            call add(wrk,wrksize,3,3,0,0,0,0,1,1,-One,m3,syma,m4,syma,rc)

            !T27.4.7 H4(a,i,j) <- T1o(a,m)aa . M4(m,i,j)
            call ccmult(wrk,wrksize,2,3,3,1,t11,1,m4,syma,h4,ssh4,rc)

            !T27.4.8 T2n(c,a,i,j) <- 1.0 . H4(c,i,j)
            call add(wrk,wrksize,3,4,1,2,a,0,syma,1,One,h4,syma,t23,1,rc)

          else
            ! if idabba /= idaabb we need to calc only contribution from
            ! 2nd part of G2

            !T27.4.7 H4(a,i,j) <- T1o(a,m)aa . M3(m,i,j)
            call ccmult(wrk,wrksize,2,3,3,1,t11,1,m3,syma,h4,ssh4,rc)

            !T27.4.8 T2n(c,a,i,j) <- 1.0 . -H4(c,i,j)
            !        (minus sign, since we skip sign in previous step)
            call add(wrk,wrksize,3,4,1,2,a,0,syma,1,-One,h4,syma,t23,1,rc)
          end if
          !parend

          ! ------- cont to T15
          !T15.4 T1n(a,i)bb <- sum(n,f-aa) [ <na||fi>abab . T1o(f,n)aa ]

          !T15.4.1 H4(i,f,n) <- H3(n,f,i)
          call map(wrk,wrksize,3,3,2,1,0,h3,syma,h4,post,rc)

          !T15.4.2 M3(i) <- H4(i,f,n) . T1o(f,n)aa
          call ccmult(wrk,wrksize,3,2,1,2,h4,syma,t11,1,m3,ssm3,rc)

          !T15.4.3 T1n(a,i)bb <- 1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,m3,syma,t14,1,rc)
        end if
        !parend

      end if

      if (m1length > 0) then

        !C.7 read M1(m,ef) = <ma||ef>bbbb
        !par
        if (myRank == idbbbb) then
          call rea(n1abeta,m1length,wrk(m1%pos0))
        else
          call reajalovy(n1abeta,m1length,wrk(m1%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)bbbb
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(j,i,e,f) = T2o(e,f,i,j)abab
        ! V4(i,j,ef)  = T2o(ef,ij)bbbb
        ! M1(m,ef)  - <ma||ef>bbbb
        ! M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
        ! H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
        ! H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
        ! H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
        ! free: M3,M4,H4

        !par
        if (myRank == idbbbb) then

          ! ------- cont to T16
          !T16.3 T1n(a,i)bb <- - sum(m,e>f-bbb) [ T2o(ef,i,m)bbbb  . <ma||ef>bbbb ]

          !T16.3.1* M3(i) = V4(i,m,ef). M1(m,ef)
          call ccmult(wrk,wrksize,4,3,1,3,v4,1,m1,syma,m3,ssm3,rc)

          !T16.3.2 t1n(a,i)bb <- -1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,-One,m3,syma,t14,1,rc)

          ! ------- cont to T2Ex
          !T2E.34 R(a,m,ij)bbbb   <= - sum(e>f-bb) [ <ma||ef>bbbb . Tau(ef,ij)bbbb ] +
          !T2E.3 T2n(ab,ij)bbbb   <- - sum(m-b)    [ T1o(b,m)bb   .  R(a,m,ij)bbbb ] +
          !T2E.4 T2n(ab,ij)bbbb   <-   sum(m-b)    [ T1o(a,m)bb   .  R(b,m,ij)bbbb ] +

          !T2E.34.1 M3(m,ij) = M1(m,ef) . V1(ef,ij)
          call ccmult(wrk,wrksize,3,4,3,2,m1,syma,v1,1,m3,ssm3,rc)

          !T2E.34.1 M4(b,ij) = T1o(b,m)bb . M3(m,ij)
          call ccmult(wrk,wrksize,2,3,3,1,t12,1,m3,ssm3,m4,ssm4,rc)

          !T2E.34.2 T2n(ab,ij)bbbb <- 1.0 . M4(b,ij) (@@ pozor na factor @@)
          call add(wrk,wrksize,3,4,1,1,a,0,syma,1,One,m4,ssm4,t22,1,rc)

          !     ------- cont to W32
          !W32.2 WIII(m,e,b,j)bbbb <- + sum(f-b) [  <mb||ef>bbbb . T1o(f,j)bb ]

          !W32.2.1 M3(m,e,f) <- M1(m,ef)
          call expand(wrk,wrksize,3,2,m1,syma,m3,rc)

          !W32.2.2 M4(m,e,j) <- M3(m,e,f) . T1o(f,j)bb
          call ccmult(wrk,wrksize,3,2,3,1,m3,syma,t12,1,m4,ssm4,rc)

          !W32.2.3 H1(W3bbbb)(m,e,j) <- +1.0 . M4(m,e,j)
          call add(wrk,wrksize,3,3,0,0,0,0,1,1,One,m4,syma,h1,syma,rc)

          !W3.2 write W3bbbb(m,e,j) to  lunw3bbbb file if size is not zero
          if (h1length > 0) call wri(lunw3bbbb,h1length,wrk(h1%pos0))

          ! ------- cont to FI3
          !F13.3 FI(a,e)bb <-  sum(m,f-bb) [ <ma||fe>bbbb . T1o(f,m)bb ]

          !F13.3.1 M4(e,f,m) <- M3(m,f,e)  (M3 from previous part is not destroyed)
          call map(wrk,wrksize,3,3,2,1,0,m3,syma,m4,post,rc)

          !F13.3.2 H4(e) <- M4(e,f,m) . T1o(f,m)bb
          call ccmult(wrk,wrksize,3,2,1,2,m4,syma,t12,1,h4,ssh4,rc)

          !F13.3.3 F1(a,e)bb <- 1.0 . H4(e)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,h4,syma,f12,1,rc)
        end if
        !parend

      end if

      if (m2length > 0) then
        !B.8 read M2(m,e,f) = <ma||e,f>abab
        if ((myRank == idaabb) .or. (myRank == idabba)) then
          call rea(n1abeta,m2length,wrk(m2%pos0))
        else
          call reajalovy(n1abeta,m2length,wrk(m2%pos0))
        end if

        ! Status:
        ! V1(ef,ij)   = Tau(ef,ij)bbbb
        ! V2(e,f,i,j) = Tau(e,f,i,j)abab
        ! V3(j,i,e,f) = T2o(e,f,i,j)abab
        ! V4(i,j,ef)  = T2o(ef,ij)bbbb
        ! M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
        ! M2(m,e,f) - <ma||ef>abab
        ! H1(m,e,j) - <ma||ej>bbbb  (free, but reserved map's)
        ! H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
        ! H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
        ! free: M3,M4,H4

        !par
        if (myRank == idaabb) then

          ! ------- cont to T16
          !T16.4 T1n(a,i)bb <- + sum(m,e,f-aab) [ T2o(e,f,m,i)abab . <ma||ef>abab ]

          !T16.4.1 M3(i) = V3(i,m,e,f) . M2(m,e,f)
          call ccmult(wrk,wrksize,4,3,1,3,v3,1,m2,syma,m3,ssm3,rc)

          !T16.4.2 t1n(a,i)bb <- +1.0 . M3(i)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,m3,syma,t14,1,rc)

          ! ------- cont to T2Ex
          !T2E.6 R2(a,m,i,j)baab <= - sum(e,f-ab) [ <ma||ef>abab . Tau(e,f,i,j)abab ]
          !T2E.6 T2n(a,b,i,j)abab <-   sum(m-a)    [ T1o(m,a)aa . R2(b,m,i,j)abab ]

          !T2E.6.1 M4(m,i,j) <- M2(m,e,f) . V2(e,f,i,j)
          call ccmult(wrk,wrksize,3,4,3,2,m2,syma,v2,1,m4,ssm4,rc)

          !T2E.6.2 M3(a,i,j) <- T1o(a,m)aa . M4(m,i,j)
          call ccmult(wrk,wrksize,2,3,3,1,t11,1,m4,syma,m3,ssm3,rc)

          !T2E.6.3 T2n(a,b,i,j)abab <- -1.0 M3(a,i,j)
          call add(wrk,wrksize,3,4,1,2,a,0,syma,1,-One,m3,ssm3,t23,1,rc)
        end if
        !parend

        ! ------- cont to W32
        !W32.3 WIII(m,e,b,j)aabb <- + sum(f-b) [  <mb||ef>abab . T1o(f,j)bb ]
        !W32.4 WIII(m,e,b,j)abba <- - sum(f-a) [  <mb||fe>abab . T1o(f,j)aa ]

        !par
        if (myRank == idaabb) then
          !W32.3.1 M4(m,e,j) = M2(m,e,f) . T1o(f,j)bb
          call ccmult(wrk,wrksize,3,2,3,1,m2,syma,t12,1,m4,ssm4,rc)

          !W32.3.2 H3(W3aabb)(m,e,j) <- 1.0 M4(m,e,j)
          call add(wrk,wrksize,3,3,0,0,0,0,1,1,One,m4,syma,h3,syma,rc)

          !W3.3 write W3aabb(m,e,j) to  lunw3baab file if size is not zero
          if (h3length > 0) call wri(lunw3aabb,h3length,wrk(h3%pos0))
        end if
        !parend

        !par
        if (myRank == idabba) then
          !W32.4.1 M4(m,e,f) <- M2(m,f,e)
          call map(wrk,wrksize,3,1,3,2,0,m2,syma,m4,post,rc)

          !W32.4.2 M3(m,e,j) <- M4(m,e,f) . T1o(f,j)aa
          call ccmult(wrk,wrksize,3,2,3,1,m4,syma,t11,1,m3,ssm3,rc)

          !W32.4.3 H2(W3abba)(m,e,j) <- -1.0 M3(m,e,j)
          call add(wrk,wrksize,3,3,0,0,0,0,1,1,-One,m3,syma,h2,syma,rc)

          !W3.4 write W3abba(m,e,j) to  lunw3baba file if size is not zero
          if (h2length > 0) call wri(lunw3abba,h2length,wrk(h2%pos0))
        end if
        !parend

        ! ------- cont to FI3
        !par
        if (myRank == idaabb) then
          !F13.4 FI(a,e)bb <-  sum(m,f-aa) [ <ma||fe>abab . T1o(f,m)aa ]

          !F13.4.1 M3(e,f,m) <- M2(m,f,e)
          call map(wrk,wrksize,3,3,2,1,0,m2,syma,m3,post,rc)

          !F13.4.2 M4(e)     <- M1(e,f,m) . T1o(f,m)aa
          call ccmult(wrk,wrksize,3,2,1,2,m3,syma,t11,1,m4,ssm4,rc)

          !F13.4.3 F1(a,e)bb <- -1.0 M4(e)
          call add(wrk,wrksize,1,2,1,1,a,0,syma,1,One,m4,syma,f12,1,rc)
        end if
      end if
      !parend

    end do

  end do

  !C.9 close n1abeta,n2abeta
  call filemanager(3,n1abeta,rc)
  call filemanager(3,n2abeta,rc)
end if
!parend

return

end subroutine sumovera
