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

subroutine defv(wrk,wrksize,deftyp,v,ssv,r,ssr,rc)
! this routine defines v mediate as:
! V_i(abc) = <ab||ic> from integrals, stored in R_i(abc)
! R _i(a,bc)bbb is a matrix, where integrals <ab|ic> are stored
! for b>=c
!
! deftyp - type of definition (see table) (I)
! v      - V (I/O)
! ssv    - overall spin of V (O)
! r      - R (I)
! ssr    - overall spin of R (I)
! rc     - return (error) code (O)
!
! Table of definitions
!
! deftyp         Operation                 Implement.
! 1      V(ab,c)aaa  = R(abc)-R(bac)         Yes
! 2      V(ab,c)bbb  = R(abc)-R(bac)         Yes
! 3      V(a,b,c)abb = R(abc)                Yes
! 4      V(a,b,c)aba = -R(bac)               Yes

use CCT3_global, only: dimm, Map_Type, nvb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, deftyp, ssr
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(inout) :: v
integer(kind=iwp), intent(out) :: ssv
type(Map_Type), intent(in) :: r
integer(kind=iwp), intent(inout) :: rc
integer(kind=iwp) :: ir1, ir2, iv, nhelp1, nhelp10, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6, nhelp7, nhelp8, nhelp9, posr1, posr2, &
                     post, posv, syma, symb, symc

!0.* def v%d,v%i
if (deftyp == 1) then
  ! case V(ab,c)aaa
  call cct3_grc0(3,1,3,3,3,0,ssr,v,post)
else if (deftyp == 2) then
  ! case V(ab,c)bbb
  call cct3_grc0(3,1,4,4,4,0,ssr,v,post)
else if (deftyp == 3) then
  ! case V(a,b,c)abb
  call cct3_grc0(3,0,3,4,4,0,ssr,v,post)
else if (deftyp == 4) then
  ! case V(a,b,c)aba
  call cct3_grc0(3,0,3,4,3,0,ssr,v,post)
else
  ! RC=1 , deftyp out of range (1-4) (Stup)
  rc = 1
  return
end if

!0.* some tests

!0.* define spin
ssv = ssr

if ((deftyp == 1) .or. (deftyp == 2)) then

  !12 case V(ab,c)aaa,bbb

  do iv=1,v%d(0,5)

    !12.* def symmetries
    syma = v%d(iv,3)
    symb = v%d(iv,4)
    symc = v%d(iv,5)

    if (syma == symb) then
      !12.1 syma=symb

      if (symb == symc) then
        !12.1.1 case syma=symb=symc

        !12.1.1.* def positions of V, R1
        posv = v%d(iv,1)
        ir1 = r%i(syma,symb,1)
        posr1 = r%d(ir1,1)

        !12.1.1.* def dimensions
        nhelp1 = nvb(syma)
        nhelp2 = nhelp1*(nhelp1+1)/2
        nhelp4 = dimm(v%d(0,1),syma)
        nhelp3 = nhelp4*(nhelp4-1)/2
        nhelp5 = nvb(syma)-nhelp4

        !12.1.1.* realize definition of V
        call defvhlp1(wrk(posr1),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5)

      else
        !12.1.2 case syma=symb/=symc

        !12.1.2.* def positions of V, R1
        posv = v%d(iv,1)
        if (syma >= symc) then
          ir1 = r%i(syma,symb,1)
        else
          ir1 = r%i(syma,symc,1)
        end if
        posr1 = r%d(ir1,1)

        !12.1.2.* def dimensions
        nhelp1 = nvb(syma)
        nhelp2 = nvb(symc)
        nhelp4 = dimm(v%d(0,1),syma)
        nhelp3 = nhelp4*(nhelp4-1)/2
        nhelp5 = dimm(v%d(0,3),symc)
        nhelp6 = nvb(syma)-nhelp4
        nhelp7 = nvb(symc)-nhelp5

        !12.1.2.* realize definition of V
        if (syma >= symc) then
          call defvhlp21(wrk(posr1),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
        else
          call defvhlp22(wrk(posr1),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
        end if

      end if

    else
      !12.2 syma > symb

      if (syma == symc) then
        !12.2.1 case syma>symb , syma=symc

        !12.2.1.* def positions of V, R1, R2
        posv = v%d(iv,1)
        ! R1 is permuted, since (a=c)>b
        ir1 = r%i(syma,symc,1)
        posr1 = r%d(ir1,1)
        ir2 = r%i(symb,syma,1)
        posr2 = r%d(ir2,1)

        !12.2.1.* def dimensions
        nhelp1 = nvb(syma)
        nhelp2 = nvb(symb)
        nhelp3 = nvb(symc)
        nhelp4 = nhelp1*(nhelp1+1)/2
        nhelp5 = dimm(v%d(0,1),syma)
        nhelp6 = dimm(v%d(0,2),symb)
        nhelp7 = dimm(v%d(0,3),symc)
        nhelp8 = nvb(syma)-nhelp5
        nhelp9 = nvb(symb)-nhelp6
        nhelp10 = nvb(symc)-nhelp7

        !12.2.1.* realize definition of V
        call defvhlp3(wrk(posr1),wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,nhelp9,nhelp10)

      else if (symb == symc) then
        !12.2.2 case syma>symb , symb=symc

        !12.2.2.* def positions of V, R1, R2
        posv = v%d(iv,1)
        ir1 = r%i(syma,symb,1)
        posr1 = r%d(ir1,1)
        ir2 = r%i(symb,syma,1)
        posr2 = r%d(ir2,1)

        !12.2.2.* def dimensions
        nhelp1 = nvb(syma)
        nhelp3 = nvb(symb)
        nhelp4 = nvb(symc)
        nhelp2 = nhelp3*(nhelp3+1)/2
        nhelp5 = dimm(v%d(0,1),syma)
        nhelp6 = dimm(v%d(0,2),symb)
        nhelp7 = dimm(v%d(0,3),symc)
        nhelp8 = nvb(syma)-nhelp5
        nhelp9 = nvb(symb)-nhelp6
        nhelp10 = nvb(symc)-nhelp7

        !12.2.2.* realize definition of V
        call defvhlp4(wrk(posr1),wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,nhelp9,nhelp10)

      else
        !12.2.3 case syma>symb , symc/=syma,symb

        !12.2.3.* def positions of V, R1, R2
        posv = v%d(iv,1)

        if (symb >= symc) then
          ir1 = r%i(syma,symb,1)
        else
          ir1 = r%i(syma,symc,1)
        end if
        posr1 = r%d(ir1,1)

        if (syma >= symc) then
          ir2 = r%i(symb,syma,1)
        else
          ir2 = r%i(symb,symc,1)
        end if
        posr2 = r%d(ir2,1)

        !12.2.3.* def dimensions
        nhelp1 = nvb(syma)
        nhelp2 = nvb(symb)
        nhelp3 = nvb(symc)
        nhelp4 = dimm(v%d(0,1),syma)
        nhelp5 = dimm(v%d(0,2),symb)
        nhelp6 = dimm(v%d(0,3),symc)
        nhelp7 = nvb(syma)-nhelp4
        nhelp8 = nvb(symb)-nhelp5
        nhelp9 = nvb(symc)-nhelp6

        !12.2.3.* realize definition of V
        if ((syma > symc) .and. (symb > symc)) then
          call defvhlp51(wrk(posr1),wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,nhelp9)
        else if ((syma > symc) .and. (symb < symc)) then
          call defvhlp52(wrk(posr1),wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,nhelp9)
        else if ((syma < symc) .and. (symb > symc)) then
          call defvhlp53(wrk(posr1),wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,nhelp9)
        else if ((syma < symc) .and. (symb < symc)) then
          call defvhlp54(wrk(posr1),wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8,nhelp9)
        end if

      end if

    end if

  end do

else if (deftyp == 3) then

  !3 case V(a,b,c)abb

  do iv=1,v%d(0,5)

    !3.* def symmetries
    syma = v%d(iv,3)
    symb = v%d(iv,4)
    symc = v%d(iv,5)

    if (symb == symc) then
      !3.1 symb=symc

      !3.1.* def positions of V, R1
      posv = v%d(iv,1)
      ir1 = r%i(syma,symb,1)
      posr1 = r%d(ir1,1)

      !3.1.* def dimensions
      nhelp1 = nvb(syma)
      nhelp2 = nvb(symb)
      nhelp3 = nhelp2*(nhelp2+1)/2
      nhelp4 = dimm(v%d(0,1),syma)
      nhelp5 = dimm(v%d(0,2),symb)
      nhelp6 = dimm(v%d(0,3),symc)
      nhelp7 = nvb(syma)-nhelp4

      !3.1.* realize definition of V
      call defvhlp7(wrk(posr1),wrk(posv),nhelp1,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)

    else
      ! case symb/=symc

      !3.1.* def positions of V, R1
      posv = v%d(iv,1)
      if (symb >= symc) then
        ir1 = r%i(syma,symb,1)
      else
        ir1 = r%i(syma,symc,1)
      end if
      posr1 = r%d(ir1,1)

      !3.1.* def dimensions
      nhelp1 = nvb(syma)
      nhelp2 = nvb(symb)
      nhelp3 = nvb(symc)
      nhelp4 = dimm(v%d(0,1),syma)
      nhelp5 = dimm(v%d(0,2),symb)
      nhelp6 = dimm(v%d(0,3),symc)
      nhelp7 = nvb(syma)-nhelp4

      !3.1.* realize definition of V
      if (symb > symc) then
        call defvhlp61(wrk(posr1),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
      else
        call defvhlp62(wrk(posr1),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7)
      end if

    end if

  end do

else if (deftyp == 4) then

  !4 case V(a,b,c)aba

  do iv=1,v%d(0,5)

    !4.* def symmetries
    syma = v%d(iv,3)
    symb = v%d(iv,4)
    symc = v%d(iv,5)

    if (syma == symc) then
      !4.1 syma=symc

      !4.1.* def positions of V, R2
      posv = v%d(iv,1)
      ir2 = r%i(symb,syma,1)
      posr2 = r%d(ir2,1)

      !4.1.* def dimensions
      nhelp1 = nvb(symb)
      nhelp2 = nvb(syma)
      nhelp3 = nhelp2*(nhelp2+1)/2
      nhelp4 = dimm(v%d(0,1),syma)
      nhelp5 = dimm(v%d(0,2),symb)
      nhelp6 = dimm(v%d(0,3),symc)
      nhelp7 = nvb(syma)-nhelp4
      nhelp8 = nvb(symc)-nhelp6

      !4.1.* realize definition of V
      call defvhlp9(wrk(posr2),wrk(posv),nhelp1,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)

    else
      ! case symb/=symc

      !4.1.* def positions of V, R2
      posv = v%d(iv,1)
      if (syma >= symc) then
        ir2 = r%i(symb,syma,1)
      else
        ir2 = r%i(symb,symc,1)
      end if
      posr2 = r%d(ir2,1)

      !4.1.* def dimensions
      nhelp1 = nvb(symb)
      nhelp2 = nvb(syma)
      nhelp3 = nvb(symc)
      nhelp4 = dimm(v%d(0,1),syma)
      nhelp5 = dimm(v%d(0,2),symb)
      nhelp6 = dimm(v%d(0,3),symc)
      nhelp7 = nvb(syma)-nhelp4
      nhelp8 = nvb(symc)-nhelp6

      !4.1.* realize definition of V
      if (syma > symc) then
        call defvhlp81(wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
      else
        call defvhlp82(wrk(posr2),wrk(posv),nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7,nhelp8)
      end if

    end if

  end do

else
  ! not implemented deftyp

end if

return

end subroutine defv
