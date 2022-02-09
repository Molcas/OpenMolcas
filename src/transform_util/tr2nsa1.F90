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

subroutine tr2NsA1(CMO,X1,nX1,X2,nX2,X3,nX3,pqUS,npqUS,pqRU,npqRU,pqTU,npqTU,lBuf)
! SECOND ORDER TWO-ELECTRON TRANSFORMATION ROUTINE
!
! THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
! (ISP,ISQ,ISR,ISS) WITH ISP >= ISQ AND ISR >= ISS.
! P,Q,R,S are SO indices.
! A,B are MO indices, counting only non-frozen and non-deleted.
! T,U are occupied MO indices, only non-frozen and non-deleted.
! INTEGRALS (AB/TU) ARE ALWAYS GENERATED
! EXCHANGE INTEGRALS (AT/BU) ARE GENERATED AS FOLLOWS:
! (AT/BU) IF ISP >= ISR
! (AT/UB) IF ISP > ISS AND ISP /= ISQ
! (TA/BU) IF ISQ > ISR AND ISP /= ISQ
! (TA/UB) IF ISQ >= ISS AND ISP /= ISQ

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "rasdim.fh"
#include "caspt2.fh"
integer(kind=iwp), intent(in) :: nX1, nX2, nX3, npqUS, npqRU, npqTU, lBuf
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(out) :: X1(nX1), X2(nX2), X3(nX3)
real(kind=wp), intent(inout) :: pqUS(npqUS), pqRU(npqRU), pqTU(npqTU)
integer(kind=iwp) :: IAD1, IAD1S, IAD2, IAD2S, IAD3, IAD3S, icc, icxc1, icxc3, icxc5, icxc7, iOpt, IOUT1, IOUT2, IOUT3, IP, &
                     IPQMX1, IPQMX2, IPQMX3, IQ, iRc, IRSST, LPQ, NORU, NOTU, NOUS, Num
#include "trafo.fh"

NOTU = NOCR*NOCS
if (ISR == ISS) NOTU = (NOCR**2+NOCR)/2
NOUS = NOCR*NBS
NORU = NBR*NOCS
icc = NOP*NOQ*NOCR*NOCS
icxc1 = NOP*NOCQ*NOR*NOCS
icxc3 = NOP*NOCQ*NOCR*NOS
icxc5 = NOCP*NOQ*NOR*NOCS
icxc7 = NOCP*NOQ*NOCR*NOS

! Check for in core or out of core transformation

! 1. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|US) ON UNIT LUHLF1
IPQMX1 = NBPQ
!vv prevent integer overflow
if (One*NBPQ*NOUS > LURPQ) then
  IPQMX1 = LURPQ/NOUS
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|US)',IPQMX1
  IAD1S = 0
  call dDAFILE(LUHLF1,0,pqUS,IPQMX1,IAD1S)
end if
IAD1 = 0
IOUT1 = 0
! 2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
IPQMX2 = NBPQ
!vv prevent integer overflow
if (One*NBPQ*NORU > LRUPQ) then
  IPQMX2 = LRUPQ/NORU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
  IAD2S = 0
  call dDAFILE(LUHLF2,0,pqRU,IPQMX2,IAD2S)
end if
IAD2 = 0
IOUT2 = 0
! 3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
IPQMX3 = NBPQ
!vv prevent integer overflow
if (One*NBPQ*NOTU > LTUPQ) then
  IPQMX3 = LTUPQ/NOTU
  !write(u6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
  IAD3S = 0
  call dDAFILE(LUHLF3,0,PQTU,IPQMX3,IAD3S)
end if
IAD3 = 0
IOUT3 = 0

! First half transformation

LPQ = 0
NPQ = 0
iRc = 0
iOpt = 1
IRSST = 1-NBRS
! Loop over p,q symmetry pair, if(ISP == ISQ) loop should be triangle
do IP=1,NBP
  Num = NBQ
  if (ISP == ISQ) Num = IP
  do IQ=1,Num
    IOUT1 = IOUT1+1
    IOUT2 = IOUT2+1
    IOUT3 = IOUT3+1
    ! Read integrals (pq,rs)
    if (LPQ == NPQ) then
      call Rdord(iRc,iOpt,ISP,ISQ,ISR,ISS,X1,lBuf,nPQ)
      if (IRC > 1) then
        write(u6,*) ' ERROR RETURN CODE IRC=',IRC
        write(u6,*) ' FROM RDORD, CALLED FROM TRA2.'
        call Abend()
      end if
      iOpt = 2
      LPQ = 0
      IRSST = 1-NBRS
    end if
    LPQ = LPQ+1
    IRSST = IRSST+NBRS
    ! Square if necessary
    if (ISR == ISS) then
      call Square(X1(IRSST),X2,1,NBS,NBS)
    else
      call dcopy_(NBRS,X1(IRSST),1,X2,1)
    end if

    !====================================================
    ! First half transformation to (pq,Us), if ISR /= ISS
    ! For exchange case3,4 (AT,UB), case7,8 (TA,UB)
    !====================================================
    if ((icxc3 /= 0) .or. (icxc7 /= 0)) then
      if (ISR /= ISS) then
        !  (pq,rs) -> (pq,Us)
        call DGEMM_('N','N',NBS,NOCR,NBR,One,X2,NBS,CMO(LMOR2),NBR,Zero,X3,NBS)
        !  (pq,Us) Sorting
        if (IOUT1 > IPQMX1) then
          IOUT1 = 1
          call dDAFILE(LUHLF1,1,pqUS,IPQMX1*NOUS,IAD1)
        end if
        call DCOPY_(NOUS,X3,1,pqUS(IOUT1),IPQMX1)
      end if
    end if

    ! First half transformation to (pq,rU)
    ! For coulomb (AB,TU), exchange case1,2 (AT,BU), case5,6 (TA,BU)

    if ((icc /= 0) .or. (icxc1 /= 0) .or. (icxc5 /= 0)) then
      ! (pq,rs) -> (pq,rU)
      call DGEMM_('T','N',NBR,NOCS,NBS,One,X2,NBS,CMO(LMOS2),NBS,Zero,X3,NBR)
      ! (pq,rU) Sorting
      if (IOUT2 > IPQMX2) then
        IOUT2 = 1
        call dDAFILE(LUHLF2,1,pqRU,IPQMX2*NORU,IAD2)
      end if
      call DCOPY_(NORU,X3,1,pqRU(IOUT2),IPQMX2)

      ! First half transformation to (pq,TU), if(ISR == ISS) then triangle
      ! (pq,rU) -> (pq,TU)
      if (icc /= 0) then
        if (ISR == ISS) then
          call DGEMM_Tri('T','N',NOCR,NOCR,NBR,One,X3,NBR,CMO(LMOR2),NBR,Zero,X2,NOCR)
        else
          call DGEMM_('T','N',NOCS,NOCR,NBR,One,X3,NBR,CMO(LMOR2),NBR,Zero,X2,NOCS)
        end if
        ! (pq,TU) Sorting
        if (IOUT3 > IPQMX3) then
          IOUT3 = 1
          call dDAFILE(LUHLF3,1,PQTU,IPQMX3*NOTU,IAD3)
        end if
        call DCOPY_(NOTU,X2,1,PQTU(IOUT3),IPQMX3)
      end if
    end if
  end do
end do
!  Store last buffer
if (IPQMX1 < NBPQ) then
  call dDAFILE(LUHLF1,1,pqUS,IPQMX1*NOUS,IAD1)
end if
if (IPQMX2 < NBPQ) then
  call dDAFILE(LUHLF2,1,pqRU,IPQMX2*NORU,IAD2)
end if
if (IPQMX3 < NBPQ) then
  call dDAFILE(LUHLF3,1,PQTU,IPQMX3*NOTU,IAD3)
end if

return

end subroutine tr2NsA1
