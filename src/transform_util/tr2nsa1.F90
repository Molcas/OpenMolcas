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

implicit real*8(a-h,o-z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "trafo.fh"
#include "intgrl.fh"
#include "SysDef.fh"
dimension CMO(NCMO)
dimension X1(nX1), X2(nX2), X3(nX3)
dimension PQTU(nPQTU), pqRU(npqRU), pqUS(npqUS)

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
if (1.0d0*NBPQ*NOUS > LURPQ) then
  IPQMX1 = LURPQ/NOUS
  !write(6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|US)',IPQMX1
  IAD1S = 0
  call dDAFILE(LUHLF1,0,pqUS,IPQMX1,IAD1S)
end if
IAD1 = 0
IOUT1 = 0
! 2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
IPQMX2 = NBPQ
!vv prevent integer overflow
if (1.0d0*NBPQ*NORU > LRUPQ) then
  IPQMX2 = LRUPQ/NORU
  !write(6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
  IAD2S = 0
  call dDAFILE(LUHLF2,0,pqRU,IPQMX2,IAD2S)
end if
IAD2 = 0
IOUT2 = 0
! 3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
IPQMX3 = NBPQ
!vv prevent integer overflow
if (1.0d0*NBPQ*NOTU > LTUPQ) then
  IPQMX3 = LTUPQ/NOTU
  !write(6,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
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
        write(6,*) ' ERROR RETURN CODE IRC=',IRC
        write(6,*) ' FROM RDORD, CALLED FROM TRA2.'
        call Abend
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
        call DGEMM_('N','N',NBS,NOCR,NBR,1.0d0,X2,NBS,CMO(LMOR2),NBR,0.0d0,X3,NBS)
        !  (pq,Us) Sorting
        if (IOUT1 > IPQMX1) then
          IOUT1 = 1
          !vv dO I=1,NOUS
          !vv   call dDAFILE(LUHLF1,1,pqUS(1+IPQMX1*(I-1)),IPQMX1,IAD1)
          !vv end do
          call dDAFILE(LUHLF1,1,pqUS,IPQMX1*NOUS,IAD1)
        end if
        call DCOPY_(NOUS,X3,1,pqUS(IOUT1),IPQMX1)
      end if
    end if

    ! First half transformation to (pq,rU)
    ! For coulomb (AB,TU), exchange case1,2 (AT,BU), case5,6 (TA,BU)

    if ((icc /= 0) .or. (icxc1 /= 0) .or. (icxc5 /= 0)) then
      ! (pq,rs) -> (pq,rU)
      call DGEMM_('T','N',NBR,NOCS,NBS,1.0d0,X2,NBS,CMO(LMOS2),NBS,0.0d0,X3,NBR)
      ! (pq,rU) Sorting
      if (IOUT2 > IPQMX2) then
        IOUT2 = 1
        !vv do I=1,NORU
        !vv   call dDAFILE(LUHLF2,1,pqRU(1+IPQMX2*(I-1)),IPQMX2,IAD2)
        !vv end do
        call dDAFILE(LUHLF2,1,pqRU,IPQMX2*NORU,IAD2)
      end if
      call DCOPY_(NORU,X3,1,pqRU(IOUT2),IPQMX2)

      ! First half transformation to (pq,TU), if(ISR == ISS) then triangle
      ! (pq,rU) -> (pq,TU)
      if (icc /= 0) then
        if (ISR == ISS) then
          call MXMT(X3,NBR,1,CMO(LMOR2),1,NBR,X2,NOCR,NBR)
        else
          call DGEMM_('T','N',NOCS,NOCR,NBR,1.0d0,X3,NBR,CMO(LMOR2),NBR,0.0d0,X2,NOCS)
        end if
        ! (pq,TU) Sorting
        if (IOUT3 > IPQMX3) then
          IOUT3 = 1
          !vv do I=1,NOTU
          !vv   call dDAFILE(LUHLF3,1,PQTU(1+IPQMX3*(I-1)),IPQMX3,IAD3)
          !vv enddo
          call dDAFILE(LUHLF3,1,PQTU,IPQMX3*NOTU,IAD3)
        end if
        call DCOPY_(NOTU,X2,1,PQTU(IOUT3),IPQMX3)
      end if
    end if
  end do
end do
!  Store last buffer
if (IPQMX1 < NBPQ) then
  !vv do i=1,NOUS
  !vv   call dDAFILE(LUHLF1,1,pqUS(1+IPQMX1*(I-1)),IPQMX1,IAD1)
  !vv end do
  call dDAFILE(LUHLF1,1,pqUS,IPQMX1*NOUS,IAD1)
end if
if (IPQMX2 < NBPQ) then
  !vv do i=1,NORU
  !vv   call dDAFILE(LUHLF2,1,pqRU(1+IPQMX2*(I-1)),IPQMX2,IAD2)
  !vv end do
  call dDAFILE(LUHLF2,1,pqRU,IPQMX2*NORU,IAD2)
end if
if (IPQMX3 < NBPQ) then
  !vv do i=1,NOTU
  !vv   call dDAFILE(LUHLF3,1,PQTU(1+IPQMX3*(I-1)),IPQMX3,IAD3)
  !vv end do
  call dDAFILE(LUHLF3,1,PQTU,IPQMX3*NOTU,IAD3)
end if

return

end subroutine tr2NsA1
