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
! Copyright (C) 1987, Bjorn O. Roos                                    *
!***********************************************************************
!--------------------------------------------*
! 1987  B. O. ROOS                           *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND, SWEDEN                 *
!--------------------------------------------*

subroutine tr2NsB(CMO,X1,X2,pqrs,TUrs,lBuf,MAXRS)
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
!
! ********** IBM-3090 RELEASE 87 09 14 **********
! Replace MXMA with DGEMM P-AA Malmqvist 1992-05-06.
!
!
! This and tr2NsB routines transform non-squared AO integrals. The
! transformed MO integrals are stored as the same as Tr2Sq
! subroutine does.

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
#include "rasdim.fh"
#include "caspt2.fh"
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(_OUT_) :: X1(*), X2(*)
real(kind=wp), intent(inout) :: PQRS(*), TURS(*)
integer(kind=iwp), intent(in) :: lBuf, MAXRS
integer(kind=iwp) :: IAD3, IAD3S, icc, iOpt, IPQ, IPQMX3, IPQST, iRc, IRS, IRSST, ISPQRS, ITU, IX2, Kread, Length, LPQ, LRS, NOTU, &
                     NP, NQ, NR, Nread, Nrest, NRS, NS, NSYMP, NT, NU, Num, NumPQ, NumRS
#include "trafo.fh"
#include "intgrl.fh"

icc = NOCP*NOCQ*NOR*NOS

if (ISP > ISR) then

  NSYMP = NSYM*(NSYM+1)/2
  NOTU = NOCP*NOCQ
  if (ISP == ISQ) NOTU = (NOCP**2+NOCP)/2

  ! SORT OF PARTIALLY TRANSFORMED INTEGRALS (TU/RS) ON UNIT LUHLF3
  IPQMX3 = NBRS
  if (NBRS*NOTU > LTUPQ) then
    IPQMX3 = LTUPQ/NOTU
    !write(u6,*)'OUT OF CORE SORT FOR INTEGRALS (TU/RS)',IPQMX3,nbrs
    IAD3S = 0
    call dDAFILE(LUHLF3,0,TURS,IPQMX3,IAD3S)
  end if
  IAD3 = 0

  ! MaxRS should be given
  IRS = 0
  LRS = 0
  NRS = 0
  Kread = 0
  Nread = NBRS/MaxRS
  Nrest = mod(NBRS,MaxRS)
  if (Nrest == 0) then
    Nrest = MaxRS
  else
    Nread = Nread+1
  end if

  if (icc /= 0) then

    ! Loop over r,s pair
    do NR=1,NBR
      NumRS = NBS
      if (ISR == ISS) NumRS = NR
      do NS=1,NumRS
        IRS = IRS+1

        ! Loop over p,q pair
        if (LRS == NRS) then
          Kread = Kread+1
          IPQ = 0
          LPQ = 0
          NPQ = 0
          iRc = 0
          iOpt = 1
          IRSST = 1-NBRS
          do NP=1,NBP
            NumPQ = NBQ
            if (ISP == ISQ) NumPQ = NP
            do NQ=1,NumPQ
              IPQ = IPQ+1

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
              ! Here, we have all the (pq,rs) for this p,q pair
              ! Copy X1 and construct (pq,rs) for this p,q pair
              Length = MaxRS
              if (Kread == Nread) Length = Nrest
              call dcopy_(Length,X1(IRSST+MaxRS*(Kread-1)),1,PQRS(IPQ),NBPQ)
              ! End of loop over p,q pair
            end do
          end do
          LRS = 0
          NRS = Length
        end if
        LRS = LRS+1

        ! Transfer for this r,s pair
        if (ISP == ISQ) then
          call Square(PQRS(NBPQ*(LRS-1)+1),X2,1,NBP,NBP)
          ! (pq,rs) -> (pU,rs)
          call DGEMM_('T','N',NBP,NOCQ,NBQ,One,X2,NBQ,CMO(LMOQ2),NBQ,Zero,X1,NBP)
          ! (pU,rs) -> (TU,rs)
          call DGEMM_Tri('T','N',NOCP,NOCP,NBP,One,X1,NBP,CMO(LMOP2),NBP,Zero,X2,NOCP)
        else
          call dcopy_(NBPQ,PQRS(NBPQ*(LRS-1)+1),1,X2,1)
          ! (pq,rs) -> (pU,rs)
          call DGEMM_('T','N',NBP,NOCQ,NBQ,One,X2,NBQ,CMO(LMOQ2),NBQ,Zero,X1,NBP)
          ! (pU,rs) -> (TU,rs)
          call DGEMM_('T','N',NOCQ,NOCP,NBP,One,X1,NBP,CMO(LMOP2),NBP,Zero,X2,NOCQ)
        end if
        ! Store buffer
        if (IRS > IPQMX3) then
          IRS = 1
          call dDAFILE(LUHLF3,1,TURS,IPQMX3*NOTU,IAD3)
        end if
        ! Sorting
        call dcopy_(NOTU,X2,1,TURS(IRS),IPQMX3)
        ! End of loop over r,s pair
      end do
    end do
    ! Store the last buffer
    if (IPQMX3 < NBRS) then
      call dDAFILE(LUHLF3,1,TURS,IPQMX3*NOTU,IAD3)
    end if
  end if

  if (icc /= 0) then
    ISPQRS = ((ISP**2-ISP)/2+ISQ-1)*NSYMP+(ISR**2-ISR)/2+ISS
    IAD2M(1,ISPQRS) = IAD13
    ITU = 0
    ! Loop over t,u pair
    do NT=1,NOCP
      Num = NOCQ
      if (ISP == ISQ) Num = NT
      do NU=1,Num
        IPQST = 1+NBRS*ITU
        ITU = ITU+1
        ! Read buffer
        if (IPQMX3 < NBRS) then
          call RBuf_tra2(LUHLF3,TURS,NBRS,IPQMX3,NOTU,ITU,IPQST,IAD3S)
        end if
        if (ISR == ISS) then
          ! Square
          call Square(TURS(IPQST),X2,1,NBR,NBR)
          ! (TU,rs) -> (TU,sB)
          call DGEMM_('T','N',NBR,NOS,NBS,One,X2,NBS,CMO(LMOS2),NBS,Zero,X1,NBR)
          ! (TU,sB) -> (TU,AB)
          call DGEMM_Tri('T','N',NOR,NOR,NBR,One,X1,NBR,CMO(LMOR2),NBR,Zero,X2,NOR)
          IX2 = (NOR*NOR+NOR)/2
        else
          call dcopy_(NBRS,TURS(IPQST),1,X2,1)
          ! (TU,rs) -> (TU,As)
          call DGEMM_('T','N',NBR,NOS,NBS,One,X2,NBS,CMO(LMOS2),NBS,Zero,X1,NBR)
          ! (TU,As) -> (TU,AB)
          call DGEMM_('T','N',NOS,NOR,NBR,One,X1,NBR,CMO(LMOR2),NBR,Zero,X2,NOS)
          IX2 = NOR*NOS
        end if
        ! Store (TU,AB) of this t,u pair

        ! WRITE THESE BLOCK OF INTEGRALS ON LUINTM

        call GADSum(X2,IX2)
        call dDAFILE(LUINTM,1,X2,IX2,IAD13)
        ! End of Loop over t,u pair
      end do
    end do
  end if
end if

return

end subroutine tr2NsB
