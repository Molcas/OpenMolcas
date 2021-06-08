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

subroutine RDINT2_MP2(IPRX)
! SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. TEST SECTION
!
! THIS SUBROUTINE READS AND CHECKS THE RESULT OF THE SECOND ORDER
! TWO-ELECTRON TRANSFORMATION PROGRAM TRA2. IT CAN BE CALLED BY
! TR2CTL IMMEDIATELY AFTER THE CALL TO TRA2

use MBPT2_Global, only: LuIntM, nBas
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IPRX
integer(kind=iwp) :: IAD131, IAD132, IAD13C, IADC, IADX1, IADX2, ISPQRS, LENGTH, LREC, LRECX, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NT, &
                     NU, NUM
real(kind=wp), allocatable :: Tmp(:)
#include "corbinf.fh"
#include "trafo.fh"

! READ ADDRESS RECORD ON UNIT LUINTM

IAD13 = 0
call iDAFILE(LUINTM,2,IADOUT,3888,IAD13)

! LOOP OVER QUADRUPLES OF SYMMETRIES (NSP,NSP,NSR,NSS)

ISPQRS = 0
do NSP=1,NSYM
  NBP = NBAS(NSP)
  NOP = NORB(NSP)
  NOCP = NOCC(NSP)
  do NSQ=1,NSP
    NBQ = NBAS(NSQ)
    NOQ = NORB(NSQ)
    NOCQ = NOCC(NSQ)
    NSPQ = Mul(NSP,NSQ)
    do NSR=1,NSYM
      NBR = NBAS(NSR)
      NOR = NORB(NSR)
      NOCR = NOCC(NSR)
      NSPQR = Mul(NSPQ,NSR)
      ISR = NSR
      do NSS=1,NSR
        NBS = NBAS(NSS)
        ISPQRS = ISPQRS+1
        if (NSPQR /= NSS) cycle
        NOS = NORB(NSS)
        NOCS = NOCC(NSS)
        if (NOCP*NOCQ*NOCR*NOCS == 0) cycle

        ! FIND ADDRESSES FOR THIS SYMMETRY BLOCK

        IADC = IADOUT(3*ISPQRS-2)
        IADX1 = IADOUT(3*ISPQRS-1)
        IADX2 = IADOUT(3*ISPQRS)
        write(u6,1000) NSP,NSQ,NSR,NSS
        if (IADC == 0) then
          write(u6,1100)
        else
          write(u6,1200) IADC
          IAD13C = IADC
        end if
        if (IADX1 == 0) then
          write(u6,1110)
        else
          write(u6,1210) IADX1
          IAD131 = IADX1
        end if
        if (IADX2 == 0) then
          write(u6,1120)
        else
          write(u6,1220) IADX2
          IAD132 = IADX2
        end if
        LREC = NOR*NOS
        if (NSR == NSS) LREC = (NOR**2+NOR)/2
        LRECX = NOR*NOS
        do NT=1,NOCP
          NUM = NOCQ
          if (NSP == NSQ) NUM = NT
          do NU=1,NUM
            if (IADC /= 0) then
              call mma_allocate(Tmp,LREC,label='Tmp')
              call dDAFILE(LUINTM,2,Tmp,LREC,IAD13C)
              if (IPRX /= 0) LENGTH = LREC
              if (IPRX == 0) LENGTH = min(LREC,10)
              write(u6,1300) NT,NU,Tmp(1:LENGTH)
              call mma_deallocate(Tmp)
            end if

            ! THE LOOP ABOVE OVER T AND U RECOVERS ONE BLOCK OF INTEGRALS (AB|TU
            ! FOR EACH PAIR T,U. TO PROCESS ONE BLOCK FOR ALL A AND B THE
            ! FOLLOWING LOOP STRUCTURE IS USED:
            !   IAB = 0
            !   do NA=1,NOR
            !     NBM = NOS
            !     if (NSR == NSS) NBM = NA
            !     do NB=1,NBM
            !       IAB = IAB+1
            !     end do
            !   end do
            ! WORK(IAB) NOW CONTAINS THE INTEGRAL (AB|TU)
            if (IADX1 /= 0) then
              call mma_allocate(Tmp,LRECX,label='Tmp')
              call dDAFILE(LUINTM,2,Tmp,LRECX,IAD131)
              if (IPRX /= 0) LENGTH = LRECX
              if (IPRX == 0) LENGTH = min(LRECX,10)
              write(u6,1310) NT,NU,Tmp(1:LENGTH)
              call mma_deallocate(Tmp)
            end if
            ! THE EXCHANGE INTEGRALS OF TYPE 1 ,(AT|BU) ARE PROCESSED AS THE
            ! COULOMB INTEGRALS. IF NST.NE.NSU THERE ARE ALSO EXCHANGE
            ! INTEGRALS OF TYPE 2, (AU|BT). THE ORDERING IS STILL T,U AND A,B
            ! BUT T IS NOW THE FOURTH INDEX AND U THE SECOND
            ! EXCHANGE INTEGRALS ARE ALWAYS QUADRATIC IN A,B
            if (IADX2 /= 0) then
              call mma_allocate(Tmp,LRECX,label='Tmp')
              call dDAFILE(LUINTM,2,Tmp,LRECX,IAD132)
              if (IPRX /= 0) LENGTH = LRECX
              if (IPRX == 0) LENGTH = min(LRECX,10)
              write(u6,1320) NT,NU,Tmp(1:LENGTH)
              call mma_deallocate(Tmp)
            end if
          end do
        end do

        ! ALL INTEGRALS FOR SYMMETRY BLOCK NSP,NSQ,NSR,NSS ARE READ

      end do
    end do
  end do
end do

return

1000 format(/1x,'SYMMETRY BLOCK',4i4)
1100 format(1x,'NO COULOMB INTEGRALS FOR THIS SYMMETRY BLOCK?')
1110 format(1x,'NO EXCHAN1 INTEGRALS FOR THIS SYMMETRY BLOCK?')
1120 format(1x,'NO EXCHAN2 INTEGRALS FOR THIS SYMMETRY BLOCK?')
1200 format(1x,'ADDRESS FOR COULOMB INTEGRALS',i8)
1210 format(1x,'ADDRESS FOR EXCHAN1 INTEGRALS',i8)
1220 format(1x,'ADDRESS FOR EXCHAN2 INTEGRALS',i8)
1300 format(/1x,'COULOMB INTEGRALS FOR TU PAIR',2i3/(1x,10f10.6))
1310 format(/1x,'EXCHAN1 INTEGRALS FOR TU PAIR',2i3/(1x,10f10.6))
1320 format(/1x,'EXCHAN2 INTEGRALS FOR TU PAIR',2i3/(1x,10f10.6))

end subroutine RDINT2_MP2
