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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine COUNT_CPF(NINTGR,NSYM,NORB,MUL)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: NINTGR
integer(kind=iwp), intent(in) :: NSYM, NORB(*), MUL(8,8)
integer(kind=iwp) :: NOP, NOQ, NOR, NORBP, NOS, NSP, NSPQ, NSPQR, NSQ, NSR, NSS, NSSM, NT, NTM, NU, NUMAX, NUMIN, NV, NX, NXM

! COUNT TWO-ELECTRON INTEGRALS
NINTGR = 0
do NSP=1,NSYM
  NOP = NORB(NSP)
  do NSQ=1,NSP
    NSPQ = MUL(NSP,NSQ)
    NOQ = NORB(NSQ)
    do NSR=1,NSP
      NSPQR = MUL(NSPQ,NSR)
      NOR = NORB(NSR)
      NSSM = NSR
      if (NSR == NSP) NSSM = NSQ
      do NSS=1,NSSM
        if (NSS /= NSPQR) cycle
        NOS = NORB(NSS)
        NORBP = NOP*NOQ*NOR*NOS
        if (NORBP == 0) cycle
        do NV=1,NOR
          NXM = NOS
          if (NSR == NSS) NXM = NV
          do NX=1,NXM
            NTM = 1
            if (NSP == NSR) NTM = NV
            do NT=NTM,NOP
              NUMIN = 1
              if ((NSP == NSR) .and. (NT == NV)) NUMIN = NX
              NUMAX = NOQ
              if (NSP == NSQ) NUMAX = NT
              do NU=NUMIN,NUMAX
                NINTGR = NINTGR+1
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

return

end subroutine COUNT_CPF
