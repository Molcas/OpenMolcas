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
! Copyright (C) 1996,1997, Anders Bernhardsson                         *
!***********************************************************************

subroutine Prec_dig(rpre,ipre,idsym)
!***********************************************************************
!                                                                      *
!     idsym, symmetry of orbital hessian of interest                   *
!     CMtx preconditioner                                              *
!                                                                      *
!     The orbital hessian is dominated of elements that couples        *
!                                                                      *
!     kappa  -> kappa      where i is occupied and p,q is general.     *
!          ip        iq                                                *
!                                                                      *
!     we therefore approximate the hessian with those diagonal         *
!     terms in the preconditioner                                      *
!                                                                      *
!     Anders Bernhardsson 96                                           *
!                                                                      *
!     active; active,general is needed for rasscf calculation          *
!     and is not coded yet (ugly bastard) (970109, AB)                 *
!***********************************************************************

use Symmetry_Info, only: Mul
use MCLR_Data, only: F0SQMO, FAMO, FIMO, ipCM, nrec
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nRS1, nRS2, nRS3, nSym, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: rpre(*)
integer(kind=iwp), intent(_OUT_) :: ipre(*)
integer(kind=iwp), intent(in) :: idsym
integer(kind=iwp) :: iB, iBB, ip, ipi, iRC, iS, jS, n2, nD, ni, nmm, nmmm, nTemp
real(kind=wp) :: Sgn
real(kind=wp), allocatable :: JInt(:), KInt(:), Scr(:), Temp1(:,:), Temp2(:), Temp3(:), Temp4(:)

!                                                                      *
!***********************************************************************
!                                                                      *
nmm = max(0,maxval(nAsh(1:nSym)+nIsh(1:nSym)))
nmmm = max(0,maxval(nBas(1:nSym)))
n2 = nMMM**2
nmmm = ((nmmm-1)/nRec+1)*nRec
nmm = nmm*nMMM
nmm = nmm**2

call mma_allocate(JInt,n2,Label='JInt')
call mma_allocate(KInt,n2,Label='KInt')
call mma_allocate(Scr,n2,Label='Scr')

ip = 1
ipi = 1
Sgn = One
do iS=1,nSym
  jS = Mul(is,iDSym)
  nD = nBas(js)-nIsh(jS)
  ni = nBas(js)**2
  call mma_allocate(Temp2,ni,Label='Temp2')
  call mma_allocate(Temp3,ni,Label='Temp3')
  call mma_allocate(Temp4,ni,Label='Temp4')
  Temp4(:) = Zero
  call mma_MaxDBLE(nTemp)
  nTemp = min(nmm,nTemp/2)
  call mma_allocate(Temp1,nTemp,2,Label='Temp1')

  if (nD /= 0) then
    do iB=1,nIsh(iS)
      Temp3(1:nD**2) = Zero
      ibb = nBas(is)*(ib-1)+ib-2

      if (iMethod == 2) then

        ! G
        !  iaib

        if (nash(js) > 0) &
          call Preciaa(ib,is,js,nd,Temp3,nbas(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),FIMO(ipCM(js)),FAMO(ipCM(js)), &
                       F0sqMO(ipCM(js)),Sgn,JInt,KInt,Scr,n2) ! OK

        ! G
        !  ipia

        if ((nbas(js)-nish(js)-nash(js))*nash(js) > 0) &
          call Preciba(ib,is,js,nd,Temp3,nbas(js),FIMO(ipCM(js)),FAMO(ipCM(js)),F0sqMO(ipCM(js)),Sgn,JInt,KInt,Scr,n2) ! OK
      end if

      ! G
      !  ipiq

      if ((nbas(js)-nish(js)-nash(js)) > 0) &
        call Precibb_td(ib,is,js,nd,Temp3,nBas(js),Temp1(:,1),Temp1(:,2),Temp2,FiMo(1+ipCM(is)+ibb),FAMO(1+ipcm(is)+ibb), &
                        FiMo(ipCM(js)),FAMO(ipcm(js)),Sgn)  ! OK

      ! Factorize G:
      !
      !     T
      ! G=LL

      if (.not. timedep) then
        call SQM(Temp3,rpre(ip),nd)
        irc = 0
        call dgetrf_(nd,nd,rpre(ip),nd,ipre(ipi),irc)
        if (irc /= 0) then
          write(u6,*) 'Error in DGETRF called from prec_dig'
          call Abend()
        end if
      else
        call SQM(Temp3,Temp4,nD)
        rpre(ip:ip+nd-1) = Temp4(1:nd**2:nd+1)
      end if
      if (TimeDep) then
        ip = ip+nD
      else
        ip = ip+nD**2
        ipi = ipi+nD
      end if

    end do   ! iB, inactive
  end if

  Temp4(1:ni) = Zero
  do iB=1,nAsh(iS)
    ibb = nBas(is)*(nish(is)+ib-1)+nish(is)+ib-2
    if (ib <= nRS1(iS)) then
      nD = nBas(js)-nRs1(js)
    else if (ib <= nRs1(iS)+nRs2(iS)) then
      nD = nBas(js)-nRs2(js)
    else
      nD = nBas(js)-nRs3(js)
    end if
    if (nD /= 0) then
      Temp3(1:nD**2) = Zero
      if (nish(js) > 0) &
        call Precaii(ib,is,js,nd,Temp3,nbas(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),F0SqMO(1+ipCM(is)+ibb),FIMO(ipCM(js)), &
                     FAMO(ipCM(js)),Sgn,JInt,KInt,Scr,n2) ! OK
      !call Precaai(ib,nd,ir,rpre(ip))
      !call Precaaa(ib,nd,ir,rpre(ip))
      if (nish(js)*nBas(js) > 0) &
        call Precabi(ib,is,js,nd,Temp3,nBas(js),FIMO(ipCM(js)),FAMO(ipCM(js)),Sgn,JInt,KInt,Scr,n2) !+/-?

      !call Precaba(ib,nd,ir,rpre(ip))
      if (nBas(js) > 0) &
        call Precabb(ib,is,js,nd,nbas(js),Temp3,Temp1(:,1),ntemp,Temp1(:,2),Temp2,F0SQMO(1+ipCM(is)+ibb),FiMo(ipCM(js)),Sgn)
      if (.not. timedep) then
        call SQM(Temp3,rpre(ip),nD)
        irc = 0
        call dgetrf_(nd,nd,rpre(ip),nd,ipre(ipi),irc)
        if (irc /= 0) then
          write(u6,*) 'Error in DGETRF called from prec_dig'
          call Abend()
        end if
      else
        ! From Triang mat
        call SQM(Temp3,Temp4,nD)
        rpre(ip:ip+nd-1) = Temp4(1:nd**2:nd+1)
      end if
      if (timedep) then
        ip = ip+nd
      else
        ip = ip+nD**2
        ipi = ipi+nD
      end if
    end if
  end do ! iB
  call mma_deallocate(Temp1)
  call mma_deallocate(Temp2)
  call mma_deallocate(Temp3)
  call mma_deallocate(Temp4)
end do ! End loop over symmetries
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Scr)
call mma_deallocate(KInt)
call mma_deallocate(JInt)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Prec_dig
