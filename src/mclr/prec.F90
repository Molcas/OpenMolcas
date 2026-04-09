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

subroutine Prec(rpre,ipre,idsym)
!***********************************************************************
!
!  idsym, symmetry of orbital hessian of interest
!  CMtx preconditioner
!
! The orbital hessian is dominated of elements that couples
!
! kappa  -> kappa      where i is occupied and p,q is general.
!      ip        iq
!
! we therefore approximate the hessian with those diagonal
! terms in the preconditioner
!
!  Anders Bernhardsson 96
!
!     active; active,general is needed for rasscf calculation
!     and is not coded yet (ugly bastard) (970109, AB)
!***********************************************************************

use Symmetry_Info, only: Mul
use MCLR_Data, only: F0SQMO, FAMO, FIMO, ipCM, nrec
use input_mclr, only: iMethod, nAsh, nBas, NewCho, nIsh, nOrb, nOrb, nRS1, nRS2, nRS3, nSym, ntAsh
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: rpre(*)
integer(kind=iwp), intent(_OUT_) :: ipre(*)
integer(kind=iwp), intent(in) :: idsym
integer(kind=iwp) :: iAdr, iB, iBB, ip, ipi, iR, iRC, iS, jS, n2, nD, ni, nmm, nmmm, nTemp
real(kind=wp) :: Sgn
real(kind=wp), allocatable :: ActInt(:), JA(:), KA(:), Scr(:), Temp1(:,:), Temp2(:), Temp3(:)

nmm = max(0,maxval(nAsh(1:nSym)+nIsh(1:nSym)))
nmmm = max(0,maxval(nBas(1:nSym)))
n2 = nMMM**2
nmmm = ((nmmm-1)/nRec+1)*nRec
nmm = nmm*nMMM
nmm = nmm**2

call mma_allocate(JA,n2,Label='JA')
call mma_allocate(KA,n2,Label='KA')
call mma_allocate(Scr,n2,Label='Scr')
if ((nRs1(iDSym) /= 0) .or. (nRs3(iDSym) /= 0)) then
  call mma_allocate(ActInt,ntAsh**4,Label='ActInt')
  call Precaaa_Pre(ActInt,JA,Scr)
end if

ip = 1
ipi = 1
iAdr = 0
do iS=1,nSym
  jS = Mul(is,iDSym)
  nD = nOrb(js)-nIsh(jS)
  ni = nBas(js)**2
  Sgn = One
  call mma_allocate(Temp2,ni,Label='Temp2')
  call mma_allocate(Temp3,ni,Label='Temp3')
  Temp3(:) = Zero
  call mma_MaxDBLE(nTemp)
  nTemp = min(nmm,nTemp/2)
  call mma_allocate(Temp1,nTemp,2,Label='Temp1')
  if (nd /= 0) then
    do iB=1,nIsh(iS)
      Temp3(1:nD**2) = Zero
      ibb = nOrb(is)*(ib-1)+ib-2
      !                                                                *
      !*****************************************************************
      !                                                                *
      if (NewCho) then  ! Cho-Fock

        call preci_cho(jS,nD,Temp3,nOrb(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),FIMO(ipCM(js)),FAMO(ipCM(js)), &
                       F0sqMO(ipCM(js)),Sgn,JA,n2,iAdr) ! OK

      else  ! Cho-MO

        if (iMethod == 2) then
          !                                                            *
          !*************************************************************
          !                                                            *
          ! G
          !  iaib

          if (nash(js) > 0) &
            call Preciaa(ib,is,js,nd,Temp3,nOrb(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),FIMO(ipCM(js)),FAMO(ipCM(js)), &
                         F0sqMO(ipCM(js)),Sgn,JA,KA,Scr,n2) ! OK
          !                                                            *
          !*************************************************************
          !                                                            *
          ! G
          !  ipia

          if ((nOrb(js)-nish(js)-nash(js))*nash(js) > 0) &
            call Preciba(ib,is,js,nd,Temp3,nOrb(js),FIMO(ipCM(js)),FAMO(ipCM(js)),F0sqMO(ipCM(js)),Sgn,JA,KA,Scr,n2) ! OK

        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! G
        !  ipiq

        if ((nOrb(js)-nish(js)-nash(js)) > 0) &
          call Precibb(ib,is,js,nd,Temp3,norb(js),Temp1(:,1),Temp1(:,2),Temp2,FiMo(1+ipCM(is)+ibb),FAMO(1+ipcm(is)+ibb), &
                       FiMo(ipCM(js)),FAMO(ipcm(js)),Sgn)  ! OK

      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Factorize G:
      !
      !     T
      ! G=LL
      !
      !write(u6,*) 'Preconditioner i =',iB
      !do i=1,min(nd,10)
      !  write(u6,'(10F12.8)') (Temp3((j-1)*nd+i-nTri_Elem(j-1)),j=1,i)
      !end do

      call SQM(Temp3,rpre(ip),nd)

      !write(u6,*) ' ====== rpre ======'
      !do i=1,nd*nd
      !  write(u6,*) i,'rpre',rpre(ip+i-1)
      !end do

      irc = 0
      call dgetrf_(nd,nd,rpre(ip),nd,ipre(ipi),irc)
      if (irc /= 0) then
        write(u6,*) 'Error in DGETRF called from prec'
        call Abend()
      end if
      ip = ip+nD**2
      ipi = ipi+nD

    end do

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iB=1,nAsh(iS)
    ibb = nOrb(is)*(nish(is)+ib-1)+nish(is)+ib-2
    if (ib <= nRs1(iS)+nRs2(is)+nRs3(is)) iR = 3
    if (ib <= nRs1(iS)+nRs2(is)) iR = 2
    if (ib <= nRs1(iS)) iR = 1
    if (ir == 1) nD = nOrb(js)-nRs1(js)
    if (ir == 2) nD = nOrb(js)-nRs2(js)
    if (ir == 3) nD = nOrb(js)-nRs3(js)
    if (nd /= 0) then
      Temp3(1:nD**2) = Zero

      !                                                                *
      !*****************************************************************
      !                                                                *
      if (NewCho) then

        !  New Cholesky code

        call Preca_cho(ib,is,js,nd,Temp3,nOrb(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),F0SqMO(1+ipCM(is)+ibb),FIMO(ipCM(js)), &
                       FAMO(ipCM(js)),Sgn,JA,n2)

      else

        if (nish(js) > 0) &
          call Precaii(ib,is,js,nd,Temp3,nOrb(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),F0SqMO(1+ipCM(is)+ibb),FIMO(ipCM(js)), &
                       FAMO(ipCM(js)),Sgn,JA,KA,Scr,n2) ! OK
        !call Precaai(ib,nd,ir,rpre(ip))
        !call Precaaa(ib,nd,ir,rpre(ip))
        if (nish(js)*nOrb(js) > 0) &
          call Precabi(ib,is,js,nd,Temp3,nOrb(js),FIMO(ipCM(js)),FAMO(ipCM(js)),Sgn,JA,KA,Temp1(:,2),n2)!+/-?
        !call Precaba(ib,nd,ir,rpre(ip))
        if (nOrb(js) > 0) &
          call Precabb_2(ib,is,js,nd,nOrb(js),Temp3,Temp1(:,1),ntemp,Temp1(:,2),Temp2,F0SQMO(1+ipCM(is)+ibb),FiMo(ipCM(js)),Sgn)

        ! symmetry not yet
        ! Eq. (C.12e)
        if ((nRs1(iS) /= 0) .or. (nRs3(iS) /= 0)) &
          call Precaaa(ib,is,js,nd,ir,Temp3,nOrb(js),FIMO(ipCM(js)),F0SqMO(ipCM(js)),Sgn,Scr,n2,ActInt) ! OK

      end if
      !                                                                *
      !*****************************************************************
      !                                                                *

      call SQM(Temp3,rpre(ip),nD)
      irc = 0
      call dgetrf_(nd,nd,rpre(ip),nd,ipre(ipi),irc)
      if (irc /= 0) then
        write(u6,*) 'Error in DGETRF called from prec'
        call Abend()
      end if
      ip = ip+nD**2
      ipi = ipi+nD
    end if
  end do
  call mma_deallocate(Temp1)
  call mma_deallocate(Temp2)
  call mma_deallocate(Temp3)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if ((nRs1(iDSym) /= 0) .or. (nRs3(iDSym) /= 0)) call mma_deallocate(ActInt)
call mma_deallocate(Scr)
call mma_deallocate(KA)
call mma_deallocate(JA)

!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Prec
