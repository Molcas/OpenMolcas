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

subroutine Prec_dig(rpre,idsym)
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

use MCLR_Data, only: nrec
use input_mclr, only: nSym, nAsh, nIsh, nBas, nRS1, nRS2, nRS3, iMethod, TimeDep

implicit none
real*8 rpre(*)
integer idsym
real*8, allocatable :: JInt(:), KInt(:), Scr(:)
real*8, allocatable :: Temp1(:,:), Temp2(:), Temp3(:), Temp4(:)

!                                                                      *
!***********************************************************************
!                                                                      *

call Prec_dig_internal(rpre)

! This is to allow type punning without an explicit interface
contains

subroutine Prec_dig_internal(rpre)

  use iso_c_binding
  use Arrays, only: FAMO, FIMO, F0SQMO
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use MCLR_Data, only: ipCM

  implicit none
  real*8, target :: rpre(*)
  integer, pointer :: ipre(:)
  integer nmm, nmmm, iS, n2, ip, jS, nD, ni, nTemp, iB, iBB, iRC, iR
  real*8 Sign

  nmm = 0
  nmmm = 0
  do iS=1,nSym
    nMM = max(nMM,nAsh(is)+nIsh(iS))
    nMMM = max(nmmM,nBas(is))
  end do
  n2 = nMMM**2
  nmmm = ((nmmm-1)/nRec+1)*nRec
  nmm = nmm*nMMM
  nmm = nmm**2

  call mma_allocate(JInt,n2,Label='JInt')
  call mma_allocate(KInt,n2,Label='KInt')
  call mma_allocate(Scr,n2,Label='Scr')

  ip = 1
  sign = 1.0d0
  do iS=1,nSym
    jS = ieor(is-1,iDSym-1)+1
    nD = nBas(js)-nIsh(jS)
    ni = nBas(js)**2
    call mma_allocate(Temp2,ni,Label='Temp2')
    call mma_allocate(Temp3,ni,Label='Temp3')
    call mma_allocate(Temp4,ni,Label='Temp4')
    Temp4(:) = 0.0d0
    call mma_MaxDBLE(nTemp)
    nTemp = min(nmm,nTemp/2)
    call mma_allocate(Temp1,nTemp,2,Label='Temp1')

    if (nD /= 0) then
      do iB=1,nIsh(iS)
        Temp3(1:nD**2) = 0.0d0
        ibb = nBas(is)*(ib-1)+ib-2

        if (iMethod == 2) then

          ! G
          !  iaib

          if (nash(js) > 0) &
            call Preciaa(ib,is,js,nd,Temp3,nbas(is),nbas(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),F0sqMO(1+ipCM(is)+ibb), &
                         FIMO(ipCM(js)),FAMO(ipCM(js)),F0sqMO(ipCM(js)),sign,JInt,KInt,Scr,n2) ! OK

          ! G
          !  ipia

          if ((nbas(js)-nish(js)-nash(js))*nash(js) > 0) &
            call Preciba(ib,is,js,nd,Temp3,nbas(js),FIMO(ipCM(js)),FAMO(ipCM(js)),F0sqMO(ipCM(js)),sign,JInt,KInt,Scr,n2) ! OK
        end if

        ! G
        !  ipiq

        if ((nbas(js)-nish(js)-nash(js)) > 0) &
          call Precibb_td(ib,is,js,nd,Temp3,nBas(js),Temp1(:,1),Temp1(:,2),Temp2,FiMo(1+ipCM(is)+ibb),FAMO(1+ipcm(is)+ibb), &
                          FiMo(ipCM(js)),FAMO(ipcm(js)),sign)  ! OK

        ! Factorize G:
        !
        !     T
        ! G=LL

        if (.not. timedep) then
          call SQM(Temp3,rpre(ip),nd)
#         ifdef RS6K
          call DGEF(rPre(ip),nD,nD,rpre(ip+nD**2))
#         else
          irc = 0
          call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
          call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
          nullify(ipre)
          if (irc /= 0) then
            write(6,*) 'Error in DGETRF called from prec_dig'
            call Abend()
          end if
#         endif
        else
          call SQM(Temp3,Temp4,nD)
          call SortOutDiagonal(Temp4,rpre(ip),nd)
        end if
        if (TimeDep) then
          ip = ip+nD
        else
          ip = ip+nD*(nd+1)
        end if

      end do   ! iB, inactive
    end if

    Temp4(1:ni) = 0.0d0
    do iB=1,nAsh(iS)
      ibb = nBas(is)*(nish(is)+ib-1)+nish(is)+ib-2
      if (ib <= nRs1(iS)+nRs2(is)+nRs3(is)) iR = 3
      if (ib <= nRs1(iS)+nRs2(is)) iR = 2
      if (ib <= nRs1(iS)) iR = 1
      if (ir == 1) nD = nBas(js)-nRs1(js)
      if (ir == 2) nD = nBas(js)-nRs2(js)
      if (ir == 3) nD = nBas(js)-nRs3(js)
      if (nD /= 0) then
        Temp3(1:nD**2) = 0.0d0
        if (nish(js) > 0) &
          call Precaii(ib,is,js,nd,ir,Temp3,nbas(is),nbas(js),FIMO(1+ipCM(is)+ibb),FAMO(1+ipCM(is)+ibb),F0SqMO(1+ipCM(is)+ibb), &
                       FIMO(ipCM(js)),FAMO(ipCM(js)),F0SqMO(ipCM(js)),sign,JInt,KInt,Scr,n2) ! OK
        !call Precaai(ib,nd,ir,rpre(ip))
        !call Precaaa(ib,nd,ir,rpre(ip))
        if (nish(js)*nBas(js) > 0) &
          call Precabi(ib,is,js,ir,nd,Temp3,nBas(js),FIMO(ipCM(js)),FAMO(ipCM(js)),F0SQMO(ipCM(js)),sign,JInt,KInt,Scr,n2) !+/-?

        !call Precaba(ib,nd,ir,rpre(ip))
        if (nBas(js) > 0) &
          call Precabb(ib,is,js,nd,nbas(js),Temp3,Temp1(:,1),ntemp,Temp1(:,2),Temp2,F0SQMO(1+ipCM(is)+ibb),FiMo(ipCM(js)), &
                       FAMO(ipcm(js)),F0SQMO(ipCM(js)),sign)
        if (.not. timedep) then
          call SQM(Temp3,rpre(ip),nD)
#         ifdef RS6K
          call DGEF(rPre(ip),nD,nd,rpre(ip+nd**2))
#         else
          irc = 0
          call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
          call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
          nullify(ipre)
          if (irc /= 0) then
            write(6,*) 'Error in DGETRF called from prec_dig'
            call Abend()
          end if
#         endif
        else
          ! From Triang mat
          call SQM(Temp3,Temp4,nD)
          call SortOutDiagonal(Temp4,rpre(ip),nd)
        end if
        if (timedep) then
          ip = ip+nd
        else
          ip = ip+nD*(nd+1)
        end if
      end if
    end do ! iB
    call mma_deallocate(Temp1)
    call mma_deallocate(Temp2)
    call mma_deallocate(Temp3)
    call mma_deallocate(Temp4)
  end do ! End loop over symmetries
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_deallocate(Scr)
  call mma_deallocate(KInt)
  call mma_deallocate(JInt)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine Prec_dig_internal

subroutine SortOutDiagonal(Matrix,diagonal,nb)

  ! Copy the diagonal elements from Matrix to the vector Diagonal

  implicit none
  integer nb
  real*8 Matrix(*), diagonal(nb)
  integer ipM, i

  ipM = 1
  do i=1,nb
    Diagonal(i) = Matrix(ipM)
    ipM = ipM+nb+1
  end do

end subroutine SortOutDiagonal

end subroutine Prec_dig
