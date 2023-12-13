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

subroutine Tr_prm_cnt(idbg,nBas_Cont,nBas_Prim)

use Basis_Info, only: dbsc, nBas, nCnttp
use Symmetry_Info, only: nIrrep
use define_af, only: iTabMx
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idbg, nBas_Cont(8), nBas_Prim(0:7)
#include "Molcas.fh"
#include "rinfo.fh"
integer(kind=iwp) :: i, ia, iBas, iBasL, ic, icnt, iCnttp, iCont, idx, iOff, ip, ipbasL, iPrim, iPrint, iSym, j, k, ka, kbias, la, &
                     nBas_Cont_Tot, ncnt, nSize, nSym, numck, numcl
logical(kind=iwp) :: New_Center, New_l, New_m, Old_Center, Old_l
integer(kind=iwp), allocatable :: icaddr(:), ihelp(:,:), mcaddr(:), numb(:), numc(:)
real(kind=wp), allocatable :: Tr(:)

! contracted basis, atomic basis functions
!
! symmetry info
!
! lant(i): number of atoms in i:th symmetry bf
! expand the coefficient matrix into symmetry basis set
! auxiliary
! icaddr(i): adresses in coeff for a symmetry adapted function

!                                                                      *
! THIS ROUTINE ONLY WORKS WITH SPHERICAL FUNCTIONS. NO CARTESIAN D:S!  *
!***********************************************************************
!                                                                      *
nSym = nIrrep
iPrint = 0
if ((iprint >= 10) .or. (idbg > 0)) then
  write(idbg,*) ' in repmat',nsym
  write(idbg,*) nSym,(nBas(i),i=0,nsym-1)
  write(idbg,*) nSym,(nrBas(i),i=1,nsym)
  write(idbg,*) nSym,(nBas_Prim(i),i=0,nsym-1)
  write(idbg,*) nSym,(nBas_Cont(i),i=1,nsym)
end if

nBas_Cont_Tot = 0
do iSym=1,nSym
  nBas_Cont_Tot = nBas_Cont_Tot+nBas_Cont(iSym)
end do
ncnt = 0
do iCnttp=1,nCnttp
  do icnt=1,dbsc(iCnttp)%nCntr
    ncnt = ncnt+1
  end do
end do
call mma_allocate(icaddr,nBas_Cont_Tot,label='icaddr')
call mma_allocate(mcaddr,nBas_Cont_Tot,label='mcaddr')
call mma_allocate(numb,nBas_Cont_Tot,label='numb')
call mma_allocate(numc,nBas_Cont_Tot,label='numc')
call mma_allocate(ihelp,ncnt,iTabMx,label='ihelp')
!                                                                      *
!***********************************************************************
!                                                                      *
! set up pointer

k = 0
ia = 0  ! center index
ka = 0  ! shell index
do iCnttp=1,nCnttp
  do icnt=1,dbsc(iCnttp)%nCntr
    ia = ia+1
    do la=1,nAngr(ia)+1
      ka = ka+1
      ihelp(ia,la) = k
      k = k+nPrimr(ka)*nBasisr(ka)
    end do
  end do
end do
if ((iPrint >= 10) .or. (idbg > 0)) then
  write(idbg,*) ' Help vector'
  ia = 0
  do iCnttp=1,nCnttp
    do icnt=1,dbsc(iCnttp)%nCntr
      ia = ia+1
      write(idbg,'(10i5)') (ihelp(ia,j),j=1,nAngr(ia)+1)
    end do
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over irreps

k = 0
do iSym=1,nSym
  numck = 1
  numcl = 0
  kbias = 0

  ! Loop over basis functions in irrep

  do iCont=1,nBas_Cont(iSym)
    k = k+1
    if (iCont > 1) then
      New_Center = icent(k) /= icent(k-1)
      New_l = lnang(k) /= lnang(k-1)
      New_m = lmag(k) /= lmag(k-1)
      if (New_m) kbias = kbias-numc(k-1)
      if (New_Center .or. New_l) kbias = 0
    else
      New_Center = .true.
      New_l = .true.
    end if
    kbias = kbias+1
    ka = 0
    ia = 0
    do iCnttp=1,nCnttp
      do iCnt=1,dbsc(iCnttp)%nCntr
        ia = ia+1
        do la=1,nAngr(ia)+1
          ka = ka+1
          Old_Center = icent(k) == ia
          Old_l = lnang(k) == (la-1)
          if (idbg > 0) write(idbg,*) ' at numck',k,ia,icent(k),la-1,lnang(k),ia,New_Center,New_l
          if (Old_Center .and. Old_l) then
            numc(k) = nBasisr(ka)
            numb(k) = nPrimr(ka)
            icaddr(k) = ihelp(ia,la)+(kbias-1)*nPrimr(ka)
            if ((k > 1) .and. (kbias == 1)) numck = numcl+numck
            mcaddr(k) = numck
            numcl = nPrimr(ka)
          end if
        end do  ! la
      end do    ! iCnt
    end do      ! iCnttp
  end do        ! iCont
end do          ! iSym
k = 0
if ((iPrint >= 10) .or. (idbg > 0)) then
  ic = 0
  ip = 0
  do iSym=1,nSym
    write(idbg,*) ' symmetry',iSym
    write(idbg,*) ' numb'
    write(idbg,'(20i4)') (numb(i+ic),i=1,nBas_Cont(iSym))
    write(idbg,*) ' numc'
    write(idbg,'(20i4)') (numc(i+ic),i=1,nBas_Cont(iSym))
    write(idbg,*) ' Pointer to contraction vector'
    write(idbg,'(20i4)') (icaddr(i+ic),i=1,nBas_Cont(iSym))
    write(idbg,*) ' mcaddr'
    write(idbg,'(20i4)') (mcaddr(i+ic),i=1,nBas_Cont(iSym))
    ic = ic+nBas_Cont(iSym)
    ip = ip+nBas_Prim(iSym-1)
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nSize = 0
do iSym=1,nSym
  nSize = nSize+nBas_Cont(iSym)*nBas_Prim(iSym-1)
end do
call mma_allocate(Tr,nSize,label='Tr')
Tr(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the transformation matrix.

ibasL = 0
iOff = 0

!---- Loop over irreps

do iSym=1,nSym

  ! loop over contracted

  do iBas=1,nBas_Cont(iSym)
    ibasL = ibasL+1
    ipbasL = mcaddr(ibasL)-1

    ! loop over uncontracted

    do iPrim=1,numb(ibasL)
      ipbasL = ipbasL+1
      idx = iOff+(iBas-1)*nBas_Prim(iSym-1)+ipbasL
      !write(u6,*) iBas,ipbasL
      Tr(idx) = rCof(icaddr(ibasL)+iPrim)

    end do     ! iPrim
  end do
  iOff = iOff+nBas_Cont(iSym)*nBas_Prim(iSym-1)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call Put_dArray('NEMO TPC',Tr,nSize)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Tr)
call mma_deallocate(icaddr)
call mma_deallocate(mcaddr)
call mma_deallocate(numb)
call mma_deallocate(numc)
call mma_deallocate(ihelp)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Tr_prm_cnt
