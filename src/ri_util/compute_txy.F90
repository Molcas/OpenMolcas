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

subroutine compute_txy(DM1,nDM,Txy,nTxy,nAuxVec,nIrrep,Diag,DMTmp,nAct)
!***********************************************************************
!                                                                      *
!     Compute the matrices needed for CD-CASSCF gradients              *
!                                                                      *
!     input : G2    = 2-body density matrix                           *
!             nDM    = size of the one-body DM                         *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use pso_stuff, only: G2, lsa, nnP
use Constants, only: One, Two, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDM, nTxy, nAuxVec, nIrrep, nAct(0:7)
real(kind=wp), intent(in) :: DM1(nDM,nAuxVec)
real(kind=wp), intent(out) :: Txy(nTxy,nAuxVec), Diag(nDM,nAuxVec), DMtmp(nTri_Elem(nDM))
integer(kind=iwp) :: i, icol, iend, iline, ista, isym, it, itu, ituvx, ituvx2, itv, itx, iu, iuv, iux, iv, iVec, ivx, ix, jend, &
                     jsta, jsym, kend, klsym, ksta, ksym, lend, lsta, lsym, nCumAct(0:7), nCumAct2(0:7), nkl, nvx, Txy_sta, Txy_sta2
real(kind=wp) :: Fac, Fac2, tmp

nCumAct(0) = 0
do i=1,nIrrep-1
  nCumAct(i) = nCumAct(i-1)+nAct(i-1)
end do

do iVec=1,nAuxVec
  !*********************************************************************
  !                                                                    *
  !   Remove Coulomb and exchange contribution from the 2DM            *
  !                                                                    *
  !*********************************************************************
  Txy_sta = 1
  Txy_sta2 = 1
  do klsym=0,nIrrep-1 ! compound symmetry
    nkl = 0

    do lSym=0,nIrrep-1
      lsta = nCumAct(lsym)+1
      lend = nCumAct(lsym)+nAct(lSym)

      ksym = Mul(lsym+1,klsym+1)-1
      if (ksym > lsym) cycle
      ksta = nCumAct(ksym)+1
      kend = nCumAct(ksym)+nAct(ksym)
      if (kSym == lSym) then
        nvx = nTri_Elem(nAct(lSym))
      else
        nvx = nAct(lSym)*nAct(ksym)
      end if
      nCumAct2(lSym) = nkl

      do jsym=0,lSym

        jsta = nCumAct(jsym)+1
        jend = nCumAct(jsym)+nAct(jSym)

        isym = Mul(jsym+1,klsym+1)-1
        if (isym > jsym) cycle
        ista = nCumAct(iSym)+1
        iend = nCumAct(iSym)+nAct(iSym)

        iline = nkl

        do ix=lsta,lend
          if (kSym == lSym) kend = ix
          do iv=ksta,kend
            ivx = iTri(ix,iv)
            if (jSym == lSym) jend = ix
            iline = iline+1
            icol = nCumAct2(jSym)
            do iu=jsta,jend
              iux = iTri(ix,iu)
              iuv = iTri(iu,iv)
              if (iSym == jSym) iend = iu
              do it=ista,iend
                itu = iTri(iu,it)
                if (itu > ivx) cycle
                itx = iTri(ix,it)
                itv = iTri(iv,it)
                ituvx = iTri(ivx,itu)
                icol = icol+1
                ituvx2 = nTri_Elem(iline-1)+icol

                Fac = One
                if (ix /= iv) Fac = Two*Fac
                if (it /= iu) Fac = Two*Fac
                Fac2 = One
                if (it == iu) Fac2 = Two

                DMTmp(ituvx2) = Fac*(Fac2*G2(ituvx,iVec))
                if (.not. lSA) then
                  ! For SA-CASSCF, don't remove Coulomb and exchange
                  if (iSym == jSym) DMTmp(ituvx2) = DMTmp(ituvx2)-Fac*(DM1(itu,iVec)*DM1(ivx,iVec))
                  if (iSym == kSym) DMTmp(ituvx2) = DMTmp(ituvx2)+Fac*(Quart*DM1(itv,iVec)*DM1(iux,iVec))
                  if (iSym == lSym) DMTmp(ituvx2) = DMTmp(ituvx2)+Fac*(Quart*DM1(itx,iVec)*DM1(iuv,iVec))
                end if

              end do
            end do
          end do
        end do

      end do
      nkl = nkl+nvx
    end do

    !*******************************************************************
    !                                                                  *
    ! Eigen-decompose the density                                      *
    !                                                                  *
    !*******************************************************************

    ! Diagonalize G2

    call unitmat(Txy(Txy_sta2:Txy_sta2+nkl**2-1,iVec),nkl)

    call NIdiag(DMTmp,Txy(Txy_sta2,iVec),nkl,nkl)

    ! Multiply by Sqrt[eigenvalue]

    do i=1,nkl
      Diag(i+Txy_sta-1,iVec) = DMTmp(nTri_Elem(i))
      tmp = sqrt(abs(DMTmp(nTri_Elem(i))))
      Txy(Txy_sta2+(i-1)*nkl:Txy_sta2+i*nkl-1,iVec) = tmp*Txy(Txy_sta2+(i-1)*nkl:Txy_sta2+i*nkl-1,iVec)
    end do

    ! Since there is no screening yet

    nnP(klsym) = nkl

    Txy_sta = Txy_sta+nkl
    Txy_sta2 = Txy_sta2+nkl**2
  end do
end do

end subroutine compute_txy
