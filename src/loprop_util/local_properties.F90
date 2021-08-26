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

subroutine Local_Properties(Coor,nAtoms,ip_sq_Mu,nElem,Sq_Temp,Origin,iCenter,Ttot_Inv,Temp,nij,nPert,ip_D,rMP,lMax,rMPq,C_o_C,EC, &
                            iANr,Standard,nBas1,nTemp,Q_Nuc,Bond_Threshold,Utility,Opt_Method,iPlot,iPrint,nSym)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nElem, ip_sq_Mu(0:nElem-1), nBas1, iCenter(nBas1), nij, nPert, ip_D(0:6), lMax, &
                                 iANr(nAtoms), nTemp, iPlot, iPrint, nSym
real(kind=wp), intent(in) ::  Coor(3,nAtoms), Origin(3,0:lMax), Ttot_Inv(nBas1**2), C_o_C(3), Q_Nuc(nAtoms), Bond_Threshold
real(kind=wp), intent(out) ::  Sq_Temp(nTemp), Temp(nTemp)
real(kind=wp), intent(inout) :: rMP(nij,0:nElem-1,0:nPert-1), rMPq(0:nElem-1), EC(3,nij)
logical(kind=iwp), intent(in) :: Standard, Utility
character(len=12), intent(in) :: Opt_Method
#include "Molcas.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAtom, ii, ii_, ij, ij_, iMu, iOffD, iOffO, ip_center, ip_Charge, iPert, iPL, ipNBFpA, iScratch_1, &
                     iScratch_2, iScratch_A, iScratch_B, iT_sets, iT_values, iWarnings, ix, ixnrMP, ixrMP, ixxrMP, iy, j, jAtom, &
                     ji_, jj, l, mElem, NBAST, nNUC, Num_Warnings, tNuc
real(kind=wp) :: A(3), Acc, B(3)
character(len=LenIn) :: CNAME(MxAtom)
real(kind=wp), parameter :: Ref(3) = [Zero, Zero, Zero]
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Binom()
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over perturbations

do iPert=0,nPert-1

  ! Transform the density matrix to the LoProp basis

  iOffD = ip_D(iPert)-1
  do i=1,nBas1
    do j=1,i-1
      ij = i*(i-1)/2+j
      ij_ = (j-1)*nBas1+i
      ji_ = (i-1)*nBas1+j
      Sq_Temp(ij_) = Half*Work(ij+iOffD)
      Sq_Temp(ji_) = Sq_Temp(ij_)
    end do
    ii = i*(i+1)/2
    ii_ = (i-1)*nBas1+i
    Sq_Temp(ii_) = Work(ii+iOffD)
  end do

  call DGEMM_('N','T',nBas1,nBas1,nBas1,One,Sq_Temp,nBas1,Ttot_Inv,nBas1,Zero,Temp,nBas1)
  call DGEMM_('N','N',nBas1,nBas1,nBas1,One,Ttot_Inv,nBas1,Temp,nBas1,Zero,Sq_Temp,nBas1)

  ! vv
  call Get_iScalar('Unique atoms',nNUC)
  call Get_cArray('Unique Atom Names',CNAME,(LenIn)*nNuc)

  call Get_iScalar('LP_nCenter',tNuc)
  ! someday this code will use symmetry

  iPL = iPrintLevel(-1)
  if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

  if ((nSym == 1) .and. (tNuc == nNuc) .and. (iPL >= 2)) then

    NBAST = nBas1
    call Allocate_iWork(ip_center,NBAST)
    call Get_iArray('Center Index',iWork(ip_center),NBAST)
    call Allocate_iWork(ipNBFpA,tNUC)
    do I=1,tNUC
      iWork(ipNBFpA+I-1) = 0
    end do
    do I=1,NBAST
      iWork(ipNBFpA+iWork(ip_center+I-1)-1) = iWork(ipNBFpA+iWork(ip_center+I-1)-1)+1
    end do
    call Allocate_Work(ip_Charge,tNuc)
    call Get_dArray('Effective nuclear charge',Work(ip_Charge),tNuc)
    call Free_Work(ip_Charge)
    call Free_Work(ipNBFpA)
    call Free_Work(ip_center)
  end if
  ! vv
  iMu = -1
  do l=0,lMax
    do ix=l,0,-1
      do iy=l-ix,0,-1
        iMu = iMu+1

        ! Compute local properties

        do iAtom=1,nAtoms
          call dcopy_(3,Coor(1,iAtom),1,A,1)
          do jAtom=1,iAtom
            call dcopy_(3,Coor(1,jAtom),1,B,1)

            ! Sum up contributions to the domain ij

            Acc = Zero
            iOffO = ip_sq_mu(iMu)-1
            do j=1,nBas1
              do i=1,nBas1
                if (((iCenter(i) == iAtom) .and. (iCenter(j) == jAtom)) .or. &
                    ((iCenter(i) == jAtom) .and. (iCenter(j) == iAtom))) then
                  ij = (j-1)*nBas1+i
                  Acc = Acc+Sq_Temp(ij)*Work(ij+iOffO)
                end if
              end do
            end do
            ij = iAtom*(iAtom-1)/2+jAtom
            rMP(ij,iMu,iPert) = -Acc

          end do   ! jAtom
        end do   ! iAtom

      end do   ! iy
    end do   ! ix
  end do   ! l
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Free_Work(ip_D(iPert))
end do   ! iPert
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify all multipole moments
!
! 1) Modify first all to a common origin (0.0,0.0,0.0)
! 2) Modify to local moments

#ifdef _DEBUGPRINT_
call xSpot('Middle  Local_Properties')
call RecPrt('rMP',' ',rMP,nij,nElem*nPert)
#endif
do iAtom=1,nAtoms
  ii = iAtom*(iAtom+1)/2

  ! Observe the loop order to make sure that element (ii) and
  ! (jj) always are processed before (ij)!

  do jAtom=iAtom,1,-1
    ij = iAtom*(iAtom-1)/2+jAtom
    jj = jAtom*(jAtom+1)/2
    !write(u6,*) 'ij=',ij

    ! Step 1

    do l=1,lMax-1
      mElem = (l+1)*(l+2)*(l+3)/6
      do iPert=0,nPert-1
        call ReExpand(rMP(1,0,iPert),nij,mElem,Origin(1,l),Origin(1,l+1),ij,l)
      end do
      if (ij == 1) call ReExpand(rMPq,1,mElem,Origin(1,l),Origin(1,l+1),1,l)
    end do
    !write(u6,*) 'End step 1a'
    mElem = (lMax+1)*(lMax+2)*(lMax+3)/6
    do iPert=0,nPert-1
      call ReExpand(rMP(1,0,iPert),nij,mElem,Origin(1,lMax),Ref,ij,lMax)
    end do
    if (ij == 1) call ReExpand(rMPq,1,mElem,Origin(1,lMax),C_o_C,1,lmax)
    !write(u6,*) 'End step 1b'

    ! Step 2
    !
    ! Establish the initial expansion centers for each domain.
    !
    ! Default values

    if (iAtom == jAtom) then
      A(1) = Coor(1,iAtom)
      A(2) = Coor(2,iAtom)
      A(3) = Coor(3,iAtom)
    else
      A(1) = (EC(1,ii)+EC(1,jj))*Half
      A(2) = (EC(2,ii)+EC(2,jj))*Half
      A(3) = (EC(3,ii)+EC(3,jj))*Half
    end if
    call dcopy_(3,A,1,EC(1,ij),1)
    do iPert=0,nPert-1
      call ReExpand(rMP(1,0,iPert),nij,mElem,Ref,A,ij,lMax)
    end do
    !write(u6,*) 'End step 2'

  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Distributes the contributions from the bonds that doesn't fulfill the requirement
! Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
! two atoms involved in the bond.

call Move_Prop(rMP,EC,lMax,nElem,nAtoms,nPert,nij,iANr,Bond_Threshold)

! Modify the expansion centers

Num_Warnings = 0
call Allocate_Work(iT_values,nij)
call Allocate_iWork(iT_sets,nij)
call Allocate_iWork(iWarnings,nij)
call iCopy(nij,[0],0,iWork(iWarnings),1)
call iCopy(nij,[0],0,iWork(iT_Sets),1)
call dCopy_(nij,[Zero],0,Work(iT_Values),1)
if (.not. Standard) then
  call Allocate_Work(iScratch_1,nij*(2*lMax+1))
  call Allocate_Work(iScratch_2,nij*(2*lMax+1))
  call Allocate_Work(iScratch_A,3*nij)
  call Allocate_Work(iScratch_B,3*nij)
  call Allocate_Work(ixrMP,nij*nElem)
  call Allocate_Work(ixxrMP,nij*nElem)
  call Allocate_Work(ixnrMP,nij*nElem)

  call Move_EC(rMP,EC,Work(iScratch_1),Work(iScratch_2),Work(ixrMP),Work(ixxrMP),Work(ixnrMP),lMax,Work(iScratch_A), &
               Work(iScratch_B),nij,nElem,Coor,nAtoms,Q_Nuc,C_o_C,nPert,Bond_Threshold,iANr,Work(iT_Values),iWork(iT_Sets), &
               iWork(iWarnings),Num_Warnings,Opt_Method,iPlot,iPrint)

  call Free_Work(ixnrMP)
  call Free_Work(ixxrMP)
  call Free_Work(ixrMP)
  call Free_Work(iScratch_B)
  call Free_Work(iScratch_A)
  call Free_Work(iScratch_2)
  call Free_Work(iScratch_1)
end if
if ((iPrint >= 1) .or. (iPlot >= 1) .or. (Num_Warnings > 0)) then
  call Print_T_Values(Work(iT_Values),iWork(iT_Sets),iANr,EC,Bond_Threshold,nAtoms,nij,Standard,iWork(iWarnings),Num_Warnings, &
                      iPrint)
end if
call Free_iWork(iWarnings)
call Free_iWork(iT_values)
call Free_iWork(iT_sets)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('EC',' ',EC,3,nij)
call RecPrt('rMP',' ',rMP,nij,nElem*nPert)
call RecPrt('rMPq',' ',rMPq,1,nElem)
call xSpot('Exit  Local_Properties')
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_logical(Utility)

end subroutine Local_Properties
