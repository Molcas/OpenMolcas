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

subroutine Local_Properties(Coor,nAtoms,sq_Mu,nElem,Sq_Temp,Origin,iCenter,Ttot_Inv,Temp,nij,nPert,D,rMP,lMax,rMPq,C_o_C,EC,iANr, &
                            Standard,nBas1,nTemp,Q_Nuc,Bond_Threshold,Opt_Method,iPlot,iPrint,nSym)

use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nElem, nBas1, iCenter(nBas1), nij, nPert, lMax, iANr(nAtoms), nTemp, iPlot, iPrint, nSym
real(kind=wp), intent(in) :: Coor(3,nAtoms), sq_Mu(nBas1**2,0:nElem-1), Origin(3,0:lMax), Ttot_Inv(nBas1**2), C_o_C(3), &
                             Q_Nuc(nAtoms), Bond_Threshold
type(Alloc1DArray_Type), intent(inout) :: D(0:6)
real(kind=wp), intent(out) :: Sq_Temp(nTemp), Temp(nTemp)
real(kind=wp), intent(inout) :: rMP(nij,0:nElem-1,0:nPert-1), rMPq(0:nElem-1), EC(3,nij)
logical(kind=iwp), intent(in) :: Standard
character(len=12), intent(in) :: Opt_Method
#include "Molcas.fh"
integer(kind=iwp) :: i, iAtom, ii, ii_, ij, ij_, iMu, iPert, ix, iy, j, jAtom, ji_, jj, l, mElem, Num_Warnings
real(kind=wp) :: A(3), Acc
integer(kind=iwp), allocatable :: T_sets(:), Warnings(:) !, center(:), Charge(:), NBFpA(:)
real(kind=wp), allocatable :: T_values(:)
!character(len=LenIn), allocatable :: CNAME(:)
real(kind=wp), parameter :: Ref(3) = [Zero,Zero,Zero]

#include "macros.fh"
unused_var(nSym)

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

  do i=1,nBas1
    do j=1,i-1
      ij = i*(i-1)/2+j
      ij_ = (j-1)*nBas1+i
      ji_ = (i-1)*nBas1+j
      Sq_Temp(ij_) = Half*D(iPert)%A(ij)
      Sq_Temp(ji_) = Sq_Temp(ij_)
    end do
    ii = i*(i+1)/2
    ii_ = (i-1)*nBas1+i
    Sq_Temp(ii_) = D(iPert)%A(ii)
  end do

  call DGEMM_('N','T',nBas1,nBas1,nBas1,One,Sq_Temp,nBas1,Ttot_Inv,nBas1,Zero,Temp,nBas1)
  call DGEMM_('N','N',nBas1,nBas1,nBas1,One,Ttot_Inv,nBas1,Temp,nBas1,Zero,Sq_Temp,nBas1)

  ! vv
  !call Get_iScalar('Unique atoms',nNuc)
  !call mma_allocate(CNAME,nNuc,label='CNAME')
  !call Get_cArray('Unique Atom Names',CNAME,(LenIn)*nNuc)
  !
  !call Get_iScalar('LP_nCenter',tNuc)
  ! someday this code will use symmetry
  !
  !iPL = iPrintLevel(-1)
  !if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
  !
  !if ((nSym == 1) .and. (tNuc == nNuc) .and. (iPL >= 2)) then
  !  call mma_allocate(center,nBas1,label='center')
  !  call Get_iArray('Center Index',center,nBas1)
  !  call mma_allocate(NBFpA,tNuc,label='NBFpA')
  !  NBFpA(:) = 0
  !  do I=1,nBas1
  !    NBFpA(center(I)) = NBFpA(center(I))+1
  !  end do
  !  call mma_allocate(Charge,tNuc,label='Charge')
  !  call Get_dArray('Effective nuclear charge',Charge,tNuc)
  !  call mma_deallocate(Charge)
  !  call mma_deallocate(NBFpA)
  !  call mma_deallocate(center)
  !  call mma_deallocate(CNAME)
  !end if
  ! vv
  iMu = -1
  do l=0,lMax
    do ix=l,0,-1
      do iy=l-ix,0,-1
        iMu = iMu+1

        ! Compute local properties

        do iAtom=1,nAtoms
          A(:) = Coor(:,iAtom)
          do jAtom=1,iAtom

            ! Sum up contributions to the domain ij

            Acc = Zero
            do j=1,nBas1
              do i=1,nBas1
                if (((iCenter(i) == iAtom) .and. (iCenter(j) == jAtom)) .or. &
                    ((iCenter(i) == jAtom) .and. (iCenter(j) == iAtom))) then
                  ij = (j-1)*nBas1+i
                  Acc = Acc+Sq_Temp(ij)*sq_mu(ij,iMu)
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
  call mma_deallocate(D(iPert)%A)
end do   ! iPert
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify all multipole moments
!
! 1) Modify first all to a common origin (0.0,0.0,0.0)
! 2) Modify to local moments

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Middle  Local_Properties'
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
    EC(:,ij) = A(:)
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
call mma_allocate(T_values,nij,label='T_values')
call mma_allocate(T_sets,nij,label='T_sets')
call mma_allocate(Warnings,nij,label='Warnings')
T_values(:) = Zero
T_sets(:) = 0
Warnings(:) = 0
if (.not. Standard) then
  call Move_EC(rMP,EC,lMax,nij,nElem,Coor,nAtoms,Q_Nuc,C_o_C,nPert,Bond_Threshold,iANr,T_Values,T_Sets,Warnings,Num_Warnings, &
               Opt_Method,iPlot)
end if
if ((iPrint >= 1) .or. (iPlot >= 1) .or. (Num_Warnings > 0)) then
  call Print_T_Values(T_Values,T_Sets,iANr,EC,Bond_Threshold,nAtoms,nij,Standard,Warnings,Num_Warnings,iPrint)
end if
call mma_deallocate(T_values)
call mma_deallocate(T_sets)
call mma_deallocate(Warnings)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('EC',' ',EC,3,nij)
call RecPrt('rMP',' ',rMP,nij,nElem*nPert)
call RecPrt('rMPq',' ',rMPq,1,nElem)
write(u6,*)
write(u6,*) 'Exit  Local_Properties'
#endif

return

end subroutine Local_Properties
