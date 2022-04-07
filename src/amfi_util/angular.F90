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

subroutine angular(Lhigh,keep,makemean,bonn,breit,sameorb,ifinite,onecartx,onecarty,onecartz,powexp,coulovlp,preXZ,preY,icheckxy, &
                   icheckz,interxyz,isgnprod)
!bs COMBINES THE RADIAL INTEGRALS WITH THE ANGULAR FACTORS
!
!bs if keep=.true. then
!bs all the integrals will be kept in memory.
!bs Perhaps, there will be the option to make the
!bs transformation to the cartesian basis-sets
!bs everytime, they are required.
!bs Therefore, the integrals are kept in memory and
!bs can be further transformed, whenever required.
!bs in order not to waste to much memory, the atomic
!bs integrals are thrown away after each l,l,l,l-block

use AMFI_global, only: AOcoeffs, ipowxyz, Lblocks, Lfirst, Llast, Lmax, Lstarter, MxcontL, MxprimL, ncontrac, noccorb, occup
use index_functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lhigh, ifinite
logical(kind=iwp), intent(in) :: keep, makemean, bonn, breit, sameorb
real(kind=wp), intent(inout) :: onecartX(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax), &
                                onecartY(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax), &
                                onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax)
real(kind=wp), intent(in) :: powexp(MxprimL,MxprimL,0:Lmax,0:Lmax,0:(Lmax+Lmax+5)), coulovlp(*)
real(kind=wp), intent(out) :: preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
integer(kind=iwp), intent(out) :: icheckxy(0:Lmax,0:Lmax,0:Lmax,0:Lmax), icheckz(0:Lmax,0:Lmax,0:Lmax,0:Lmax), &
                                  interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax), isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
integer(kind=iwp) :: iangfirst, icont4, indx, indy, indz, inter2, inter3, inter4, isignM(-Lmax:Lmax), jblock, l1, l2, l3, l4, &
                     Lleftmax, Lleftmin, locstar, Lrightmax, Lrightmin, M1, m1upper, M2, m2upper, M3, M4, mblock, mblockx, &
                     mblocky, mblockz, mxangint, ncont, numbcart
real(kind=wp) :: preroots(2,0:Lmax), scratch(0:2*Lmax+1)
logical(kind=iwp) :: cleaner, NFINI
integer(kind=iwp), allocatable :: mcombcart(:,:,:,:,:), mcombina(:,:,:,:,:)
real(kind=wp), allocatable :: clebsch(:,:,:,:), ConOO(:), ConSO(:), CartOO(:), CartSO(:), AngOO(:), AngSO(:)
!bs NFINI means not finite nucleus

!bs ####################################################################
!bs   some preparation of factors needed later on..                    #
!bs ####################################################################
!bs calculate some prefactors that will be needed quite often
call mma_allocate(clebsch,[1,3],[1,2],[-Lmax,Lmax],[0,Lmax],label='clebsch')
call prefac(Lmax,preroots,clebsch)
if (ifinite /= 2) then
  !bs clean array for one electron integrals
  onecartX(:,:,:,:) = Zero
  onecartY(:,:,:,:) = Zero
  onecartZ(:,:,:,:) = Zero
  NFINI = .true.
else
  NFINI = .false.
end if

!bs generate an array with sign for (even/odd) m-values
if (mod(Lmax,2) == 0) then
  isignM(-Lmax:Lmax:2) = 1
  isignM(-Lmax+1:Lmax-1:2) = -1
else
  isignM(-Lmax:Lmax:2) = -1
  isignM(-Lmax+1:Lmax-1:2) = 1
end if
!bs ####################################################################
!bs   prefactors preXZ und preY include the factors 1/root(2)
!bs   for the +/- linear combinations of spherical harmonics
!bs ####################################################################
preXZ(:,:,:,:) = Quart
preXZ(:,:,:,0) = preXZ(:,:,:,0)*sqrt(Two)
preXZ(:,:,0,:) = preXZ(:,:,0,:)*sqrt(Two)
preXZ(:,0,:,:) = preXZ(:,0,:,:)*sqrt(Two)
preXZ(0,:,:,:) = preXZ(0,:,:,:)*sqrt(Two)
preY(:,:,:,:) = preXZ(:,:,:,:)
!bs ####################################################################
!bs   additional (-) signs from the (-i) factors  in the
!bs   (-) linear combinations   (see tosigX(Y,Z))
!bs ####################################################################
!bs   + - - -   =>   minus
preXZ(0:,:-1,:-1,:-1) = -preXZ(0:,:-1,:-1,:-1)
!bs   - + - -   =>   minus
preXZ(:-1,0:,:-1,:-1) = -preXZ(:-1,0:,:-1,:-1)
!bs   + + + -   =>   minus
preXZ(0:,0:,0:,:-1) = -preXZ(0:,0:,0:,:-1)
!bs   + + - +   =>   minus
preXZ(0:,0:,:-1,0:) = -preXZ(0:,0:,:-1,0:)
!bs   + + - -   =>   minus
preY(0:,0:,:-1,:-1) = -preY(0:,0:,:-1,:-1)
!bs   - - + +   =>   minus
preY(:-1,:-1,0:,0:) = -preY(:-1,:-1,0:,0:)
call genprexyz13(icheckxy)
call genprexyz14(icheckz,interxyz)
call genprexyz15a(icheckxy,icheckz,interxyz)
!bs ####################################################################
!bs   isgnprod gives the sign due to powers (-1)**M  this are again
!bs   angular m-values
!bs ####################################################################
do M4=-Lmax,Lmax
  if (M4 > 0) then
    inter4 = isignM(M4)
  else
    inter4 = 1
  end if
  do M3=-Lmax,Lmax
    if (M3 > 0) then
      inter3 = inter4*isignM(M3)
    else
      inter3 = inter4
    end if
    do M2=-Lmax,Lmax
      if (M2 > 0) then
        inter2 = inter3*isignM(M2)
      else
        inter2 = inter3
      end if
      do M1=-Lmax,Lmax
        if (M1 > 0) then
          isgnprod(m1,m2,m3,m4) = inter2*isignM(M1)
        else
          isgnprod(m1,m2,m3,m4) = inter2
        end if
      end do
    end do
  end do
end do
!bs ####################################################################
!bs   some preparation of factors needed later on..  finished          #
!bs ####################################################################

! set some counters
call mma_allocate(mcombina,[1,2],[-Lmax,Lmax],[-Lmax,Lmax],[-Lmax,Lmax],[-Lmax,Lmax],label='mcombina')
call mma_allocate(mcombcart,[1,2],[-Lmax,Lmax],[-Lmax,Lmax],[-Lmax,Lmax],[-Lmax,Lmax],label='mcombcart')
!bs counter for total number of cartesian integrals
numbcart = 0
!bs same orbit integrals integrals  on carteXSO carteYSO and carteSO
!bs other orbit integrals  on carteXOO carteYOO and carteOO
iangfirst = 0 ! first block of angular integrals
!bs loop over all possible < l1 l2, l3 l4 > blocks
!BS write(u6,'(A)') '   L1   L2   L3   L4'
do l1=0,Lhigh   ! improving is probably possible...
  do l2=0,Lhigh
    do l3=0,l1
      do l4=0,l2
        !bs check parity
        if (mod(l1+l2+l3+l4,2) == 0) then
          !bs check that Lleft and Lright do not always differ by more than one
          !bs a difference of two means two spin flips and is therefore not allowed
          Lleftmax = l1+l2
          Lrightmax = l3+l4
          Lleftmin = abs(l1-l2)
          Lrightmin = abs(l3-l4)
          if (((Lrightmin-Lleftmax <= 1) .and. (Lrightmax-Lleftmin > -1)) .or. &
              ((Lleftmin-Lrightmax <= 1) .and. (Lleftmax-Lrightmin > -1))) then
            !bs additional check for mean-field
            if (((l1 == l3) .and. (l2 == l4)) .or. ((l1 == l2) .and. (l3 == l4))) then
              if (l1+l3 /= 0) then
                !BS write(u6,'(4I5)') l1,l2,l3,l4
                !BS now I determine the size of the angular integral arrays
                jblock = 0
                do m1=-l1,l1
                  do m2=-l2,l2
                    do m3=-l3,l3
                      m4 = m1+m2-m3+1
                      if (abs(m4) <= l4) then
                        if ((.not. makemean) .or. ((l1 == l3) .and. (l2 == l4) .and. (abs(m2) == abs(m4))) .or. &
                            ((l1 == l2) .and. (l3 == l4) .and. ((abs(m1) == abs(m2)) .or. (abs(m3) == abs(m4))))) then
                          jblock = jblock+1
                        end if
                      end if
                    end do
                  end do
                end do
                do m1=0,l1
                  do m2=-l2,l2
                    do m3=-l3,l3
                      m4 = m1+m2-m3
                      if ((.not. makemean) .or. ((l1 == l3) .and. (l2 == l4) .and. (abs(m2) == abs(m4))) .or. &
                          ((l1 == l2) .and. (l3 == l4) .and. ((abs(m1) == abs(m2)) .or. (abs(m3) == abs(m4))))) then
                        if ((m1 /= 0) .or. (m2 /= 0) .or. (m3 /= 0)) then !  all m=0 make no sense
                          if (abs(m4) <= l4) then
                            jblock = jblock+1
                          end if
                        end if
                      end if
                    end do
                  end do
                end do
                !BS done !!
                !bs number of contracted integrals for each block
                ncont = ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
                mxangint = jblock*ncont
                !bs determine the size icont4 for the radial integrals
                call gencoulDIM(l1,l2,l3,l4,makemean,icont4)

                call mma_allocate(ANGSO,mxangint,Label='AngSO')
                call mma_allocate(ANGOO,mxangint,Label='AngOO')
                call mma_allocate(CartSO,nCont,Label='CartSO')
                call mma_allocate(CartOO,nCont,Label='CartOO')
                call mma_allocate(ConSO,iCont4,Label='ConSO')
                call mma_allocate(ConOO,iCont4,Label='ConOO')

                call gencoul(l1,l2,l3,l4,makemean,bonn,breit,sameorb,conSO,conOO,icont4,powexp,coulovlp)
                ! gen and trans integrals
                !bs local counter for integral adresses
                mblock = 0 ! counter of (m,m,m,m)-blocks for (l1,l2,l3,l4)
                !bs if keep is set to false, the angular integrals are
                !bs thrown away after each block of l-values
                !bs which means integrals start at address 0
                if (.not. keep) iangfirst = 0
                locstar = iangfirst ! local starting adress counter
                ! col 1 will hold type of integrals (1,2,3)
                ! col 2 will hold number of block
                mcombina(:,-l1:l1,-l2:l2,-l3:l3,-l4:l4) = 0
                do m1=-l1,l1
                  do m2=-l2,l2
                    do m3=-l3,l3
                      !bs m4 is more or less fixed by m1-3
                      !#######################################################################
                      !#######################################################################
                      !########## the L- -type block to be combined with sigma+ ##############
                      !#######################################################################
                      !#######################################################################
                      m4 = m1+m2-m3+1
                      if (abs(m4) <= l4) then !the L- -block to  combine with sigma+
                        !bs not all m-combinations are needed for the mean-field
                        if ((.not. makemean) .or. ((l1 == l3) .and. (l2 == l4) .and. (abs(m2) == abs(m4))) .or. &
                            ((l1 == l2) .and. (l3 == l4) .and. ((abs(m1) == abs(m2)) .or. (abs(m3) == abs(m4))))) then
                          mcombina(1,m1,m2,m3,m4) = 1
                          mblock = mblock+1
                          if (locstar+ncont > mxangint) then
                            write(u6,*) 'not enough space allocated for angular integrals'
                            write(u6,*) 'increase mxangint to at least ',locstar+ncont
                            call Abend()
                          end if
                          !bs mkangLmin = make_angular_integrals_for_L- type operator
                          !bs really generates  the angular prefactors and combines them with
                          !bs the radial integrals
                          call mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,AngSO(1+locstar),AngOO(1+locstar),Lfirst,Llast,Lblocks, &
                                         ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),ConSO(Lstarter(1)), &
                                         ConSO(Lstarter(2)),ConSO(Lstarter(3)),ConSO(Lstarter(4)),ConOO(Lstarter(1)), &
                                         ConOO(Lstarter(2)),ConOO(Lstarter(3)),ConOO(Lstarter(4)),preroots,clebsch,scratch,bonn, &
                                         breit,sameorb)
                          locstar = locstar+ncont ! increase starting address
                          mcombina(2,m1,m2,m3,m4) = mblock  ! set the block number
                          !#######################################################################
                          !#######################################################################
                          !########## the L+ -type block to be combined with sigma- ##############
                          !#######################################################################
                          !#######################################################################

                          ! these integrals are obtained by changing the signs of the m-values.
                          ! As the integrals are the same, the pointer points to the same integrals...

                          mcombina(1,-m1,-m2,-m3,-m4) = 3
                          mcombina(2,-m1,-m2,-m3,-m4) = mblock
                        end if
                      end if
                    end do
                  end do
                end do
                !#######################################################################
                !#######################################################################
                !########## the L0 -type block to be combined with sigma0 ##############
                !#######################################################################
                !#######################################################################
                do m1=0,l1
                  do m2=-l2,l2
                    do m3=-l3,l3
                      !bs m4 is more or less fixed by m1-3
                      m4 = m1+m2-m3 ! the L0-block to be combined with sigma0
                      !bs not all m-combinations are needed for the mean-field
                      if ((.not. makemean) .or. ((l1 == l3) .and. (l2 == l4) .and. (abs(m2) == abs(m4))) .or. &
                          ((l1 == l2) .and. (l3 == l4) .and. ((abs(m1) == abs(m2)) .or. (abs(m3) == abs(m4))))) then

                        if ((m1 /= 0) .or. (m2 /= 0) .or. (m3 /= 0)) then !all m=0 make no sense
                          if (abs(m4) <= l4) then
                            mcombina(1,m1,m2,m3,m4) = 2
                            mblock = mblock+1
                            if (locstar+ncont > mxangint) then
                              write(u6,*) 'not enough space allocated for angular integrals'
                              write(u6,*) 'increase mxangint to at least ',locstar+ncont
                              call Abend()
                            end if
                            call mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,angSO(1+locstar),AngOO(1+locstar),Lfirst,Llast,Lblocks, &
                                         ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),ConSO(Lstarter(1)), &
                                         ConSO(Lstarter(2)),ConSO(Lstarter(3)),ConSO(Lstarter(4)),ConOO(Lstarter(1)), &
                                         ConOO(Lstarter(2)),ConOO(Lstarter(3)),ConOO(Lstarter(4)),preroots,clebsch,scratch,bonn, &
                                         breit,sameorb)
                            locstar = locstar+ncont
                            mcombina(2,m1,m2,m3,m4) = mblock
                          end if
                        end if
                      end if
                    end do
                  end do
                end do
                !bs ###################################################################
                !bs ###################################################################
                !bs   transformation to l,m dependent integrals is finished
                !bs ###################################################################

                !bs ###################################################################
                !bs   begin transformation to cartesian integrals
                !bs ###################################################################
                !bs ###################################################################
                !bs check out, which combinations of m-values will
                !bs contribute to cartesian integrals
                ! col 1 will hold the type  x=1 y=2 z=3
                ! col 2 will hold the block number
                mcombcart(:,-l1:l1,-l2:l2,-l3:l3,-l4:l4) = 0
                mblockx = 0
                mblocky = 0
                mblockz = 0
                do m3=-l3,l3
                  do m4=-l4,l4
                    !bs if the l-values are the same : triangular matrix over m-values
                    !bs is sufficient
                    if (l1 == l3) then
                      m1upper = m3
                    else
                      m1upper = l1
                    end if
                    if (makemean) m1upper = l1
                    !bs if the l-values are the same : triangular matrix over m-values
                    !bs is sufficient
                    if (l2 == l4) then
                      m2upper = m4
                    else
                      m2upper = l2
                    end if
                    if (makemean) m2upper = l2
                    do m1=-l1,m1upper
                      if ((l1 == l3) .and. (m1 == m3)) then ! clean real zeros by symmetry
                        !bs this a problem of the spin-other-orbit integrals, as they are by
                        !bs formula not antisymmetric in the indices for particle 1.
                        cleaner = .true.
                      else
                        cleaner = .false.
                      end if
                      do m2=-l2,m2upper
                        !bs not all m-combinations are needed for the mean-field
                        if ((.not. makemean) .or. ((l1 == l3) .and. (l2 == l4) .and. (m2 == m4)) .or. &
                            ((l1 == l2) .and. (l3 == l4) .and. ((m1 == m2) .or. (m3 == m4)))) then

                          indx = ipowxyz(1,m1,l1)+ipowxyz(1,m2,l2)+ipowxyz(1,m3,l3)+ipowxyz(1,m4,l4)
                          indy = ipowxyz(2,m1,l1)+ipowxyz(2,m2,l2)+ipowxyz(2,m3,l3)+ipowxyz(2,m4,l4)
                          indz = ipowxyz(3,m1,l1)+ipowxyz(3,m2,l2)+ipowxyz(3,m3,l3)+ipowxyz(3,m4,l4)
                          indx = mod(indx,2)
                          indy = mod(indy,2)
                          indz = mod(indz,2)
                          if ((indx == 0) .and. (indy == 1) .and. (indz == 1) .and. &
                              (icheckxy(abs(m1),abs(m2),abs(m3),abs(m4)) > 0)) then
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !++++++++++++++++      SIGMA X      ++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            ! Y*Z ->  transforms like  L_x (B1)
                            !bs integrals for sigma_x
                            mblockx = mblockx+1
                            mcombcart(1,m1,m2,m3,m4) = 1
                            mcombcart(2,m1,m2,m3,m4) = mblockx
                            call tosigX(m1,m2,m3,m4,AngSO(1+iangfirst),mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3), &
                                        ncontrac(l4),CartSO,preXZ,interxyz(1,abs(m1),abs(m2),abs(m3),abs(m4)),isgnprod,cleaner)

                            if ((.not. bonn) .and. (.not. breit)) &
                              call tosigX(m1,m2,m3,m4,AngOO(1+iangfirst),mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3), &
                                          ncontrac(l4),cartOO,preXZ,interxyz(1,abs(m1),abs(m2),abs(m3),abs(m4)),isgnprod,cleaner)
                            if (makemean) then ! generate mean-field-contributions
                              !#######################################################################
                              !############  mean-field-part #########################################
                              !#######################################################################
                              if ((l1 == l3) .and. (l2 == l4)) then
                                if ((m2 == m4) .and. (m1 < m3) .and. (abs(m1+m3) == 1) .and. (l1 /= 0)) then
                                  call two2mean13(CartSO,occup(1,l2),AOcoeffs(1,1,l2),onecartx(1,1,iTri(m1+l1+1,m3+l3+1),l1), &
                                                  ncontrac(l1),ncontrac(l2),noccorb(l2))
                                end if
                              end if

                              if (NFINI) then
                                if ((l1 == l2) .and. (l3 == l4)) then
                                  if ((m1 == m2) .and. (l3 /= 0) .and. (l3 /= l1)) then
                                    if ((m3 < m4) .and. (abs(m4+m3) == 1)) then
                                      !bs for the "Bonn-approach"   exchange cartexOO by cartexSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean34a(cartSO,cartSO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartx(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      else
                                        if (NFINI) call two2mean34a(cartSO,cartOO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartx(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      end if
                                    end if
                                    if ((m3 > m4) .and. (abs(m4+m3) == 1)) then
                                      !bs for the "Bonn-approach"   exchange cartexOO by cartexSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean34b(CartSO,CartSO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartx(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      else
                                        if (NFINI) call two2mean34b(CartSO,CartOO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartx(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      end if
                                    end if
                                  else if ((m3 == m4) .and. (l1 /= 0)) then
                                    if (m1 < m2 .and. abs(m1+m2) == 1) then
                                      !bs for the "Bonn-approach"   exchange cartexOO by cartexSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean12a(CartSO,CartSO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartx(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      else
                                        if (NFINI) call two2mean12a(CartSO,cartOO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartx(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      end if
                                    end if
                                    if ((m1 > m2) .and. (abs(m1+m2) == 1)) then
                                      !bs for the "Bonn-approach"   exchange cartexOO by cartexSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean12b(cartSO,CartSO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartx(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      else
                                        if (NFINI) call two2mean12b(CartSO,CartOO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartx(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      end if
                                    end if
                                  end if
                                end if
                              end if ! If (NFINI) Then
                              !#######################################################################
                              !############  mean-field-part #########################################
                              !#######################################################################
                            end if
                          else if ((indx == 1) .and. (indy == 0) .and. (indz == 1) .and. &
                                   (icheckxy(abs(m1),abs(m2),abs(m3),abs(m4)) > 0)) then
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !++++++++++++++++      SIGMA Y      ++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            ! X*Z transforms like L_y  (B2)
                            !bs integrals for sigma_y
                            mblocky = mblocky+1
                            mcombcart(1,m1,m2,m3,m4) = 2
                            mcombcart(2,m1,m2,m3,m4) = mblocky
                            call tosigY(m1,m2,m3,m4,AngSO(1+iangfirst),mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3), &
                                        ncontrac(l4),CartSO,preY,interxyz(1,abs(m1),abs(m2),abs(m3),abs(m4)),isgnprod,cleaner)

                            if ((.not. bonn) .and. (.not. breit)) &
                              call tosigY(m1,m2,m3,m4,AngOO(1+iangfirst),mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3), &
                                          ncontrac(l4),cartOO,preY,interxyz(1,abs(m1),abs(m2),abs(m3),abs(m4)),isgnprod,cleaner)
                            if (makemean) then ! generate mean-field-contributions
                              !#######################################################################
                              !############  mean-field-part #########################################
                              !#######################################################################
                              if ((l1 == l3) .and. (l2 == l4)) then
                                if ((m2 == m4) .and. (m1 < m3) .and. (abs(m3-m1) == 1) .and. (l1 /= 0)) then
                                  call two2mean13(CartSO,occup(1,l2),AOcoeffs(1,1,l2),onecartY(1,1,iTri(m1+l1+1,m3+l3+1),l1), &
                                                  ncontrac(l1),ncontrac(l2),noccorb(l2))
                                end if
                              end if

                              if (NFINI) then
                                if ((l1 == l2) .and. (l3 == l4)) then
                                  if ((m1 == m2) .and. (l3 /= 0) .and. (l3 /= l1)) then
                                    if ((m3 < m4) .and. (abs(m3-m4) == 1)) then
                                      !bs for the "Bonn-approach"   exchange carteYOO by carteYSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean34a(CartSO,CartSO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartY(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      else
                                        if (NFINI) call two2mean34a(CartSO,CartOO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartY(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      end if
                                    end if
                                    if ((m3 > m4) .and. (abs(m3-m4) == 1)) then
                                      !bs for the "Bonn-approach"   exchange carteYOO by carteYSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean34b(CartSO,CartSO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartY(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      else
                                        if (NFINI) call two2mean34b(CartSO,CartOO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartY(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      end if
                                    end if
                                  else if ((m3 == m4) .and. (l1 /= 0)) then
                                    if ((m1 < m2) .and. (abs(m1-m2) == 1)) then
                                      !bs for the "Bonn-approach"   exchange carteOO by carteSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean12a(CartSO,CartSO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartY(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      else
                                        if (NFINI) call two2mean12a(CartSO,CartOO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartY(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      end if
                                    end if
                                    if ((m1 > m2) .and. (abs(m1-m2) == 1)) then
                                      !bs for the "Bonn-approach"   exchange carteYOO by carteYSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean12b(CartSO,CartSO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartY(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      else
                                        if (NFINI) call two2mean12b(CartSO,CartOO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartY(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      end if
                                    end if
                                  end if
                                end if
                              end if ! If (NFINI) Then
                              !#######################################################################
                              !############  mean-field-part #########################################
                              !#######################################################################
                            end if
                          else if ((indx == 1) .and. (indy == 1) .and. (indz == 0) .and. &
                                   (icheckz(abs(m1),abs(m2),abs(m3),abs(m4)) > 0)) then
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !++++++++++++++++      SIGMA Z      ++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            ! X*Y transforms like L_z  (A2)
                            !bs integrals for sigma_z
                            mblockz = mblockz+1
                            mcombcart(1,m1,m2,m3,m4) = 3
                            mcombcart(2,m1,m2,m3,m4) = mblockz
                            call tosigZ(m1,m2,m3,m4,angSO(1+iangfirst),mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3), &
                                        ncontrac(l4),CartSO,preXZ,interxyz(1,abs(m1),abs(m2),abs(m3),abs(m4)),isgnprod,cleaner)

                            if ((.not. bonn) .and. (.not. breit)) &
                              call tosigZ(m1,m2,m3,m4,AngOO(1+iangfirst),mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3), &
                                          ncontrac(l4),CartOO,preXZ,interxyz(1,abs(m1),abs(m2),abs(m3),abs(m4)),isgnprod,cleaner)
                            if (makemean) then ! generate mean-field-contributions
                              !#######################################################################
                              !############  mean-field-part #########################################
                              !#######################################################################
                              if ((l1 == l3) .and. (l2 == l4)) then
                                if ((m2 == m4) .and. (m1 < m3) .and. (m1 == -m3) .and. (l1 /= 0)) then
                                  call two2mean13(CartSO,occup(1,l2),AOcoeffs(1,1,l2),onecartz(1,1,iTri(m1+l1+1,m3+l3+1),l1), &
                                                  ncontrac(l1),ncontrac(l2),noccorb(l2))
                                end if
                              end if

                              if (NFINI) then
                                if ((l1 == l2) .and. (l3 == l4)) then
                                  if ((m1 == m2) .and. (l3 /= 0) .and. (l3 /= l1)) then
                                    if ((m3 < m4) .and. (m3 == -m4)) then
                                      !bs for the "Bonn-approach"   exchange carteOO by carteSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean34a(CartSO,CartSO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartz(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      else
                                        if (NFINI) call two2mean34a(CartSO,CartOO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartz(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      end if
                                    end if
                                    if ((m3 > m4) .and. (m3 == -m4)) then
                                      !bs for the "Bonn-approach"   exchange carteOO by carteSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean34b(CartSO,CartSO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartz(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      else
                                        if (NFINI) call two2mean34b(CartSO,CartOO,occup(1,l1),AOcoeffs(1,1,l1), &
                                                                    onecartz(1,1,iTri(m3+l3+1,m4+l4+1),l3),ncontrac(l3), &
                                                                    ncontrac(l1),noccorb(l2),sameorb)
                                      end if
                                    end if
                                  else if ((m3 == m4) .and. (l1 /= 0)) then
                                    if (m1 < m2 .and. m1 == -m2) then
                                      !bs for the "Bonn-approach"   exchange carteOO by carteSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean12a(CartSO,CartSO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartz(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      else
                                        if (NFINI) call two2mean12a(CartSO,CartOO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartz(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      end if
                                    end if
                                    if ((m1 > m2) .and. (m1 == -m2)) then
                                      !bs for the "Bonn-approach"   exchange carteOO by carteSO
                                      if (bonn .or. breit) then
                                        if (NFINI) call two2mean12b(cartSO,CartSO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartz(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      else
                                        if (NFINI) call two2mean12b(cartSO,cartOO,occup(1,l3),AOcoeffs(1,1,l3), &
                                                                    onecartz(1,1,iTri(m1+l1+1,m2+l2+1),l1),ncontrac(l1), &
                                                                    ncontrac(l3),noccorb(l3),sameorb)
                                      end if
                                    end if
                                  end if
                                end if
                              end if ! If (NFINI) Then
                              !#######################################################################
                              !############  mean-field-part #########################################
                              !#######################################################################
                            end if
                          end if
                          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end if ! for check of significance for meanfield.
                      end do
                    end do
                  end do
                end do
                numbcart = numbcart+(mblockx+mblocky+mblockz)*ncont
                !bs just controlling if x and y integrals have the same number of blocks
                if (mblockx /= mblocky) then
                  write(u6,*) 'numbers of integrals for sigma_x and sigma_y not equal!'
                  write(u6,'(A12,4I3,2(A3,I5))') 'l1,l2,l3,l4 ',l1,l2,l3,l4,' X:',mblockx,' Y:',mblocky
                  write(u6,*) ' check the ipowxyz-array'
                  call Abend()
                end if
                !bs start adresses for the next <ll|ll> block of integrals
                call mma_deallocate(AngSO)
                call mma_deallocate(AngOO)
                call mma_deallocate(CartSO)
                call mma_deallocate(CartOO)
                call mma_deallocate(ConSO)
                call mma_deallocate(ConOO)
              end if
            end if
          end if
        end if
      end do
    end do
  end do
end do
call mma_deallocate(clebsch)
call mma_deallocate(mcombina)
call mma_deallocate(mcombcart)

return

end subroutine angular
