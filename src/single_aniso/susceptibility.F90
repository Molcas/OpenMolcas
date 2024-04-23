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

subroutine SUSCEPTIBILITY(NSS,ESO,S_SO,DIPSO,nT,nTempMagn,T,tmin,tmax,XTexp,zJ,tinput,chiT_theta,doplot,iPrint,mem)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! this routine calculates the magnetic susceptibility and all related to it values.
! the units are cgsemu: cm^3*k/mol
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, Ten, cm_s, hPlanck, kBoltzmann, mBohr, rNAVO
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: nss, nT, nTempMagn, iprint, mem
real(kind=wp), intent(in) :: eso(nss), T(nT+nTempMagn), tmin, tmax, XTexp(nT+nTempMagn), zJ
complex(kind=wp), intent(in) :: s_so(3,nss,nss), dipso(3,nss,nss)
logical(kind=iwp), intent(in) :: tinput, doplot
real(kind=wp), intent(out) :: chit_theta(nT+nTempMagn)
integer(kind=iwp) :: ic, im, info, it, j, jc, jm, jT, mem_local
real(kind=wp) :: a_dir(3,3), a_inv(3,3), det, gtens(3), maxes(3,3), WT(3), XMM(3,3), XSM(3,3), XSS(3,3), xxm, XZJ(3,3), Zst, ZT(3,3)
character(len=50) :: label
real(kind=wp), allocatable :: chi_theta_1(:), chit(:), chit_tens(:,:,:), chit_theta_tens(:,:,:), zstat1(:)
real(kind=wp), parameter :: coeff_X = rNAVO*mBohr**2/kBoltzmann/Ten, &
                            boltz_k = kBoltzmann/(cm_s*hPlanck) ! boltzmann constant in cm-1*K-1
real(kind=wp), external :: dev

! constants used in this subrutine
mem_local = 0

write(u6,*)
write(u6,'(A)') repeat('%',95)
write(u6,'(35X,A)') 'CALCULATION OF THE MAGNETIC SUSCEPTIBILITY'
write(u6,'(A)') repeat('%',95)

write(u6,*)
if (TINPUT) then
  write(u6,'(5X,A)') 'Temperature dependence of the magnetic susceptibility calculated according '
  write(u6,'(5X,A)') 'to experimental values provided in the input'
else
  write(u6,'(5X,A)') 'Temperature dependence of the magnetic susceptibility calculated in'
  write(u6,'(5x,I3,a,f4.1,a,f6.1,a)') nT,' points, equally distributed in temperature range',Tmin,' ---',Tmax,' K.'
end if
write(u6,'(5X,A)') 'The algorithm employed for XT=f(T) in this section is based on the zero magnetic field limit.'
if (doplot) write(u6,'(5X,A)') 'The GNUPLOT script and correponding images are generated in $WorkDir'
write(u6,*)
!-----------------------------------------------------------------------
if (iPrint > 2) then
  write(u6,*) 'Tmin       =',tmin
  write(u6,*) 'Tmax       =',tmax
  write(u6,*) '  nT       =',nT
  write(u6,*) '  nTempMagn=',nTempMagn
  write(u6,*) 'Temperature:',T(:)
  write(u6,*) 'ESO:',ESO(:)
  call prMom('SUSC: input  S_SO(l,i,j):',S_SO,nss)
  call prMom('SUSC: input DIPSO(l,i,j):',DIPSO,nss)
  call atens(DIPSO,nss,gtens,maxes,2)
end if
!-----------------------------------------------------------------------
! ***********************************************************
! *     Powder averaged high-T limit of susceptibility      *
! *                   and g-tensor                          *
! ***********************************************************
!ee = Zero
!ee_2 = Zero
!do Iss=1,NSS
!  do Jss=1,NSS
!    do ic=1,3
!      ee = ee+real(conjg(DipSO(ic,Iss,Jss))*DipSO(ic,Iss,Jss))
!      if ((Iss <= 4) .and. (Jss <=4)) ee_2 = ee_2+real(conjg(DipSO(ic,Iss,Jss))*DipSO(ic,Iss,Jss))
!    end do ! ic
!  end do ! Jss
!end do ! Iss
!chiT_high = coeff_X*ee/real(NSS,kind=wp)/Three
!SS1 = real(IGSM-1)**2*Quart+real(IGSM-1)*Half
!chiT_high_2 = coeff_X*ee_2/Twelve
!g_high = sqrt(chiT_high/(coeff_X*SS1/Three))
!g_high_2 = sqrt(chiT_high_2/(coeff_X*SS1/Three))
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
call mma_allocate(chiT,nT+nTempMagn,'chiT')
call mma_allocate(Zstat1,nT+nTempMagn,'Zstat1')
call mma_allocate(chi_theta_1,nT+nTempMagn,'chi_theta_1')
call mma_allocate(chiT_tens,3,3,nT+nTempMagn,'chiT_tens')
mem_local = mem_local+size(chiT)*RtoB
mem_local = mem_local+size(Zstat1)*RtoB
mem_local = mem_local+size(chi_theta_1)*RtoB
mem_local = mem_local+size(chiT_tens)*RtoB

if (zJ == 0) then
  if (iPrint > 2) then
    write(u6,*) 'SUSC:  zJ = 0'
    write(u6,*) 'SUSC:  memory allocated (local):'
    write(u6,*) 'mem_local=',mem_local
    write(u6,*) 'SUSC:  memory allocated (total):'
    write(u6,*) 'mem_total=',mem+mem_local
  end if

  do iT=1,nT+nTempMagn
    ! compute XT tensor for this temperature:
    call chi(DipSO,DipSO,Eso,Nss,T(iT),Zst,XMM)
    if (iPrint > 2) then
      write(u6,'(A,9F12.6)') 'XMM:',XMM(:,:)
      write(u6,'(A,9F12.6)') 'chiT:',coeff_X*(XMM(1,1)+XMM(2,2)+XMM(3,3))/Three,Zst
    end if

    !XMM(:,:) = coeff_X*XMM(:,:)
    chiT_tens(:,:,iT) = coeff_X*XMM(:,:)

    !chiT_tens(:,:,iT) = XMM(:,:)
    ! compute the powder XT for this temperature:
    chiT(iT) = coeff_X*(XMM(1,1)+XMM(2,2)+XMM(3,3))/Three

    if (abs(chit(iT)) < 1.0e-20_wp) chit(iT) = 1.0e-20_wp

    chit_theta(iT) = chiT(iT)
    chi_theta_1(iT) = T(iT)/chiT(iT)
    Zstat1(iT) = Zst
  end do
else  ! i.e. when zJ /= Zero
  if (iPrint > 2) write(u6,*) 'SUSC:  zJ \= 0'
  ! allocate matrices:
  call mma_allocate(chiT_theta_tens,3,3,nT+nTempMagn,'XTT')
  ! initialize:
  mem_local = mem_local+size(chiT_theta_tens)*RtoB
  if (iPrint > 2) then
    write(u6,*) 'SUSC:  memory allocated (local):'
    write(u6,*) 'mem_local=',mem_local
    write(u6,*) 'SUSC:  memory allocated (total):'
    write(u6,*) 'mem_total=',mem+mem_local
  end if

  do iT=1,nT+nTempMagn
    ! compute tensors:
    call chi(DipSO,DipSO,Eso,Nss,T(it),Zst,XMM)
    call chi(S_SO,DipSO,Eso,Nss,T(it),Zst,XSM)
    call chi(S_SO,S_SO,Eso,Nss,T(it),Zst,XSS)
    ! form the A_dir matrix:
    ! A_dir(:,:) = 1(:,:)*kB*T(iT)-zJ*XSS(:,:)
    call unitmat(A_dir,3)
    A_dir(:,:) = Boltz_k*T(iT)*A_dir(:,:)-zJ*XSS(:,:)
    ! invert it:
    call REVERSE(A_dir,A_inv,DET)

    do ic=1,3
      do jc=1,3
        ! compute the trace of the product XSM*A*XSM
        xxm = Zero
        do im=1,3
          do jm=1,3
            xxm = xxm+XSM(im,ic)*A_inv(im,jm)*XSM(jm,jc)
          end do
        end do
        ! add this contribution to the total susceptibility tensor
        XZJ(ic,jc) = XMM(ic,jc)+zJ*xxm
      end do ! jc
    end do ! ic

    ! Scale the tensors by coeff_X factor:
    XMM(:,:) = coeff_X*XMM(:,:)
    XZJ(:,:) = coeff_X*XZJ(:,:)
    ! place the tensors in the corresponding part of the "big" arrays:
    chiT_tens(:,:,iT) = XMM(:,:)
    chiT_theta_tens(:,:,iT) = XZJ(:,:)
    ! compute powder:
    chiT(iT) = (XMM(1,1)+XMM(2,2)+XMM(3,3))/Three
    chiT_theta(iT) = (XZJ(1,1)+XZJ(2,2)+XZJ(3,3))/Three

    if (abs(chit(iT)) < 1.0e-20_wp) chit(iT) = 1.0e-20_wp
    if (abs(chiT_theta(iT)) < 1.0e-20_wp) chiT_theta(iT) = 1.0e-20_wp

    chi_theta_1(iT) = T(iT)/chiT_theta(iT)
    Zstat1(iT) = Zst
  end do ! iT
end if ! zJ

!  WRITING SOME OF THE OUTPUT...

write(u6,'(/)')
write(u6,'(A)') '     |     T      | Statistical |    CHI*T    |    CHI*T    |     CHI     |    1/CHI    |'
write(u6,'(A)') '     |            |  Sum (Z)    |    (zJ=0)   |             |             |             |'
write(u6,'(A)') '     |            |             |             |             |             |             |'
write(u6,'(A)') '-----|----------------------------------------------------------------------------------|'
write(u6,'(A)') 'Units|   kelvin   |    ---      |  cm3*K/mol  |  cm3*K/mol  |   cm3/mol   |   mol/cm3   |'
write(u6,'(A)') '-----|----------------------------------------------------------------------------------|'

do iT=1,nT
  jT = iT+nTempMagn
  write(u6,'(A,F11.6,A,ES12.5,A,F12.8,A,F12.8,A,ES12.5,A,ES12.5,A)') '     |',T(jT),' |',zstat1(jT),' |',chiT(jT),' |', &
                                                                     chiT_theta(jT),' |',chiT_theta(jT)/T(jT),' |', &
                                                                     chi_theta_1(jT),' |'
end do
write(u6,'(A)') '-----|----------------------------------------------------------------------------------|'

if (TINPUT) then
  write(u6,'(/)')
  write(u6,'(5X,A      )') 'STANDARD DEVIATION OF THE CALCULATED MAGNETIC SUSCEPTIBILITY'
  write(u6,'(5X,A,F12.7)') 'FROM EXPERIMENTAL VALUES PROVIDED IN THE INPUT FILE IS:', &
                           dev(nT,chit_theta(1+nTempMagn:),XTexp(1+nTempMagn:))

end if

!-------------------------  PLOTs -------------------------------------!
write(label,'(A)') 'no_field'
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T(1+nTempMagn:),chit_theta(1+nTempMagn:),XTexp(1+nTempMagn:))
  else
    call plot_XT_no_Exp(label,nT,T(1+nTempMagn:),chit_theta(1+nTempMagn:))
  end if
end if
!---------------------- END PLOTs -------------------------------------!

! print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
if (zJ == Zero) then

  write(u6,'(/)')
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(31X,A,22x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ = 0,  in cm3*K/mol','|'
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(A   )') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |'// &
                     '               Main Axes             |'
  do iT=1,nT
    jT = iT+nTempMagn
    info = 0
    call DIAG_R2(chiT_tens(:,:,jT),3,info,wt,zt)
    write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y ---------'// &
                    ' z ----|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(chiT_tens(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                    (zt(j,1),j=1,3),'|'
    write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(chiT_tens(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                             (zt(j,2),j=1,3),'|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(chiT_tens(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                    (zt(j,3),j=1,3),'|'
  end do
  write(u6,'(2A)') repeat('-',110),'|'

else ! zJ /= Zero

  write(u6,'(/)')
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(31X,A,F7.4,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ = ',zJ,',  in cm3*K/mol','|'
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
  do iT=1,nT
    jT = iT+nTempMagn
    info = 0
    call DIAG_R2(chiT_theta_tens(:,:,jT),3,info,wt,zt)
    write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y ---------'// &
                    ' z ----|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(chiT_theta_tens(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                    (zt(j,1),j=1,3),'|'
    write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(chiT_theta_tens(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                             (zt(j,2),j=1,3),'|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(chiT_theta_tens(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                    (zt(j,3),j=1,3),'|'
  end do
  write(u6,'(2A)') repeat('-',110),'|'
  call mma_deallocate(chiT_theta_tens)
end if ! zJ

!write(u6,'(/)')
!write(u6,'(10A)') (('------------'),K=1,10)
!write(u6,'(5X,A)') 'MAGNETIC SUSCEPTIBILITY IN THE DIRECTION OF THE MAIN MAGNETIC AXES'
!write(u6,'(10A)') (('------------'),K=1,10)
!write(u6,*)
!write(u6,'(7X,A)') 'T(K)         gx          gy          gz'
!write(u6,*)
!do iT=1,nT
!  write(u6,'(4X,F6.1,6X,3(F8.4,4X))') T(iT),(ChiT_main(iT,ic),ic=1,3)
!end do
! saving some information for tests:

call Add_Info('Temperature',T,nT+nTempMagn,5)
call Add_Info('CHIT',chiT,nT+nTempMagn,5)
call Add_Info('CHIT_THETA',chiT_theta,nT+nTempMagn,5)

call mma_deallocate(Zstat1)
call mma_deallocate(chiT_tens)
call mma_deallocate(chiT)
call mma_deallocate(chi_theta_1)

return

end subroutine SUSCEPTIBILITY
