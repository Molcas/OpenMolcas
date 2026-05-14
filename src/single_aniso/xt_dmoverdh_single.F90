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

subroutine XT_dMoverdH_single(nss,nTempMagn,nT,nM,Tmin,Tmax,XTexp,eso,T,zJ,Xfield,dM,sM,XT_no_field,tinput,smagn,mem,DoPlot)
! chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Nine, Ten, mBohr, rNAVO
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: nss, nTempMagn, nT, NM, mem
real(kind=wp), intent(in) :: Tmin, Tmax, XTexp(nT+nTempMagn), eso(nss), T(nT+nTempMagn), zJ, Xfield, XT_no_field(nT+nTempMagn)
complex(kind=wp), intent(in) :: dM(3,nss,nss), sM(3,nss,nss)
logical(kind=iwp), intent(in) :: tinput, smagn, doplot
integer(kind=iwp) :: iM, Info, iT, j, jT, mem_local, nDirX, nTempTotal
real(kind=wp) :: dHX(3), dHY(3), dHZ(3), hp, WT(3), Xfield_1, Xfield_2, Xfield_3, Xfield_4, Xfield_5, Xfield_6, Xfield_7, ZT(3,3)
character(len=50) :: label
real(kind=wp), allocatable :: MT1(:,:), MT2(:,:), MT3(:,:), MT4(:,:), MT5(:,:), MT6(:,:), MT7(:,:), ST1(:,:), ST2(:,:), ST3(:,:), &
                              ST4(:,:), ST5(:,:), ST6(:,:), ST7(:,:), WM1(:), WM2(:), WM3(:), WM4(:), WM5(:), WM6(:), WM7(:), &
                              XTM_dMdH(:), XTM_MH(:), XTtens_dMdH(:,:,:), XTtens_MH(:,:,:), ZT1(:), ZT2(:), ZT3(:), ZT4(:), &
                              ZT5(:), ZT6(:), ZT7(:)
real(kind=wp), parameter :: cm3tomB = rNAVO*mBohr/Ten, & ! in cm3 * mol-1 * T
                            THRS = 1.0e-13_wp
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, l
#  define _DBG_ .true.
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: DBG = _DBG_, m_paranoid = .true.
real(kind=wp), external :: dev

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(mem)
#endif

! threshold for convergence of average spin, in case (zJ /= 0)
mem_local = 0

!ccc-------------------------------------------------------cccc
#ifdef _DEBUGPRINT_
write(u6,'(A, 10I5)') '       nss =',nss
write(u6,'(A,   I5)') ' nTempMagn =',nTempMagn
write(u6,'(A,   I5)') '        nT =',nT
write(u6,'(A,   I5)') '        nM =',nM
write(u6,'(A,F12.5)') '      Tmin =',Tmin
write(u6,'(A,F12.5)') '      Tmax =',Tmax
write(u6,'(A,F12.5)') '        zJ =',zJ
write(u6,'(A,F12.5)') '    Xfield =',Xfield
write(u6,'(A,ES12.5)') '      THRS =',THRS
write(u6,*) '     smagn =',smagn
write(u6,*) '    tinput =',tinput
if (tinput) then
  do iT=1,nTempMagn
    write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
  end do
  do iT=1,nT
    jT = iT+nTempMagn
    write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',jT,')=',T(jT),' chiT_exp(',iT,')=',XTexp(jT)
  end do

else
  do iT=1,nT+nTempMagn
    write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
  end do
end if
write(u6,'(A)') 'SPIN-ORBIT ENERGY'
write(u6,'(10F12.5)') (ESO(i),i=1,NSS)
call xFlush(u6)
#endif
!ccc-------------------------------------------------------cccc
write(u6,*)
write(u6,'(A)') repeat('%',95)
write(u6,'(16X,A)') 'CALCULATION OF THE FIELD-DEPENDENT MAGNETIC SUSCEPTIBILITY'
write(u6,'(18X,A)') 'within true (dM/dH) and "experimentalists" (M/H) models'
write(u6,'(A)') repeat('%',95)
write(u6,*)
write(u6,'(2x,A,F10.6,A)') 'Magnetic field strength:',Xfield,' tesla.'
write(u6,'(2x,A,F10.6,A)') 'dM/dH is computed numerically using 7 point stencil formula, with h=0.00001'

if (tinput) then
  write(u6,'(2x,a)') 'Temperature dependence of the magnetic susceptibility and'
  write(u6,'(2x,a)') 'high-field magnetization will be calculated according to '
  write(u6,'(2x,a)') 'experimental values provided by the user in file "chitexp.input".'
else
  write(u6,'(2x,a,i3,a)') 'Temperature dependence of the magnetic susceptibility will be calculated in',nT,' points, '
  write(u6,'(2x,a,f4.1,a,f6.1,a)') 'equally distributed in temperature range ',tmin,' ---',tmax,' K.'
end if
! allocate memory:
call mma_allocate(WM1,nM,'WM1')
call mma_allocate(WM2,nM,'WM2')
call mma_allocate(WM3,nM,'WM3')
call mma_allocate(WM4,nM,'WM4')
call mma_allocate(WM5,nM,'WM5')
call mma_allocate(WM6,nM,'WM6')
call mma_allocate(WM7,nM,'WM7')
mem_local = mem_local+size(WM1)*RtoB
mem_local = mem_local+size(WM2)*RtoB
mem_local = mem_local+size(WM3)*RtoB
mem_local = mem_local+size(WM4)*RtoB
mem_local = mem_local+size(WM5)*RtoB
mem_local = mem_local+size(WM6)*RtoB
mem_local = mem_local+size(WM7)*RtoB

call mma_allocate(ZT1,nT+nTempMagn,'ZT1')
call mma_allocate(ZT2,nT+nTempMagn,'ZT2')
call mma_allocate(ZT3,nT+nTempMagn,'ZT3')
call mma_allocate(ZT4,nT+nTempMagn,'ZT4')
call mma_allocate(ZT5,nT+nTempMagn,'ZT5')
call mma_allocate(ZT6,nT+nTempMagn,'ZT6')
call mma_allocate(ZT7,nT+nTempMagn,'ZT7')
mem_local = mem_local+size(ZT1)*RtoB
mem_local = mem_local+size(ZT2)*RtoB
mem_local = mem_local+size(ZT3)*RtoB
mem_local = mem_local+size(ZT4)*RtoB
mem_local = mem_local+size(ZT5)*RtoB
mem_local = mem_local+size(ZT6)*RtoB
mem_local = mem_local+size(ZT7)*RtoB

call mma_allocate(MT1,3,nT+nTempMagn,'MT0')
call mma_allocate(MT2,3,nT+nTempMagn,'MT1')
call mma_allocate(MT3,3,nT+nTempMagn,'MT3')
call mma_allocate(MT4,3,nT+nTempMagn,'MT4')
call mma_allocate(MT5,3,nT+nTempMagn,'MT5')
call mma_allocate(MT6,3,nT+nTempMagn,'MT6')
call mma_allocate(MT7,3,nT+nTempMagn,'MT7')
mem_local = mem_local+size(MT1)*RtoB
mem_local = mem_local+size(MT2)*RtoB
mem_local = mem_local+size(MT3)*RtoB
mem_local = mem_local+size(MT4)*RtoB
mem_local = mem_local+size(MT5)*RtoB
mem_local = mem_local+size(MT6)*RtoB
mem_local = mem_local+size(MT7)*RtoB

call mma_allocate(ST1,3,nT+nTempMagn,'ST1')
call mma_allocate(ST2,3,nT+nTempMagn,'ST2')
call mma_allocate(ST3,3,nT+nTempMagn,'ST3')
call mma_allocate(ST4,3,nT+nTempMagn,'ST4')
call mma_allocate(ST5,3,nT+nTempMagn,'ST5')
call mma_allocate(ST6,3,nT+nTempMagn,'ST6')
call mma_allocate(ST7,3,nT+nTempMagn,'ST7')
mem_local = mem_local+size(ST1)*RtoB
mem_local = mem_local+size(ST2)*RtoB
mem_local = mem_local+size(ST3)*RtoB
mem_local = mem_local+size(ST4)*RtoB
mem_local = mem_local+size(ST5)*RtoB
mem_local = mem_local+size(ST6)*RtoB
mem_local = mem_local+size(ST7)*RtoB

call mma_allocate(XTM_MH,nT+nTempMagn,'XTM_MH')
call mma_allocate(XTM_dMdH,nT+nTempMagn,'XTM_dMdH')
call mma_allocate(XTtens_MH,3,3,nT+nTempMagn,'XTtens_MH')
call mma_allocate(XTtens_dMdH,3,3,nT+nTempMagn,'XTtens_dMdH')
mem_local = mem_local+size(XTM_MH)*RtoB
mem_local = mem_local+size(XTM_dMdH)*RtoB
mem_local = mem_local+size(XTtens_MH)*RtoB
mem_local = mem_local+size(XTtens_dMdH)*RtoB
#ifdef _DEBUGPRINT_
write(u6,*) 'XTMG:  memory allocated (local):'
write(u6,*) 'mem_local=',mem_local
write(u6,*) 'XTMG:  memory allocated (total):'
write(u6,*) 'mem_total=',mem+mem_local
#endif

nDirX = 3

! field along X
dHX(1) = One
dHY(1) = Zero
dHZ(1) = Zero
! field along Y
dHX(2) = Zero
dHY(2) = One
dHZ(2) = Zero
! field along Z
dHX(3) = Zero
dHY(3) = Zero
dHZ(3) = One

XTtens_MH(:,:,:) = Zero
XTtens_dMdH(:,:,:) = Zero

hp = 0.0001_wp
Xfield_1 = Xfield-Three*hp
Xfield_2 = Xfield-Two*hp
Xfield_3 = Xfield-hp
Xfield_4 = xField
Xfield_5 = Xfield+hp
Xfield_6 = Xfield+Two*hp
Xfield_7 = Xfield+Three*hp
#ifdef _DEBUGPRINT_
write(u6,*) 'XTMG:  Xfield: ',Xfield_1,Xfield_2,Xfield_3,Xfield_4,Xfield_5,Xfield_6,Xfield_7
#endif

nTempTotal = nT+nTempMagn

! ///  opening the loop over different directions of the magnetic field
do iM=1,nDirX

  ! compute magnetization:
  ! seven points numerical field perturbation:
  ! Xfield =XF -3h
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_1,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM1,ZT1,ST1,MT1,m_paranoid,DBG)
  ! Xfield =XF -2h
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_2,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM2,ZT2,ST2,MT2,m_paranoid,DBG)
  ! Xfield =XF -h
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_3,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM3,ZT3,ST3,MT3,m_paranoid,DBG)
  ! Xfield =XF
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_4,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM4,ZT4,ST4,MT4,m_paranoid,DBG)
  ! Xfield =XF +h
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_5,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM5,ZT5,ST5,MT5,m_paranoid,DBG)
  ! Xfield =XF +2h
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_6,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM6,ZT6,ST6,MT6,m_paranoid,DBG)
  ! Xfield =XF +3h
  call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),XField_7,ESO,zJ,THRS,DM,SM,nTempTotal,T,smagn,WM7,ZT7,ST7,MT7,m_paranoid,DBG)

# ifdef _DEBUGPRINT_
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex1(',iT,')=',(MT1(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex2(',iT,')=',(MT2(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex3(',iT,')=',(MT3(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex4(',iT,')=',(MT4(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex5(',iT,')=',(MT5(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex6(',iT,')=',(MT6(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex7(',iT,')=',(MT7(l,iT),l=1,3)
  end do
  do iT=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex 7-1 (',iT,')=',(MT7(l,iT)-MT1(l,iT),l=1,3)
  end do
# endif
  ! computing the AVERAGE MOMENTS calculated at different temperatures (T(i))
  do iT=1,nTempTotal

    ! dM/dH model
    if (iM == 1) then
      XTtens_dMdH(iM,1,iT) = (-MT1(1,iT)+Nine*MT2(1,iT)-45.0_wp*MT3(1,iT)+45.0_wp*MT5(1,iT)-Nine*MT6(1,iT)+MT7(1,iT))*T(iT)* &
                             cm3tomB/(60.0_wp*hp)
    else if (iM == 2) then
      XTtens_dMdH(iM,2,iT) = (-MT1(2,iT)+Nine*MT2(2,iT)-45.0_wp*MT3(2,iT)+45.0_wp*MT5(2,iT)-Nine*MT6(2,iT)+MT7(2,iT))*T(iT)* &
                             cm3tomB/(60.0_wp*hp)
    else if (iM == 3) then
      XTtens_dMdH(iM,3,iT) = (-MT1(3,iT)+Nine*MT2(3,iT)-45.0_wp*MT3(3,iT)+45.0_wp*MT5(3,iT)-Nine*MT6(3,iT)+MT7(3,iT))*T(iT)* &
                             cm3tomB/(60.0_wp*hp)
    end if

    ! M/H model
    if (iM == 1) then
      XTtens_MH(iM,1,iT) = MT4(1,iT)*T(iT)*cm3tomB/Xfield
    else if (iM == 2) then
      XTtens_MH(iM,2,iT) = MT4(2,iT)*T(iT)*cm3tomB/Xfield
    else if (iM == 3) then
      XTtens_MH(iM,3,iT) = MT4(3,iT)*T(iT)*cm3tomB/Xfield
    end if
  end do
! ///  closing the loops over field directions
end do ! iM (nDirX)

! computing the XT as tensor's average:
XTM_dMdH(:) = (XTtens_dMdH(1,1,:)+XTtens_dMdH(2,2,:)+XTtens_dMdH(3,3,:))/Three
XTM_MH(:) = (XTtens_MH(1,1,:)+XTtens_MH(2,2,:)+XTtens_MH(3,3,:))/Three
call Add_Info('T_dMdH            ',T,nTempTotal,5)
call Add_Info('XTM_dMdH          ',XTM_dMdH,nTempTotal,5)
call Add_Info('XTM_MH            ',XTM_MH,nTempTotal,5)

!-----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
!-----------------------------------------------------------------------
write(u6,*)
call xFlush(u6)
write(u6,'(A)') '--------------------------------------------------------------------------|'
write(u6,'(A)') '     |     T      | Statistical |   CHI*T     |   CHI*T     |   CHI*T     |'
write(u6,'(A,F6.3,A,F6.3,A)') '     |            |  Sum (Z)    | H = 0.000 T | H =',Xfield,' T | H =',XField,' T |'
write(u6,'(A)') '     |            |             |             | X = dM/dH   | X = M/H     |'
write(u6,'(A)') '-----|--------------------------------------------------------------------|'
write(u6,'(A)') 'Units|   kelvin   |    ---      |  cm3*K/mol  |  cm3*K/mol  |  cm3*K/mol  |'
write(u6,'(A)') '-----|--------------------------------------------------------------------|'

call xFlush(u6)
do iT=1,nT
  jT = iT+nTempMagn
  write(u6,'(A,F11.6,A,ES12.5,A,F12.8,A,F12.8,A,F12.7,A)') '     |',T(jT),' |',ZT3(jT),' |',XT_no_field(jT),' |',XTM_dMdH(jT), &
                                                           ' |',XTM_MH(jT),' |'
end do
write(u6,'(A)') '-----|--------------------------------------------------------------------|'
call xFlush(u6)

! calcualtion of the standard deviation:
if (tinput) then
  write(u6,'(a,5x, f20.14)') 'ST.DEV: X= dM/dH:',dev(nT,XTM_dMdH(1+nTempMagn:),XTexp(1+nTempMagn:))
  write(u6,'(a,5x, f20.14)') 'ST.DEV: X= M/H:',dev(nT,XTM_MH(1+nTempMagn:),XTexp(1+nTempMagn:))
  write(u6,'(A)') '-----|--------------------------------------------------------------------|'
end if !tinput

!-------------------------  PLOTs -------------------------------------!
write(label,'(A)') 'with_field_M_over_H'
call xFlush(u6)
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T(1+nTempMagn:),XTM_MH(1+nTempMagn:),XTexp(1+nTempMagn:))
  else
    call plot_XT_no_Exp(label,nT,T(1+nTempMagn:),XTM_MH(nTempMagn:))
  end if
end if

write(label,'(A)') 'with_field_dM_over_dH'
call xFlush(u6)
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T(1+nTempMagn:),XTM_dMdH(1+nTempMagn:),XTexp(1+nTempMagn:))
  else
    call plot_XT_no_Exp(label,nT,T(1+nTempMagn:),XTM_dMdH(1+nTempMagn:))
  end if
end if
!---------------------- END PLOTs -------------------------------------!

! print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
write(u6,'(/)')
write(u6,'(2A)') repeat('-',110),'|'
write(u6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR the  X=dM/dH  model,  in cm3*K/mol','|'
write(u6,'(2A)') repeat('-',110),'|'
write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
do iT=1,nT
  jT = iT+nTempMagn
  info = 0
  call DIAG_R2(XTtens_dMdH(:,:,jT),3,info,wt,zt)
  write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- z ----|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(XTtens_dMdH(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                  (zt(j,1),j=1,3),'|'
  write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(XTtens_dMdH(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                           (zt(j,2),j=1,3),'|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(XTtens_dMdH(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                  (zt(j,3),j=1,3),'|'
end do
write(u6,'(2A)') repeat('-',110),'|'

write(u6,'(/)')
write(u6,'(2A)') repeat('-',110),'|'

write(u6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR the  X = M/H  model,  in cm3*K/mol','|'
write(u6,'(2A)') repeat('-',110),'|'
write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
do iT=1,nT
  jT = iT+nTempMagn
  info = 0
  call DIAG_R2(XTtens_MH(:,:,jT),3,info,wt,zt)
  write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- z ----|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(XTtens_MH(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                  (zt(j,1),j=1,3),'|'
  write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(XTtens_MH(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                           (zt(j,2),j=1,3),'|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(XTtens_MH(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                  (zt(j,3),j=1,3),'|'
end do
write(u6,'(2A)') repeat('-',110),'|'

call mma_deallocate(WM1)
call mma_deallocate(WM2)
call mma_deallocate(WM3)
call mma_deallocate(WM4)
call mma_deallocate(WM5)
call mma_deallocate(WM6)
call mma_deallocate(WM7)

call mma_deallocate(ZT1)
call mma_deallocate(ZT2)
call mma_deallocate(ZT3)
call mma_deallocate(ZT4)
call mma_deallocate(ZT5)
call mma_deallocate(ZT6)
call mma_deallocate(ZT7)

call mma_deallocate(MT1)
call mma_deallocate(MT2)
call mma_deallocate(MT3)
call mma_deallocate(MT4)
call mma_deallocate(MT5)
call mma_deallocate(MT6)
call mma_deallocate(MT7)

call mma_deallocate(ST1)
call mma_deallocate(ST2)
call mma_deallocate(ST3)
call mma_deallocate(ST4)
call mma_deallocate(ST5)
call mma_deallocate(ST6)
call mma_deallocate(ST7)

call mma_deallocate(XTM_MH)
call mma_deallocate(XTM_dMdH)
call mma_deallocate(XTtens_MH)
call mma_deallocate(XTtens_dMdH)

return

end subroutine XT_dMoverdH_single
