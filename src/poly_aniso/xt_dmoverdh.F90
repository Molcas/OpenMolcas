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

subroutine XT_dMoverdH(exch,nLoc,nCenter,nneq,neqv,neq,nss,nexch,nTempMagn,nT,NM,iopt,mem,Tmin,Tmax,chit_exp,eso,w,T,RROT,zJ, &
                       Xfield,EM,THRS,XT_no_field,dipso,s_so,dipexch,s_exch,tinput,smagn,m_paranoid,m_accurate,doplot)
! this Subroutine computes the XT as M/H as it is observed for most of experiments.
! The M is averaged over the grid for each temperature point.
!       chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]

use Constants, only: Zero, One, Three, Ten, mBohr, rNAVO
use Definitions, only: wp, u6

implicit none
#include "mgrid.fh"
#include "stdalloc.fh"
integer, intent(in) :: exch, nLoc, nCenter, nneq, neqv
integer, intent(in) :: nTempMagn, nT, NM, iopt, mem
integer, intent(in) :: neq(nneq), nss(nneq), nexch(nneq)
real(kind=8), intent(in) :: Tmin, Tmax
real(kind=8), intent(in) :: eso(nneq,nLoc)
real(kind=8), intent(in) :: W(exch)
real(kind=8), intent(in) :: zJ
real(kind=8), intent(in) :: RROT(nneq,neqv,3,3)
real(kind=8), intent(in) :: T(nT+nTempMagn)
real(kind=8), intent(in) :: chit_exp(nT)
real(kind=8), intent(in) :: XT_no_field(nT+nTempMagn)
real(kind=8), intent(in) :: Xfield
real(kind=8), intent(in) :: EM
real(kind=8), intent(in) :: THRS
complex(kind=8), intent(in) :: dipso(nneq,3,nLoc,nLoc)
complex(kind=8), intent(in) :: s_so(nneq,3,nLoc,nLoc)
complex(kind=8), intent(in) :: dipexch(3,exch,exch)
complex(kind=8), intent(in) :: s_exch(3,exch,exch)
logical, intent(in) :: tinput, smagn, m_paranoid, doplot
!ccc local variables ccc
integer :: nDirX
integer :: iM, iT, jT, i, j, l, n, isite !,nP
real(kind=8), allocatable :: WEX0(:)
real(kind=8), allocatable :: WEX1(:)
real(kind=8), allocatable :: WEX2(:)   ! Zeeman exchange energies
real(kind=8), allocatable :: WL0(:,:)    ! Zeeman local energies
real(kind=8), allocatable :: WL1(:,:)
real(kind=8), allocatable :: WL2(:,:)    ! Zeeman local energies
! Zeeman local reduced energies, using only NEXCH states;
real(kind=8), allocatable :: WR0(:,:)
real(kind=8), allocatable :: WR1(:,:)
real(kind=8), allocatable :: WR2(:,:)
! local statistical sum, Boltzmann distribution
real(kind=8), allocatable :: ZL0(:,:)
real(kind=8), allocatable :: ZL1(:,:)
real(kind=8), allocatable :: ZL2(:,:)
! local statistical sum, Boltzmann distribution, using only NEXCH states
real(kind=8), allocatable :: ZR0(:,:)
real(kind=8), allocatable :: ZR1(:,:)
real(kind=8), allocatable :: ZR2(:,:)
! spin magnetisation, from the local sites, using ALL states
real(kind=8), allocatable :: SL0(:,:,:)
real(kind=8), allocatable :: SL1(:,:,:)
real(kind=8), allocatable :: SL2(:,:,:)
! spin magnetisation, from the local sites, using only NEXCH states
real(kind=8), allocatable :: SR0(:,:,:)
real(kind=8), allocatable :: SR1(:,:,:)
real(kind=8), allocatable :: SR2(:,:,:)
! magnetisation, from local sites, using ALL states
real(kind=8), allocatable :: ML0(:,:,:)
real(kind=8), allocatable :: ML1(:,:,:)
real(kind=8), allocatable :: ML2(:,:,:)
! magnetisation, from local sites, using only NEXCH states
real(kind=8), allocatable :: MR0(:,:,:)
real(kind=8), allocatable :: MR1(:,:,:)
real(kind=8), allocatable :: MR2(:,:,:)
! spin magnetisation, from the exchange block
real(kind=8), allocatable :: SEX0(:,:)
real(kind=8), allocatable :: SEX1(:,:)
real(kind=8), allocatable :: SEX2(:,:)
! magnetisation, form the exchange block
real(kind=8), allocatable :: MEX0(:,:)
real(kind=8), allocatable :: MEX1(:,:)
real(kind=8), allocatable :: MEX2(:,:)
! exchange statistical sum, Boltzmann distribution
real(kind=8), allocatable :: ZEX0(:)
real(kind=8), allocatable :: ZEX1(:)
real(kind=8), allocatable :: ZEX2(:)
! total vectors in general coordinate system:
real(kind=8), allocatable :: ZRT0(:,:)
real(kind=8), allocatable :: ZLT0(:,:)
real(kind=8), allocatable :: ZRT1(:,:)
real(kind=8), allocatable :: ZLT1(:,:)
real(kind=8), allocatable :: ZRT2(:,:)
real(kind=8), allocatable :: ZLT2(:,:)
real(kind=8), allocatable :: MRT0(:,:,:)
real(kind=8), allocatable :: MLT0(:,:,:)
real(kind=8), allocatable :: SRT0(:,:,:)
real(kind=8), allocatable :: SLT0(:,:,:)
real(kind=8), allocatable :: MRT1(:,:,:)
real(kind=8), allocatable :: MLT1(:,:,:)
real(kind=8), allocatable :: SRT1(:,:,:)
real(kind=8), allocatable :: SLT1(:,:,:)
real(kind=8), allocatable :: MRT2(:,:,:)
real(kind=8), allocatable :: MLT2(:,:,:)
real(kind=8), allocatable :: SRT2(:,:,:)
real(kind=8), allocatable :: SLT2(:,:,:)
! data for total system:
! total statistical sum, Boltzmann distribution
real(kind=8), allocatable :: ZT0(:)
real(kind=8), allocatable :: ZT1(:)
real(kind=8), allocatable :: ZT2(:)
real(kind=8), allocatable :: MT0(:,:) ! total magnetisation
real(kind=8), allocatable :: MT1(:,:) ! total magnetisation
real(kind=8), allocatable :: ST0(:,:) ! total spin magnetisation,
real(kind=8), allocatable :: ST1(:,:) ! total spin magnetisation,
real(kind=8), allocatable :: ST2(:,:) ! total spin magnetisation,
real(kind=8), allocatable :: MT2(:,:) ! total magnetisation
! standard deviation data:
real(kind=8) :: dev
real(kind=8), allocatable :: XTM_MH(:)
real(kind=8), allocatable :: XTM_dMdH(:)
real(kind=8), allocatable :: XTtens_MH(:,:,:)
real(kind=8), allocatable :: XTtens_dMdH(:,:,:)
integer :: nTempTotal
real(kind=8) :: Xfield_1, Xfield_2, dltXF, F1, F2
!real(kind=8) :: XTS(nT+nTempMagn)
! Zeeman energy and M vector
!real(kind=8) :: dirX(nDir), dirY(nDir), dirZ(nDir)
!real(kind=8) :: dir_weight(nDirZee,3)
real(kind=8) :: dHX(3)
real(kind=8) :: dHY(3)
real(kind=8) :: dHZ(3)
integer :: RtoB, mem_local
integer :: info, ibuf, ibuf1, ibuf3
real(kind=8) :: wt(3), zt(3,3)
real(kind=8) :: dnrm2_
external :: dev, dnrm2_
logical :: DBG, m_accurate
character(len=50) :: label
real(kind=8), parameter :: cm3tomB = rNAVO*mBohr/Ten ! in cm3 * mol-1 * T

DBG = .false.
!m_paranoid = .true. ! .false.
!ccc-------------------------------------------------------cccc
if (DBG) then
  write(u6,'(A,   I5)') '      exch =',exch
  write(u6,'(A,   I5)') '      nLoc =',nLoc
  write(u6,'(A,   I5)') '   nCenter =',nCenter
  write(u6,'(A,   I5)') '      nneq =',nneq
  write(u6,'(A,   I5)') '      neqv =',neqv
  write(u6,'(A, 10I5)') '    neq(i) =',(neq(i),i=1,nneq)
  write(u6,'(A, 10I5)') '  nexch(i) =',(nexch(i),i=1,nneq)
  write(u6,'(A, 10I5)') '    nss(i) =',(nss(i),i=1,nneq)
  write(u6,'(A,   I5)') ' nTempMagn =',nTempMagn
  write(u6,'(A,   I5)') '        nT =',nT
  write(u6,'(A,   I5)') '        nM =',nM
  write(u6,'(A,   I5)') '      iopt =',iopt
  write(u6,'(A,F12.5)') '      Tmin =',Tmin
  write(u6,'(A,F12.5)') '      Tmax =',Tmax
  write(u6,'(A,F12.5)') '        zJ =',zJ
  write(u6,'(A,F12.5)') '    Xfield =',Xfield
  write(u6,'(A,F12.5)') '        EM =',EM
  write(u6,'(A,ES12.5)') '      THRS =',THRS
  write(u6,*) 'm_accurate =',m_accurate
  write(u6,*) '     smagn =',smagn
  write(u6,*) 'm_paranoid =',m_paranoid
  write(u6,*) '    tinput =',tinput
  if (tinput) then
    do iT=1,nTempMagn
      write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
    end do
    do iT=1,nT
      jT = iT+nTempMagn
      write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',jT,')=',T(jT),' chiT_exp(',iT,')=',chit_exp(iT)
    end do

  else
    do iT=1,nT+nTempMagn
      write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
    end do
  end if
  write(u6,'(A)') 'local ESO:'
  do i=1,nneq
    write(u6,'(A,i3)') 'site:',i
    write(u6,'(10F12.5)') (eso(i,j),j=1,nss(i))
  end do
  write(u6,'(A)') 'EXCHANGE ENERGY'
  write(u6,'(10F12.5)') (W(i),i=1,exch)
  write(u6,'(A)') 'Rotation Matrices:'
  do i=1,nneq
    write(u6,'(A,i3)') 'site:',i
    do j=1,neq(i)
      do l=1,3
        write(u6,'(10F12.5)') (RROT(i,j,l,n),n=1,3)
      end do
    end do
  end do
  call xFlush(u6)
end if !DBG
!ccc-------------------------------------------------------cccc
write(u6,*)
write(u6,'(100A)') (('%'),J=1,95)
write(u6,'(16X,A)') 'CALCULATION OF THE FIELD-DEPENDENT MAGNETIC SUSCEPTIBILITY'
write(u6,'(18X,A)') 'within true (dM/dH) and "experimentalists" (M/H) models'
write(u6,'(100A)') (('%'),J=1,95)
write(u6,*)
write(u6,'(2x,A,F10.6,A)') 'Magnetic field strength:',Xfield,' tesla.'
if (tinput) then
  write(u6,'(2x,a)') 'Temperature dependence of the magnetic susceptibility and'
  write(u6,'(2x,a)') 'high-field magnetization will be calculated according to '
  write(u6,'(2x,a)') 'experimental values provided by the user in file "chitexp.input".'
else
  write(u6,'(2x,a,i3,a)') 'Temperature dependence of the magnetic susceptibility will be calculated in',nT,' points, '
  write(u6,'(2x,a,f4.1,a,f6.1,a)') 'equally distributed in temperature range ',tmin,' ---',tmax,' K.'
end if

!=======================================================================
mem_local = 0
RtoB = 8
call mma_allocate(WEX0,nM,'WEX0')
call mma_allocate(WEX1,nM,'WEX1')
call mma_allocate(WEX2,nM,'WEX2')
call dcopy_(nM,[Zero],0,WEX0,1)
call dcopy_(nM,[Zero],0,WEX1,1)
call dcopy_(nM,[Zero],0,WEX2,1)
mem_local = mem_local+3*nM*RtoB

call mma_allocate(WL0,nneq,nLoc,'WL0')
call mma_allocate(WL1,nneq,nLoc,'WL1')
call mma_allocate(WL2,nneq,nLoc,'WL2')
call mma_allocate(WR0,nneq,nLoc,'WR0')
call mma_allocate(WR1,nneq,nLoc,'WR1')
call mma_allocate(WR2,nneq,nLoc,'WR2')
call dcopy_(nneq*nLoc,[Zero],0,WL0,1)
call dcopy_(nneq*nLoc,[Zero],0,WL1,1)
call dcopy_(nneq*nLoc,[Zero],0,WL2,1)
call dcopy_(nneq*nLoc,[Zero],0,WR0,1)
call dcopy_(nneq*nLoc,[Zero],0,WR1,1)
call dcopy_(nneq*nLoc,[Zero],0,WR2,1)
mem_local = mem_local+6*nneq*nLoc*RtoB

call mma_allocate(ZL0,nneq,(nT+nTempMagn),'ZL0')
call mma_allocate(ZR0,nneq,(nT+nTempMagn),'ZR0')
call mma_allocate(ZL1,nneq,(nT+nTempMagn),'ZL1')
call mma_allocate(ZR1,nneq,(nT+nTempMagn),'ZR1')
call mma_allocate(ZL2,nneq,(nT+nTempMagn),'ZL2')
call mma_allocate(ZR2,nneq,(nT+nTempMagn),'ZR2')
call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZL0,1)
call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZR0,1)
call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZL1,1)
call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZR1,1)
call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZL2,1)
call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZR2,1)
mem_local = mem_local+6*nneq*(nT+nTempMagn)*RtoB

call mma_allocate(SL0,nneq,3,(nT+nTempMagn),'SL0')
call mma_allocate(SR0,nneq,3,(nT+nTempMagn),'SR0')
call mma_allocate(SL1,nneq,3,(nT+nTempMagn),'SL1')
call mma_allocate(SR1,nneq,3,(nT+nTempMagn),'SR1')
call mma_allocate(SL2,nneq,3,(nT+nTempMagn),'SL2')
call mma_allocate(SR2,nneq,3,(nT+nTempMagn),'SR2')
call mma_allocate(ML0,nneq,3,(nT+nTempMagn),'ML0')
call mma_allocate(MR0,nneq,3,(nT+nTempMagn),'MR0')
call mma_allocate(ML1,nneq,3,(nT+nTempMagn),'ML1')
call mma_allocate(MR1,nneq,3,(nT+nTempMagn),'MR1')
call mma_allocate(ML2,nneq,3,(nT+nTempMagn),'ML2')
call mma_allocate(MR2,nneq,3,(nT+nTempMagn),'MR2')
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SL0,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SR0,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SL1,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SR1,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SL2,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SR2,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,ML0,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,MR0,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,ML1,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,MR1,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,ML2,1)
call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,MR2,1)
mem_local = mem_local+12*3*nneq*(nT+nTempMagn)*RtoB

call mma_allocate(ZEX0,(nT+nTempMagn),'ZEX0')
call mma_allocate(ZEX1,(nT+nTempMagn),'ZEX1')
call mma_allocate(ZEX2,(nT+nTempMagn),'ZEX2')
call dcopy_((nT+nTempMagn),[Zero],0,ZEX0,1)
call dcopy_((nT+nTempMagn),[Zero],0,ZEX1,1)
call dcopy_((nT+nTempMagn),[Zero],0,ZEX2,1)
mem_local = mem_local+3*(nT+nTempMagn)*RtoB

call mma_allocate(ZT0,(nT+nTempMagn),'ZT0')
call mma_allocate(ZT1,(nT+nTempMagn),'ZT1')
call mma_allocate(ZT2,(nT+nTempMagn),'ZT2')
call dcopy_((nT+nTempMagn),[Zero],0,ZEX0,1)
call dcopy_((nT+nTempMagn),[Zero],0,ZEX1,1)
call dcopy_((nT+nTempMagn),[Zero],0,ZEX2,1)
mem_local = mem_local+3*(nT+nTempMagn)*RtoB

call mma_allocate(SEX0,3,(nT+nTempMagn),'SEX0')
call mma_allocate(SEX1,3,(nT+nTempMagn),'SEX1')
call mma_allocate(SEX2,3,(nT+nTempMagn),'SEX2')
call mma_allocate(MEX0,3,(nT+nTempMagn),'MEX0')
call mma_allocate(MEX1,3,(nT+nTempMagn),'MEX1')
call mma_allocate(MEX2,3,(nT+nTempMagn),'MEX2')
call dcopy_(3*(nT+nTempMagn),[Zero],0,SEX0,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,SEX1,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,SEX2,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,MEX0,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,MEX1,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,MEX2,1)
mem_local = mem_local+6*3*(nT+nTempMagn)*RtoB

call mma_allocate(ST0,3,(nT+nTempMagn),'ST0')
call mma_allocate(ST1,3,(nT+nTempMagn),'ST1')
call mma_allocate(ST2,3,(nT+nTempMagn),'ST2')
call mma_allocate(MT0,3,(nT+nTempMagn),'MT0')
call mma_allocate(MT1,3,(nT+nTempMagn),'MT1')
call mma_allocate(MT2,3,(nT+nTempMagn),'MT2')
call dcopy_(3*(nT+nTempMagn),[Zero],0,ST0,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,ST1,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,ST2,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,MT0,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,MT1,1)
call dcopy_(3*(nT+nTempMagn),[Zero],0,MT2,1)
mem_local = mem_local+6*3*(nT+nTempMagn)*RtoB

call mma_allocate(XTM_MH,(nT+nTempMagn),'XTM_MH')
call mma_allocate(XTM_dMdH,(nT+nTempMagn),'XTM_dMdH')
call mma_allocate(XTtens_MH,3,3,(nT+nTempMagn),'XTtens_MH')
call mma_allocate(XTtens_dMdH,3,3,(nT+nTempMagn),'XTtens_dMdH')
call dcopy_((nT+nTempMagn),[Zero],0,XTM_MH,1)
call dcopy_((nT+nTempMagn),[Zero],0,XTM_dMdH,1)
call dcopy_(3*3*(nT+nTempMagn),[Zero],0,XTtens_MH,1)
call dcopy_(3*3*(nT+nTempMagn),[Zero],0,XTtens_dMdH,1)
mem_local = mem_local+20*(nT+nTempMagn)*RtoB

call mma_allocate(ZRT0,nCenter,(nT+nTempMagn),'ZRT0')
call mma_allocate(ZLT0,nCenter,(nT+nTempMagn),'ZLT0')
call mma_allocate(ZRT1,nCenter,(nT+nTempMagn),'ZRT1')
call mma_allocate(ZLT1,nCenter,(nT+nTempMagn),'ZLT1')
call mma_allocate(ZRT2,nCenter,(nT+nTempMagn),'ZRT2')
call mma_allocate(ZLT2,nCenter,(nT+nTempMagn),'ZLT2')
call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZRT0,1)
call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZLT0,1)
call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZRT1,1)
call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZLT1,1)
call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZRT2,1)
call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZLT2,1)
mem_local = mem_local+6*3*nCenter*(nT+nTempMagn)*RtoB

call mma_allocate(MRT0,nCenter,3,(nT+nTempMagn),'MRT0')
call mma_allocate(MLT0,nCenter,3,(nT+nTempMagn),'MLT0')
call mma_allocate(SRT0,nCenter,3,(nT+nTempMagn),'SRT0')
call mma_allocate(SLT0,nCenter,3,(nT+nTempMagn),'SLT0')
call mma_allocate(MRT1,nCenter,3,(nT+nTempMagn),'MRT1')
call mma_allocate(MLT1,nCenter,3,(nT+nTempMagn),'MLT1')
call mma_allocate(SRT1,nCenter,3,(nT+nTempMagn),'SRT1')
call mma_allocate(SLT1,nCenter,3,(nT+nTempMagn),'SLT1')
call mma_allocate(MRT2,nCenter,3,(nT+nTempMagn),'MRT2')
call mma_allocate(MLT2,nCenter,3,(nT+nTempMagn),'MLT2')
call mma_allocate(SRT2,nCenter,3,(nT+nTempMagn),'SRT2')
call mma_allocate(SLT2,nCenter,3,(nT+nTempMagn),'SLT2')
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MRT0,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MLT0,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SRT0,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SLT0,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MRT1,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MLT1,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SRT1,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SLT1,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MRT2,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MLT2,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SRT2,1)
call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SLT2,1)
mem_local = mem_local+12*3*nCenter*(nT+nTempMagn)*RtoB

if (dbg) write(u6,*) 'XT_MH:  memory allocated (local):'
if (dbg) write(u6,*) 'mem_local=',mem_local
if (dbg) write(u6,*) 'XT_MH:  memory allocated (total):'
if (dbg) write(u6,*) 'mem_total=',mem+mem_local

!=======================================================================
nDirX = 3

dHX(1) = One
dHY(1) = Zero
dHZ(1) = Zero

dHX(2) = Zero
dHY(2) = One
dHZ(2) = Zero

dHX(3) = Zero
dHY(3) = Zero
dHZ(3) = One

Xfield_1 = Xfield*0.99999_wp
Xfield_2 = Xfield*1.00001_wp
dltXF = Xfield_2-Xfield_1

nTempTotal = nT+nTempMagn

! ///  opening the loop over different directions of the magnetic field
do iM=1,nDirX
  call dcopy_(nM,[Zero],0,WEX0,1)
  call dcopy_(nM,[Zero],0,WEX1,1)
  call dcopy_(nM,[Zero],0,WEX2,1)
  call dcopy_((nT+nTempMagn),[Zero],0,ZEX0,1)
  call dcopy_((nT+nTempMagn),[Zero],0,ZEX1,1)
  call dcopy_((nT+nTempMagn),[Zero],0,ZEX2,1)
  call dcopy_(3*(nT+nTempMagn),[Zero],0,SEX0,1)
  call dcopy_(3*(nT+nTempMagn),[Zero],0,SEX1,1)
  call dcopy_(3*(nT+nTempMagn),[Zero],0,SEX2,1)
  call dcopy_(3*(nT+nTempMagn),[Zero],0,MEX0,1)
  call dcopy_(3*(nT+nTempMagn),[Zero],0,MEX1,1)
  call dcopy_(3*(nT+nTempMagn),[Zero],0,MEX2,1)
  call dcopy_(nneq*nLoc,[Zero],0,WL0,1)
  call dcopy_(nneq*nLoc,[Zero],0,WL1,1)
  call dcopy_(nneq*nLoc,[Zero],0,WL2,1)
  call dcopy_(nneq*nLoc,[Zero],0,WR0,1)
  call dcopy_(nneq*nLoc,[Zero],0,WR1,1)
  call dcopy_(nneq*nLoc,[Zero],0,WR2,1)
  call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZL0,1)
  call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZR0,1)
  call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZL1,1)
  call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZR1,1)
  call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZL2,1)
  call dcopy_(nneq*(nT+nTempMagn),[Zero],0,ZR2,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SL0,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SR0,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SL1,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SR1,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SL2,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,SR2,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,ML0,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,MR0,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,ML1,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,MR1,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,ML2,1)
  call dcopy_(3*nneq*(nT+nTempMagn),[Zero],0,MR2,1)

  call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZRT0,1)
  call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZLT0,1)
  call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZRT1,1)
  call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZLT1,1)
  call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZRT2,1)
  call dcopy_(nCenter*(nT+nTempMagn),[Zero],0,ZLT2,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MRT0,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MLT0,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SRT0,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SLT0,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MRT1,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MLT1,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SRT1,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SLT1,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MRT2,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,MLT2,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SRT2,1)
  call dcopy_(3*nCenter*(nT+nTempMagn),[Zero],0,SLT2,1)
  ! exchange magnetization:
  call MAGN(EXCH,NM,dHX(iM),dHY(iM),dHZ(iM),XField,W,zJ,THRS,DIPEXCH(1:3,1:exch,1:exch),S_EXCH(1:3,1:exch,1:exch),nTempTotal, &
            T(1:nTempTotal),smagn,Wex0(1:nM),Zex0(1:nTempTotal),Sex0(1:3,1:nTempTotal),Mex0(1:3,1:nTempTotal),m_paranoid,DBG)
  call MAGN(EXCH,NM,dHX(iM),dHY(iM),dHZ(iM),XField_1,W,zJ,THRS,DIPEXCH(1:3,1:exch,1:exch),S_EXCH(1:3,1:exch,1:exch),nTempTotal, &
            T(1:nTempTotal),smagn,Wex1(1:nM),Zex1(1:nTempTotal),Sex1(1:3,1:nTempTotal),Mex1(1:3,1:nTempTotal),m_paranoid,DBG)
  call MAGN(EXCH,NM,dHX(iM),dHY(iM),dHZ(iM),XField_2,W,zJ,THRS,DIPEXCH(1:3,1:exch,1:exch),S_EXCH(1:3,1:exch,1:exch),nTempTotal, &
            T(1:nTempTotal),smagn,Wex2(1:nM),Zex2(1:nTempTotal),Sex2(1:3,1:nTempTotal),Mex2(1:3,1:nTempTotal),m_paranoid,DBG)
  if (DBG) then
    do i=1,nTempTotal
      write(u6,'(A,i3,A,9ES22.14)') 'Mex0(',i,')=',(Mex0(l,i),l=1,3)
    end do
    do i=1,nTempTotal
      write(u6,'(A,i3,A,9ES22.14)') 'Mex1(',i,')=',(Mex1(l,i),l=1,3)
    end do
    do i=1,nTempTotal
      write(u6,'(A,i3,A,9ES22.14)') 'Mex2(',i,')=',(Mex2(l,i),l=1,3)
    end do
    do i=1,nTempTotal
      write(u6,'(A,i3,A,9ES22.14)') 'Mex 2-1 (',i,')=',(Mex2(l,i)-Mex1(l,i),l=1,3)
    end do
  end if !DBG

  call Add_Info('dM/dH   Mex0',[dnrm2_(3*nTempTotal,Mex0,1)],1,8)
  call Add_Info('dM/dH   Sex0',[dnrm2_(3*nTempTotal,Sex0,1)],1,8)
  call Add_Info('dM/dH   Zex0',[dnrm2_(nTempTotal,Zex0,1)],1,8)
  call Add_Info('dM/dH   Wex0',[dnrm2_(nM,Wex0,1)],1,8)

  call Add_Info('dM/dH   Mex1',[dnrm2_(3*nTempTotal,Mex1,1)],1,8)
  call Add_Info('dM/dH   Sex1',[dnrm2_(3*nTempTotal,Sex1,1)],1,8)
  call Add_Info('dM/dH   Zex1',[dnrm2_(nTempTotal,Zex1,1)],1,8)
  call Add_Info('dM/dH   Wex1',[dnrm2_(nM,Wex1,1)],1,8)

  call Add_Info('dM/dH   Mex2',[dnrm2_(3*nTempTotal,Mex2,1)],1,8)
  call Add_Info('dM/dH   Sex2',[dnrm2_(3*nTempTotal,Sex2,1)],1,8)
  call Add_Info('dM/dH   Zex2',[dnrm2_(nTempTotal,Zex2,1)],1,8)
  call Add_Info('dM/dH   Wex2',[dnrm2_(nM,Wex2,1)],1,8)

  ! compute local magnetizations:
  if (m_accurate) then
    do i=1,nneq
      ! all states:
      if (NSS(i) > NEXCH(i)) then
        ! this check is to avoid the unnecessary computation, in cases when no local excited states are present
        call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField,ESO(i,1:NSS(i)),zJ,THRS,DIPSO(i,1:3,1:NSS(i),1:NSS(i)), &
                  S_SO(i,1:3,1:NSS(i),1:NSS(i)),nTempTotal,T(1:(nT+nTempMagn)),smagn,WL0(i,1:NSS(i)),ZL0(i,1:(nT+nTempMagn)), &
                  SL0(i,1:3,1:(nT+nTempMagn)),ML0(i,1:3,1:(nT+nTempMagn)),m_paranoid,DBG)

        call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField_1,ESO(i,1:NSS(i)),zJ,THRS,DIPSO(i,1:3,1:NSS(i),1:NSS(i)), &
                  S_SO(i,1:3,1:NSS(i),1:NSS(i)),nTempTotal,T(1:(nT+nTempMagn)),smagn,WL1(i,1:NSS(i)),ZL1(i,1:(nT+nTempMagn)), &
                  SL1(i,1:3,1:(nT+nTempMagn)),ML1(i,1:3,1:(nT+nTempMagn)),m_paranoid,DBG)

        call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),Xfield_2,ESO(i,1:NSS(i)),zJ,THRS,DIPSO(i,1:3,1:NSS(i),1:NSS(i)), &
                  S_SO(i,1:3,1:NSS(i),1:NSS(i)),nTempTotal,T(1:(nT+nTempMagn)),smagn,WL2(i,1:NSS(i)),ZL2(i,1:(nT+nTempMagn)), &
                  SL2(i,1:3,1:(nT+nTempMagn)),ML2(i,1:3,1:(nT+nTempMagn)),m_paranoid,DBG)
! only local "exchange states":
        call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField,ESO(i,1:NEXCH(i)),zJ,THRS,DIPSO(i,1:3,1:NEXCH(i),1:NEXCH(i)), &
                  S_SO(i,1:3,1:NEXCH(i),1:NEXCH(i)),nTempTotal,T(1:(nT+nTempMagn)),smagn,WR0(i,1:Nexch(i)), &
                  ZR0(i,1:(nT+nTempMagn)),SR0(i,1:3,1:(nT+nTempMagn)),MR0(i,1:3,1:(nT+nTempMagn)),m_paranoid,DBG)

        call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField_1,ESO(i,1:NEXCH(i)),zJ,THRS,DIPSO(i,1:3,1:NEXCH(i),1:NEXCH(i)), &
                  S_SO(i,1:3,1:NEXCH(i),1:NEXCH(i)),nTempTotal,T(1:(nT+nTempMagn)),smagn,WR1(i,1:Nexch(i)), &
                  ZR1(i,1:(nT+nTempMagn)),SR1(i,1:3,1:(nT+nTempMagn)),MR1(i,1:3,1:(nT+nTempMagn)),m_paranoid,DBG)

        call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField_2,ESO(i,1:NEXCH(i)),zJ,THRS,DIPSO(i,1:3,1:NEXCH(i),1:NEXCH(i)), &
                  S_SO(i,1:3,1:NEXCH(i),1:NEXCH(i)),nTempTotal,T(1:(nT+nTempMagn)),smagn,WR2(i,1:Nexch(i)), &
                  ZR2(i,1:(nT+nTempMagn)),SR2(i,1:3,1:(nT+nTempMagn)),MR2(i,1:3,1:(nT+nTempMagn)),m_paranoid,DBG)
        if (DBG) then
          do iT=1,nTempTotal
            write(u6,'(A,i1,A,i3,A,3ES22.14)') 'ML(',i,',L,',iT,')=',(ML2(i,l,iT)-ML1(i,l,iT),l=1,3)
            write(u6,'(A,i1,A,i3,A,3ES22.14)') 'MR(',i,',L,',iT,')=',(MR2(i,l,iT)-MR1(i,l,iT),l=1,3)
          end do
        end if ! DBG
      end if ! NSS(i) > NEXCH(i)
    end do ! i=1,nneq

    ibuf1 = nTempTotal*nneq
    ibuf3 = 3*ibuf1
    call Add_Info('dM/dH    ML0',[dnrm2_(ibuf3,ML0,1)],1,8)
    call Add_Info('dM/dH    SL0',[dnrm2_(ibuf3,SL0,1)],1,8)
    call Add_Info('dM/dH    ZL0',[dnrm2_(ibuf1,ZL0,1)],1,8)
    call Add_Info('dM/dH    WL0',[dnrm2_(nLoc*nneq,WL0,1)],1,8)

    call Add_Info('dM/dH    MR0',[dnrm2_(ibuf3,MR0,1)],1,8)
    call Add_Info('dM/dH    SR0',[dnrm2_(ibuf3,SR0,1)],1,8)
    call Add_Info('dM/dH    ZR0',[dnrm2_(ibuf1,ZR0,1)],1,8)
    call Add_Info('dM/dH    WR0',[dnrm2_(nLoc*nneq,WR0,1)],1,8)

    call Add_Info('dM/dH    ML1',[dnrm2_(ibuf3,ML1,1)],1,8)
    call Add_Info('dM/dH    SL1',[dnrm2_(ibuf3,SL1,1)],1,8)
    call Add_Info('dM/dH    ZL1',[dnrm2_(ibuf1,ZL1,1)],1,8)
    call Add_Info('dM/dH    WL1',[dnrm2_(nLoc*nneq,WL1,1)],1,8)

    call Add_Info('dM/dH    MR1',[dnrm2_(ibuf3,MR1,1)],1,8)
    call Add_Info('dM/dH    SR1',[dnrm2_(ibuf3,SR1,1)],1,8)
    call Add_Info('dM/dH    ZR1',[dnrm2_(ibuf1,ZR1,1)],1,8)
    call Add_Info('dM/dH    WR1',[dnrm2_(nLoc*nneq,WR1,1)],1,8)

    call Add_Info('dM/dH    ML2',[dnrm2_(ibuf3,ML2,1)],1,8)
    call Add_Info('dM/dH    SL2',[dnrm2_(ibuf3,SL2,1)],1,8)
    call Add_Info('dM/dH    ZL2',[dnrm2_(ibuf1,ZL2,1)],1,8)
    call Add_Info('dM/dH    WL2',[dnrm2_(nLoc*nneq,WL2,1)],1,8)

    call Add_Info('dM/dH    MR2',[dnrm2_(ibuf3,MR2,1)],1,8)
    call Add_Info('dM/dH    SR2',[dnrm2_(ibuf3,SR2,1)],1,8)
    call Add_Info('dM/dH    ZR2',[dnrm2_(ibuf1,ZR2,1)],1,8)
    call Add_Info('dM/dH    WR2',[dnrm2_(nLoc*nneq,WR2,1)],1,8)

    ! expand the basis and rotate local vectors to the general
    ! coordinate system:

    isite = 0
    do i=1,NNEQ
      do j=1,NEQ(i)
        isite = isite+1
        ! statistical distributions
        do iT=1,nTempTotal
          ZLT0(isite,iT) = ZL0(i,iT)
          ZRT0(isite,iT) = ZR0(i,iT)
          ZLT1(isite,iT) = ZL1(i,iT)
          ZRT1(isite,iT) = ZR1(i,iT)
          ZLT2(isite,iT) = ZL2(i,iT)
          ZRT2(isite,iT) = ZR2(i,iT)
        end do
        ! magnetizations:
        !    use R_rot matrices, which have determinant +1.
        !  note that  R_lg matrices have arbitrary determinant.
        do iT=1,nTempTotal
          do l=1,3
            do n=1,3
              MLT0(isite,l,iT) = MLT0(isite,l,iT)+rrot(i,j,l,n)*ML0(i,n,iT)
              SLT0(isite,l,iT) = SLT0(isite,l,iT)+rrot(i,j,l,n)*SL0(i,n,iT)
              MRT0(isite,l,iT) = MRT0(isite,l,iT)+rrot(i,j,l,n)*MR0(i,n,iT)
              SRT0(isite,l,iT) = SRT0(isite,l,iT)+rrot(i,j,l,n)*SR0(i,n,iT)

              MLT1(isite,l,iT) = MLT1(isite,l,iT)+rrot(i,j,l,n)*ML1(i,n,iT)
              SLT1(isite,l,iT) = SLT1(isite,l,iT)+rrot(i,j,l,n)*SL1(i,n,iT)
              MRT1(isite,l,iT) = MRT1(isite,l,iT)+rrot(i,j,l,n)*MR1(i,n,iT)
              SRT1(isite,l,iT) = SRT1(isite,l,iT)+rrot(i,j,l,n)*SR1(i,n,iT)

              MLT2(isite,l,iT) = MLT2(isite,l,iT)+rrot(i,j,l,n)*ML2(i,n,iT)
              SLT2(isite,l,iT) = SLT2(isite,l,iT)+rrot(i,j,l,n)*SL2(i,n,iT)
              MRT2(isite,l,iT) = MRT2(isite,l,iT)+rrot(i,j,l,n)*MR2(i,n,iT)
              SRT2(isite,l,iT) = SRT2(isite,l,iT)+rrot(i,j,l,n)*SR2(i,n,iT)
            end do
          end do
        end do

      end do ! j, neq(i)
    end do ! i, nneq
    ibuf = 0
    ibuf = 3*(nT+nTempMagn)*nCenter
    call Add_Info('dM/dH    MLT0',[dnrm2_(ibuf,MLT0,1)],1,8)
    call Add_Info('dM/dH    SLT0',[dnrm2_(ibuf,SLT0,1)],1,8)
    call Add_Info('dM/dH    MRT0',[dnrm2_(ibuf,MRT0,1)],1,8)
    call Add_Info('dM/dH    SRT0',[dnrm2_(ibuf,SRT0,1)],1,8)

    call Add_Info('dM/dH    MLT1',[dnrm2_(ibuf,MLT1,1)],1,8)
    call Add_Info('dM/dH    SLT1',[dnrm2_(ibuf,SLT1,1)],1,8)
    call Add_Info('dM/dH    MRT1',[dnrm2_(ibuf,MRT1,1)],1,8)
    call Add_Info('dM/dH    SRT1',[dnrm2_(ibuf,SRT1,1)],1,8)

    call Add_Info('dM/dH    MLT2',[dnrm2_(ibuf,MLT2,1)],1,8)
    call Add_Info('dM/dH    SLT2',[dnrm2_(ibuf,SLT2,1)],1,8)
    call Add_Info('dM/dH    MRT2',[dnrm2_(ibuf,MRT2,1)],1,8)
    call Add_Info('dM/dH    SRT2',[dnrm2_(ibuf,SRT2,1)],1,8)
    ! compute the total magnetizations according to the derived formulas:
    if (dbg) write(u6,*) 'check point MLT0'
    do iT=1,nTempTotal
      if (smagn) then
        call MSUM(nCenter,Sex0(1:3,iT),Zex0(iT),SLT0(1:nCenter,1:3,iT),ZLT0(1:nCenter,iT),SRT0(1:nCenter,1:3,iT), &
                  ZRT0(1:nCenter,iT),iopt,ST0(1:3,iT),ZT0(iT))
        call MSUM(nCenter,Sex1(1:3,iT),Zex1(iT),SLT1(1:nCenter,1:3,iT),ZLT1(1:nCenter,iT),SRT1(1:nCenter,1:3,iT), &
                  ZRT1(1:nCenter,iT),iopt,ST1(1:3,iT),ZT1(iT))
        call MSUM(nCenter,Sex2(1:3,iT),Zex2(iT),SLT2(1:nCenter,1:3,iT),ZLT2(1:nCenter,iT),SRT2(1:nCenter,1:3,iT), &
                  ZRT2(1:nCenter,iT),iopt,ST2(1:3,iT),ZT2(iT))
      end if

      if (dbg) write(u6,*) 'check point MSUM'
      call MSUM(nCenter,Mex0(1:3,iT),Zex0(iT),MLT0(1:nCenter,1:3,iT),ZLT0(1:nCenter,iT),MRT0(1:nCenter,1:3,iT),ZRT0(1:nCenter,iT), &
                iopt,MT0(1:3,iT),ZT0(iT))
      call MSUM(nCenter,Mex1(1:3,iT),Zex1(iT),MLT1(1:nCenter,1:3,iT),ZLT1(1:nCenter,iT),MRT1(1:nCenter,1:3,iT),ZRT1(1:nCenter,iT), &
                iopt,MT1(1:3,iT),ZT1(iT))
      call MSUM(nCenter,Mex2(1:3,iT),Zex2(iT),MLT2(1:nCenter,1:3,iT),ZLT2(1:nCenter,iT),MRT2(1:nCenter,1:3,iT),ZRT2(1:nCenter,iT), &
                iopt,MT2(1:3,iT),ZT2(iT))
    end do

  end if ! m_accurate

  ibuf = 0
  ibuf = 3*(nT+nTempMagn)
  call Add_Info('dM/dH    ST0',[dnrm2_(ibuf,ST0,1)],1,6)
  call Add_Info('dM/dH    ST1',[dnrm2_(ibuf,ST1,1)],1,6)
  call Add_Info('dM/dH    ST2',[dnrm2_(ibuf,ST2,1)],1,6)
  call Add_Info('dM/dH    MT0',[dnrm2_(ibuf,MT0,1)],1,6)
  call Add_Info('dM/dH    MT1',[dnrm2_(ibuf,MT1,1)],1,6)
  call Add_Info('dM/dH    MT2',[dnrm2_(ibuf,MT2,1)],1,6)
  ! computing the AVERAGE MOMENTS calculated at different temperatures (T(i))
  do iT=1,nTempTotal
    F1 = T(iT)*cm3tomB/dltXF
    F2 = T(iT)*cm3tomB/Xfield

    ! dM/dH model
    XTtens_dMdH(iM,1,iT) = (MT2(1,iT)-MT1(1,iT))*F1
    XTtens_dMdH(iM,2,iT) = (MT2(2,iT)-MT1(2,iT))*F1
    XTtens_dMdH(iM,3,iT) = (MT2(3,iT)-MT1(3,iT))*F1

    ! M/H model
    XTtens_MH(iM,1,iT) = MT0(1,iT)*F2
    XTtens_MH(iM,2,iT) = MT0(2,iT)*F2
    XTtens_MH(iM,3,iT) = MT0(3,iT)*F2
  end do
  ! ///  closing the loops over field directions
end do ! iM (nDirX)
! computing the XT as tensor's average:
do iT=1,nTempTotal
  XTM_dMdH(iT) = (XTtens_dMdH(1,1,iT)+XTtens_dMdH(2,2,iT)+XTtens_dMdH(3,3,iT))/Three

  XTM_MH(iT) = (XTtens_MH(1,1,iT)+XTtens_MH(2,2,iT)+XTtens_MH(3,3,iT))/Three
end do !iT

ibuf = 0
ibuf = 9*(nT+nTempMagn)
call Add_Info('dM/dH    XTtens_dMdH',[dnrm2_(ibuf,XTtens_dMdH,1)],1,5)
call Add_Info('dM/dH    XTtens_MH',[dnrm2_(ibuf,XTtens_MH,1)],1,6)
ibuf = 0
ibuf = nT+nTempMagn
call Add_Info('dM/dH    XTM_dMdH',[dnrm2_(ibuf,XTM_dMdH,1)],1,5)
call Add_Info('dM/dH    XTM_MH',[dnrm2_(ibuf,XTM_MH,1)],1,6)
! ----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
! ----------------------------------------------------------------------
write(u6,*)
write(u6,'(A)') '--------------------------------------------------------------------------|'
write(u6,'(A)') '     |     T      | Statistical |   CHI*T     |   CHI*T     |   CHI*T     |'
write(u6,'(A,F6.3,A,F6.3,A)') '     |            |  Sum (Z)    | H = 0.000 T | H =',Xfield,' T | H =',XField,' T |'
write(u6,'(A)') '     |            |             |             | X = dM/dH   | X = M/H     |'
write(u6,'(A)') '-----|--------------------------------------------------------------------|'
write(u6,'(A)') 'Units|   kelvin   |    ---      |  cm3*K/mol  |  cm3*K/mol  |  cm3*K/mol  |'
write(u6,'(A)') '-----|--------------------------------------------------------------------|'

do iT=1,nT
  jT = iT+nTempMagn
  write(u6,'(A,F11.6,A,ES12.5,A,F12.8,A,F12.8,A,F12.7,A)') '     |',T(jT),' |',ZT0(jT),' |',XT_no_field(jT),' |',XTM_dMdH(jT), &
                                                           ' |',XTM_MH(jT),' |'
end do
write(u6,'(A)') '-----|--------------------------------------------------------------------|'

!  calcualtion of the standard deviation:
if (tinput) then
  write(u6,'(a,5x, f20.14)') 'ST.DEV: X= dM/dH:',dev(nT,XTM_dMdH((1+nTempMagn):(nT+nTempMagn)), &
                             chit_exp((1+nTempMagn):(nT+nTempMagn)))
  write(u6,'(a,5x, f20.14)') 'ST.DEV: X=  M/H :',dev(nT,XTM_MH((1+nTempMagn):(nT+nTempMagn)),chit_exp((1+nTempMagn):(nT+nTempMagn)))
  write(u6,'(A)') '-----|--------------------------------------------------------------------|'
end if !tinput

!-------------------------  PLOTs -------------------------------------!
write(label,'(A)') 'with_field_M_over_H'
call xFlush(u6)
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T((1+nTempMagn):(nT+nTempMagn)),XTM_MH((1+nTempMagn):(nT+nTempMagn)), &
                          chit_exp((1+nTempMagn):(nT+nTempMagn)))
  else
    call plot_XT_no_Exp(label,nT,T((1+nTempMagn):(nT+nTempMagn)),XTM_MH((1+nTempMagn):(nT+nTempMagn)))
  end if
end if

write(label,'(A)') 'with_field_dM_over_dH'
call xFlush(u6)
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T((1+nTempMagn):(nT+nTempMagn)),XTM_dMdH((1+nTempMagn):(nT+nTempMagn)), &
                          chit_exp((1+nTempMagn):(nT+nTempMagn)))
  else
    call plot_XT_no_Exp(label,nT,T((1+nTempMagn):(nT+nTempMagn)),XTM_dMdH((1+nTempMagn):(nT+nTempMagn)))
  end if
end if
!------------------------- END PLOTs ----------------------------------!

! print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
write(u6,'(/)')
write(u6,'(111A)') ('-',i=1,110),'|'
write(u6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR the  X=dM/dH  model,  in cm3*K/mol','|'
write(u6,'(111A)') ('-',i=1,110),'|'
write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
do iT=1,nT
  jT = iT+nTempMagn
  info = 0
  call DIAG_R2(XTtens_dMdH(1:3,1:3,jT),3,info,wt,zt)
  write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- z ----|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(XTtens_dMdH(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                  (zt(j,1),j=1,3),'|'
  write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(XTtens_dMdH(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                           (zt(j,2),j=1,3),'|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(XTtens_dMdH(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                  (zt(j,3),j=1,3),'|'
end do
write(u6,'(111A)') ('-',i=1,110),'|'

write(u6,'(/)')
write(u6,'(111A)') ('-',i=1,110),'|'

write(u6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR the  X = M/H  model,  in cm3*K/mol','|'
write(u6,'(111A)') ('-',i=1,110),'|'
write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
do iT=1,nT
  jT = iT+nTempMagn
  info = 0
  call DIAG_R2(XTtens_MH(1:3,1:3,jT),3,info,wt,zt)
  write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- z ----|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(XTtens_MH(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                  (zt(j,1),j=1,3),'|'
  write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(XTtens_MH(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                           (zt(j,2),j=1,3),'|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(XTtens_MH(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                  (zt(j,3),j=1,3),'|'
end do
write(u6,'(111A)') ('-',i=1,110),'|'

!=======================================================================
call mma_deallocate(WEX0)
call mma_deallocate(WEX1)
call mma_deallocate(WEX2)
call mma_deallocate(WL0)
call mma_deallocate(WL1)
call mma_deallocate(WL2)
call mma_deallocate(WR0)
call mma_deallocate(WR1)
call mma_deallocate(WR2)

call mma_deallocate(ZL0)
call mma_deallocate(ZR0)
call mma_deallocate(ZL1)
call mma_deallocate(ZR1)
call mma_deallocate(ZL2)
call mma_deallocate(ZR2)
call mma_deallocate(SL0)
call mma_deallocate(SR0)
call mma_deallocate(SL1)
call mma_deallocate(SR1)
call mma_deallocate(SL2)
call mma_deallocate(SR2)
call mma_deallocate(ML0)
call mma_deallocate(MR0)
call mma_deallocate(ML1)
call mma_deallocate(MR1)
call mma_deallocate(ML2)
call mma_deallocate(MR2)

call mma_deallocate(ZEX0)
call mma_deallocate(ZEX1)
call mma_deallocate(ZEX2)
call mma_deallocate(ZT0)
call mma_deallocate(ZT1)
call mma_deallocate(ZT2)
call mma_deallocate(SEX0)
call mma_deallocate(SEX1)
call mma_deallocate(SEX2)
call mma_deallocate(MEX0)
call mma_deallocate(MEX1)
call mma_deallocate(MEX2)
call mma_deallocate(ST0)
call mma_deallocate(ST1)
call mma_deallocate(ST2)
call mma_deallocate(MT0)
call mma_deallocate(MT1)
call mma_deallocate(MT2)
call mma_deallocate(XTM_MH)
call mma_deallocate(XTM_dMdH)
call mma_deallocate(XTtens_MH)
call mma_deallocate(XTtens_dMdH)

call mma_deallocate(ZRT0)
call mma_deallocate(ZLT0)
call mma_deallocate(ZRT1)
call mma_deallocate(ZLT1)
call mma_deallocate(ZRT2)
call mma_deallocate(ZLT2)
call mma_deallocate(MRT0)
call mma_deallocate(MLT0)
call mma_deallocate(SRT0)
call mma_deallocate(SLT0)
call mma_deallocate(MRT1)
call mma_deallocate(MLT1)
call mma_deallocate(SRT1)
call mma_deallocate(SLT1)
call mma_deallocate(MRT2)
call mma_deallocate(MLT2)
call mma_deallocate(SRT2)
call mma_deallocate(SLT2)

return

end subroutine XT_dMoverdH
