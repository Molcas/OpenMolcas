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
! Copyright (C) 2015, Kamal Sharkas                                    *
!               2019, Thomas J. Duignan                                *
!               2021, Rulin Feng                                       *
!**********************************************************************

! Note: The hyperfine code is based on the analogous
! pre-existing G-tensor functionality
subroutine HFCTS(PROP,USOR,USOI,ENSOR,NSS,ENERGY,JBNUM,DIPSOM,ESO,XYZCHR,BOLTZ_K)

use rassi_aux, only: ipglob
use Cntrl, only: EPRATHR, ICOMP, IFACALFC, IFACALFCON, IFACALFCSDON, IFACALPSO, IFACALSD, IFACALSDON, IFATCALSA, IFGTSHSA, &
                 IFSONCINI, MLTPLT, MULTIP, NPROP, NSTATE, NTP, PNAME, TMAXP, TMINP
use hfc_logical, only: MAG_X2C
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart, OneHalf, cZero, auTocm, c_in_au, gElectron
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSS, JBNUM(NSTATE)
real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP), USOR(NSS,NSS), USOI(NSS,NSS), ENSOR(NSS), ENERGY(NSTATE), ESO(NSS), Boltz_k
complex(kind=wp), intent(in) :: DIPSOm(3,NSS,NSS)
character, intent(in) :: xyzchr(3)
integer(kind=iwp) :: I, IAMFI1, IAMFI2, IAMFI3, IAMFI4, IAMFI5, IAMFI6, IC, ICEN, ICOUNT, IERR, iFinal, IFUNCT, IMLTPL, IPROP, &
                     ISO, ISS, IStart, ISTATE, IT, IXYZ, J, JC, JOB, JSS, JSTATE, JXYZ, KDGN, KPROP, KXYZ, MPLET, MPLET1, MPLET2, &
                     MSPROJ, MSPROJ1, MSPROJ2
real(kind=wp) :: ACNT, Alpha, Alpha2, AMFI1, AMFI2, AMFI3, AMFI4, AMFI5, AMFI6, CG0, CGM, CGP, CGX, CGY, CurieT(3,3), DCLEBS, &
                 DiamT(3,3), DLT_E, DLTTA, EDIFF, EEX, EEY, EEZ, EVI(3), EVR(3), FACT, FEGVAL, GSENERGY, GTENS(3,3), GTIJ, &
                 GTOTAL(3,3), HFC_1(3,3), HFC_2(3,3), HFC_3(3,3), p_Boltz, S1, S2, SM1, SM2, TMPMAT(3,3), TMPVEC(3,3), Zstat
character(len=8) :: DMPPROP, PSOPROP, SDPROP
integer(kind=iwp), allocatable :: MAPST(:), MAPSP(:), MAPMS(:)
real(kind=wp), allocatable :: LXI(:,:), LYI(:,:), LZI(:,:), MXI(:,:), MXR(:,:), MYI(:,:), MYR(:,:), MZI(:,:), MZR(:,:), &
                              PNMR(:,:,:), PNMRC(:,:,:), PNMRCPS(:,:,:,:), PNMRD(:,:,:), PNMRT(:,:,:), SOPRI(:,:), SOPRR(:,:), &
                              TMPf(:), ZXYZI(:,:,:), ZXYZR(:,:,:)
complex(kind=wp), allocatable :: DIMSO(:,:,:,:), DIPSOf(:,:,:), DIPSOfc(:,:,:), DIPSOfcsd(:,:,:), DIPSOfpso(:,:,:), &
                                 DIPSOfsd(:,:,:), GCONT(:,:,:), SPNSFS(:,:,:), Z(:,:), ZEKL(:,:,:,:)
logical(kind=iwp), allocatable :: ISGS(:)
real(kind=wp), parameter :: THRSH = 1.0e-10_wp

!AU2J = auTokJ*1.0e3_wp
!J2CM = auTocm/AU2J
!AU2JTM = (AU2J/auToT)*rNAVO
ALPHA = One/c_in_au
ALPHA2 = ALPHA*ALPHA
!AU2REDR = 200.0_wp*Debye

!coeff_chi = 0.1_wp*rNAVO/kBoltzmann*mBohr**2
FEGVAL = -gElectron
!BOLTZ = kBoltzmann/AU2J
!Rmu0 = 4.0e-7_wp*Pi

if (IFSONCINI) then
  write(u6,*)
  write(u6,*) '  Soncini pNMR Tensor and A-Matrix Approach II'
  write(u6,*) '  ============================================='
  write(u6,*) '  1st order degenerate perturbation theory'
  write(u6,*) '  within isolated kramers doublets.'
  write(u6,*) '  > spatial degeneracy'
  write(u6,*) '  > strong spin-orbit coupling'
  write(u6,*)
else
  write(u6,*)
  write(u6,*) '  A-Matrix Approach II'
  write(u6,*) '  ========================================='
  write(u6,*) '  1st order degenerate perturbation theory'
  write(u6,*) '  within isolated kramers doublets.'
  write(u6,*) '  > spatial degeneracy'
  write(u6,*) '  > strong spin-orbit coupling'
  write(u6,*)
end if

! Mapping from spin states to spin-free state and to spin:
call mma_allocate(MAPST,NSS,Label='MAPST')
call mma_allocate(MAPSP,NSS,Label='MAPSP')
call mma_allocate(MAPMS,NSS,Label='MAPMS')
ISS = 0
do ISTATE=1,NSTATE
  JOB = JBNUM(ISTATE)
  MPLET = MLTPLT(JOB)
  do MSPROJ=-MPLET+1,MPLET-1,2
    ISS = ISS+1
    MAPST(ISS) = ISTATE
    MAPSP(ISS) = MPLET
    MAPMS(ISS) = MSPROJ
  end do
end do

call mma_allocate(ISGS,NSS,Label='ISGS')
call mma_allocate(PNMRT,NTP,3,3,Label='PNMRT')
call mma_allocate(PNMR,NTP,3,3,Label='PNMR')
call mma_allocate(PNMRC,NTP,3,3,Label='PNMRC')
call mma_allocate(PNMRD,NTP,3,3,Label='PNMRD')
call mma_allocate(PNMRCPS,NTP,NSS,3,3,Label='PNMRCPS')
call mma_allocate(TMPf,NTP,Label='TMPf')
call mma_allocate(DIPSOf,3,NSS,NSS,Label='DIPSOf')
call mma_allocate(DIPSOfc,3,NSS,NSS,Label='DIPSOfc')
call mma_allocate(DIPSOfcsd,3,NSS,NSS,Label='DIPSOfcsd')
call mma_allocate(DIPSOfsd,3,NSS,NSS,Label='DIPSOfsd')
call mma_allocate(DIPSOfpso,3,NSS,NSS,Label='DIPSOfpso')
call mma_allocate(DIMSO,3,3,NSS,NSS,Label='DIMSO')
call mma_allocate(SPNSFS,3,NSS,NSS,Label='SPNSFS')
call mma_allocate(Z,NSS,NSS,Label='Z')
call mma_allocate(ZEKL,2,2,3,NSTATE,Label='ZEKL')
call mma_allocate(GCONT,3,3,NSTATE,Label='GCONT')

do IPROP=1,NPROP
  if ((PNAME(IPROP)(1:3) == 'ASD') .and. (ICOMP(IPROP) == 1)) then

    ! Get the center number
    read(PNAME(IPROP),'(4x,i4)') ICEN

    write(u6,*) '  ========================================='
    write(u6,*) '  A(Total)-Matrix for center:',ICEN
    write(u6,*) '  ========================================='

    write(SDPROP,'(a4,i4)') 'ASD ',ICEN
    write(u6,*) 'Looking for ',SDPROP
    write(PSOPROP,'(a4,i4)') 'PSOP',ICEN
    write(u6,*) 'Looking for ',PSOPROP

    ! Identify which properties are ASD matrix elements:

    !ccccccccccccccccccccccccccccccccccccccc
    ! Testing - use overlap matrix
    !ccccccccccccccccccccccccccccccccccccccc
    !SDPROP = 'MLTPL  0'
    !write(u6,*) 'Looking for overlap matrix ',SDPROP
    !do KPROP=1,NPROP
    !  if (PNAME(KPROP) == SDPROP) then
    !    IAMFI1 = KPROP
    !    IAMFI2 = 0
    !    IAMFI3 = 0
    !    IAMFI4 = 0
    !    IAMFI5 = 0
    !    IAMFI6 = 0
    !  end if
    !end do
    !ccccccccccccccccccccccccccccccccccccccc
    ! End testing - use overlap matrix
    !ccccccccccccccccccccccccccccccccccccccc

    ! These will hold the entire expression
    ! For g: L+2S
    ! For A: ?

    ! Identify which properties are Orbital Paramagnetic (PSOP) matrix elements:
    call Allocate_and_Load_PSOP()

    ! Labeled AMFI for now
    ! 1,2,3,4,5,6 -> xx,xy,xz,yy,yz,zz
    IAMFI1 = 0
    IAMFI2 = 0
    IAMFI3 = 0
    IAMFI4 = 0
    IAMFI5 = 0
    IAMFI6 = 0
    do KPROP=1,NPROP
      if ((PNAME(KPROP)(1:3) == SDPROP(1:3)) .and. (PNAME(KPROP)(5:8) == SDPROP(5:8))) then
        if (ICOMP(KPROP) == 1) IAMFI1 = KPROP
        if (ICOMP(KPROP) == 2) IAMFI2 = KPROP
        if (ICOMP(KPROP) == 3) IAMFI3 = KPROP
        if (ICOMP(KPROP) == 4) IAMFI4 = KPROP
        if (ICOMP(KPROP) == 5) IAMFI5 = KPROP
        if (ICOMP(KPROP) == 6) IAMFI6 = KPROP
      end if
    end do

    call Allocate_Z()

    do ISS=1,NSS
      ISTATE = MAPST(ISS)
      MPLET1 = MAPSP(ISS)
      MSPROJ1 = MAPMS(ISS)
      S1 = Half*real(MPLET1-1,kind=wp)
      SM1 = Half*real(MSPROJ1,kind=wp)
      do JSS=1,NSS
        JSTATE = MAPST(JSS)
        MPLET2 = MAPSP(JSS)
        MSPROJ2 = MAPMS(JSS)
        S2 = Half*real(MPLET2-1,kind=wp)
        SM2 = Half*real(MSPROJ2,kind=wp)
        AMFI1 = Zero
        AMFI2 = Zero
        AMFI3 = Zero
        AMFI4 = Zero
        AMFI5 = Zero
        AMFI6 = Zero

        ! SD
        if (IAMFI1 /= 0) AMFI1 = PROP(ISTATE,JSTATE,IAMFI1)
        if (IAMFI2 /= 0) AMFI2 = PROP(ISTATE,JSTATE,IAMFI2)
        if (IAMFI3 /= 0) AMFI3 = PROP(ISTATE,JSTATE,IAMFI3)
        if (IAMFI4 /= 0) AMFI4 = PROP(ISTATE,JSTATE,IAMFI4)
        if (IAMFI5 /= 0) AMFI5 = PROP(ISTATE,JSTATE,IAMFI5)
        if (IAMFI6 /= 0) AMFI6 = PROP(ISTATE,JSTATE,IAMFI6)

        ! Note that the 6th element contains r*r
        ! This is 4*pi * the contact term we want
        ! We need 8*pi/3 * < delta > so we just want
        ! a factor of 2/3
        ACNT = Zero

        if (IFACALFC) then
          ACNT = AMFI6/OneHalf
          !write(u6,*) 'ACNT',ACNT
        else if ((ISS == 1) .and. (JSS == 1)) then
          write(u6,*) '********************'
          write(u6,*) '* Skipping FC Part *'
          write(u6,*) '********************'
        end if

        ! Since this is a traceless tensor, generate the 6th
        ! element (3z^2-r^2)
        AMFI6 = -AMFI1-AMFI4

        if (IFACALSD) then
          AMFI1 = -AMFI1
          AMFI2 = -AMFI2
          AMFI3 = -AMFI3
          AMFI4 = -AMFI4
          AMFI5 = -AMFI5
          AMFI6 = -AMFI6
        else
          if ((ISS == 1) .and. (JSS == 1)) then
            write(u6,*) '********************'
            write(u6,*) '* Skipping SD Part *'
            write(u6,*) '********************'
          end if
          AMFI1 = Zero
          AMFI2 = Zero
          AMFI3 = Zero
          AMFI4 = Zero
          AMFI5 = Zero
          AMFI6 = Zero
        end if

        AMFI1 = AMFI1+ACNT
        AMFI4 = AMFI4+ACNT
        AMFI6 = AMFI6+ACNT

        if (.not. MAG_X2C) then
          call ADD_INFO('ASDFC1',[AMFI1],1,5)
          call ADD_INFO('ASDFC2',[AMFI2],1,5)
          call ADD_INFO('ASDFC3',[AMFI3],1,5)
          call ADD_INFO('ASDFC4',[AMFI4],1,5)
          call ADD_INFO('ASDFC5',[AMFI5],1,5)
          call ADD_INFO('ASDFC6',[AMFI6],1,5)
        end if

        !ccccccccccccccccccccccccccccccccccccccc
        ! Testing - use overlap matrix
        !ccccccccccccccccccccccccccccccccccccccc
        !AMFI1 = PROP(ISTATE,JSTATE,IAMFI1)
        !AMFI2 = Zero
        !AMFI3 = Zero
        !AMFI4 = AMFI1
        !AMFI5 = Zero
        !AMFI6 = AMFI1
        !ccccccccccccccccccccccccccccccccccccccc
        ! END Testing - use overlap matrix
        !ccccccccccccccccccccccccccccccccccccccc

        ! WIGNER-ECKART THEOREM:
        FACT = One/sqrt(real(MPLET1,kind=wp))
        if (MPLET1 == MPLET2-2) FACT = -FACT
        CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
        CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
        CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
        CGX = sqrt(Half)*(CGM-CGP)
        CGY = sqrt(Half)*(CGM+CGP)

        ZXYZR(ISS,JSS,1) = CGX*AMFI1+CG0*AMFI3
        ZXYZI(ISS,JSS,1) = CGY*AMFI2
        ZXYZR(ISS,JSS,2) = CGX*AMFI2+CG0*AMFI5
        ZXYZI(ISS,JSS,2) = CGY*AMFI4
        ZXYZR(ISS,JSS,3) = CGX*AMFI3+CG0*AMFI6
        ZXYZI(ISS,JSS,3) = CGY*AMFI5
      end do
    end do

    ZXYZI(:,:,1) = ZXYZI(:,:,1)+LXI(:,:)
    ZXYZI(:,:,2) = ZXYZI(:,:,2)+LYI(:,:)
    ZXYZI(:,:,3) = ZXYZI(:,:,3)+LZI(:,:)

    call Deallocate_PSOP()

    ! SVC 20090926 Experimental
    ! Add analysis of different contributions

    ! Establish which spin components of SFS belong to the ground state

    GSENERGY = ENERGY(1)
    do ISTATE=2,NSTATE
      if (ENERGY(ISTATE) < GSENERGY) GSENERGY = ENERGY(ISTATE)
    end do

    IMLTPL = 0
    do ISTATE=1,NSTATE
      ISGS(IMLTPL+1:IMLTPL+MLTPLT(JBNUM(ISTATE))) = abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp
      IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
    end do

    ! Analyze the different contributions to the GS Kramers doublet
    ! Zeeman matrix elements.  There are 4 different ME's: <1|Ze|1>,
    ! <1|Ze|2>, <2|Ze|1>, and <2|Ze|2>, stored in ZEKL.  Contributions
    ! of SFS i,j to SOS k,l (k,l=1,2): <k|Ze|l> = Sum(i,j) U(i,k)*
    ! <i|Ze|j> U(j,l).  This sum is decomposed into parts belonging to
    ! each SFS state i as follows:
    ! -> GS's contain only MEs with themselves and other GS's
    ! -> ES's contain MEs with themselves, the GS's (2x) and other ES's
    !    The ME's with the GS's are counted twice as they do not belong to
    !    any GS's (they contain only ME's within their own GS group)
    !    The contributions with other ES's are split between the ES's,
    !    counting them double (<i|Ze|j> and <j|Ze|i>) and divide by two later.

    ZEKL(:,:,:,:) = cZero

    IMLTPL = 0
    do ISTATE=1,NSTATE

      ISTART = IMLTPL+1
      IFINAL = IMLTPL+MLTPLT(JBNUM(ISTATE))

      if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

        ! Contribution of the GS spin components
        do IXYZ=1,3
          do ISS=ISTART,IFINAL
            do JSS=1,NSS
              if (ISGS(JSS)) then
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
              end if
            end do
          end do
        end do

      else

        ! Contributions of the ES spin components
        do IXYZ=1,3
          do ISS=ISTART,IFINAL
            do JSS=1,NSS
              call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
              call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
              if (ISGS(JSS)) then
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
              end if
              !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
            end do
          end do
        end do

      end if

      !do IXYZ=1,3
      !  write(u6,720) 'ZEKL',IXYZ,ISTATE,ZEKL(:,:,IXYZ,ISTATE)
      !end do

      IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
    end do

    ! We now have decomposed the <k|Ze|l> into terms belonging to either
    ! a GS or an ES for each k,l=1,2 and p=x,y,z stored in ZEKL(k,l,p,SFS)
    ! Now, these new decomposed terms of the ZEKL ME's are combined to
    ! form the G tensor.  Consider e.g. that <k|Ze|l> is decomposed into
    ! <k|GS|l> + <k|ES1|l> + <k|ES2|l>, then the contributions to G are given as:
    ! -> G_pq/2 = <k|Ze_p|l> <l|Ze_q|k>
    !         = (<k|GS_p|l> + <k|ES1_p|l> + <k|ES2_p|l>)
    !         * (<l|GS_q|k> + <l|ES1_q|k> + <l|ES2_q|k>)

    ! from GS: (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2
    ! from ES1: 2*((<k|ES1_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_q|l>/2 * <l|ES1_p|k>/2)/2)
    !           + (<k|ES1_p|l>/2 * <l|ES1_q|k>/2)/2 + (<k|ES1_q|l>/2 * <l|ES2_p|k>/2)/2
    !           + (<k|ES1_p|l>/2 * <l|ES2_q|k>/2)/2 + (<k|ES2_q|l>/2 * <l|ES1_p|k>/2)/2
    ! In the end, the outer division by 2 cancels on both sides, and the
    ! inner divisions by two combine to a division by 4.

    GCONT(:,:,:) = cZero
    GTOTAL(:,:) = Zero

    do ISTATE=1,NSTATE
      do IXYZ=1,3
        do JXYZ=1,3

          if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

            ! Contributions for the GS's
            do JSTATE=1,NSTATE
              if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) then
                do I=1,2
                  do J=1,2
                    GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                       ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                  end do
                end do
              end if
            end do

          else

            ! Contributions for the ES's
            do JSTATE=1,NSTATE
              do I=1,2
                do J=1,2
                  GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                     ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                  if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) &
                    GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                       ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                end do
              end do
            end do

          end if

        end do
      end do

      GTOTAL(:,:) = GTOTAL(:,:)+real(GCONT(:,:,ISTATE),kind=wp)
    end do

    ! Continue original calculation of G tensor (=gg^*)

    call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,1),ZXYZI(:,:,1))
    call PRCMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1))
    call MULMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1),eex,Z)
    DIPSOf(1,:,:) = Z(:,:)
    call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,2),ZXYZI(:,:,2))
    call MULMAT(NSS,ZXYZR(:,:,2),ZXYZI(:,:,2),eey,Z)
    DIPSOf(2,:,:) = Z(:,:)
    call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,3),ZXYZI(:,:,3))
    call MULMAT(NSS,ZXYZR(:,:,3),ZXYZI(:,:,3),eez,Z)
    DIPSOf(3,:,:) = Z(:,:)

    if (IFSONCINI) then

      DIMSO(:,:,:,:) = Zero
      write(DMPPROP,'(a4,i4)') 'DMP ',ICEN
      write(u6,*) 'Looking for ',DMPPROP
      do KPROP=1,NPROP
        if (PNAME(KPROP) == DMPPROP) then
          call mma_Allocate(SOPRR,NSS**2,1,Label='SOPRR')
          call mma_Allocate(SOPRI,NSS**2,1,Label='SOPRI')
          SOPRR(:,1) = Zero
          SOPRI(:,1) = Zero

          call SMMAT(PROP,SOPRR,NSS,KPROP,0)
          call ZTRNSF(NSS,USOR,USOI,SOPRR,SOPRI)
          !CALL PRCMAT(NSS,SOPRR,SOPRI)
          call MULMAT(NSS,SOPRR,SOPRI,eex,Z)
          ic = 0
          jc = 0
          !write(u6,*) 'ICOMP(KPROP)',ICOMP(KPROP)
          select case (ICOMP(KPROP))
            case (1)
              ic = 1
              jc = 1
            case (2)
              ic = 1
              jc = 2
            case (3)
              ic = 1
              jc = 3
            case (4)
              ic = 2
              jc = 1
            case (5)
              ic = 2
              jc = 2
            case (6)
              ic = 2
              jc = 3
            case (7)
              ic = 3
              jc = 1
            case (8)
              ic = 3
              jc = 2
            case (9)
              ic = 3
              jc = 3
          end select
          DIMSO(ic,jc,:,:) = 0
          do ISS=1,NSS
            !DIMSO(ICOMP(KPROP),ISS,ISS) = Z(ISS,ISS)
            DIMSO(ic,jc,ISS,ISS) = Z(ISS,ISS)
          end do
          call mma_deallocate(SOPRR)
          call mma_deallocate(SOPRI)
        end if

      end do
      PNMRT(:,:,:) = Zero
      PNMR(:,:,:) = Zero
      PNMRC(:,:,:) = Zero
      PNMRD(:,:,:) = Zero
      PNMRCPS(:,:,:,:) = Zero

      do iT=1,NTP
        if (iT == 1) then
          TMPf(iT) = TMINP+1.0e-4_wp
        else
          DLTTA = (TMAXP-TMINP)/real(NTP-1,kind=wp)
          TMPf(iT) = TMINP+DLTTA*real(iT-1,kind=wp)
        end if
        Zstat = Zero
        do Iss=1,Nss
          p_Boltz = exp(-ESO(Iss)/Boltz_k/TMPf(iT))
          !write(u6,*) 'p_Boltz',p_Boltz
          Zstat = Zstat+p_Boltz
          DiamT(:,:) = Zero
          CurieT(:,:) = Zero
          HFC_2(:,:) = Zero
          HFC_3(:,:) = Zero
          do Jss=1,Nss
            dlt_E = Eso(Iss)-Eso(Jss)
            HFC_1(:,:) = Zero
            do ic=1,3
              do jc=1,3

                HFC_1(ic,jc) = real(DIPSOm(ic,Iss,Jss)*conjg(DIPSOf(jc,Iss,Jss)))/(Boltz_k*TMPf(iT))

                !HFC_1(ic,jc) = real(DIPSOf(ic,Iss,Jss)*conjg(DIPSOm(jc,Iss,Jss)))/(Boltz_k*TMPf(iT))

                !HFC_3(ic,jc) = 10.e-6_wp*real(DIMSO(ic,jc,Iss,Jss))/(ALPHA2*auTocm)

                !if(ABS(dlt_E) < 1.0e-3_wp) then
                if (abs(dlt_E) < 10.97_wp) then ! what is this number?
                  CurieT(ic,jc) = CurieT(ic,jc)+HFC_1(ic,jc)
                  HFC_2(ic,jc) = HFC_2(ic,jc)+HFC_1(ic,jc)
                  !HFC_3(ic,jc) = HFC_3(ic,jc)+Zero*HFC_1(ic,jc)
                  !DiamT(ic,jc) = DiamT(ic,jc)+Zero*HFC_1(ic,jc)
                else

                  HFC_2(ic,jc) = HFC_2(ic,jc)-(real(conjg(DIPSOm(ic,Iss,Jss))*DIPSOf(jc,Iss,Jss))+ &
                                               real(DIPSOm(ic,Iss,Jss)*conjg(DIPSOf(jc,Iss,Jss))))/dlt_E
                  !HFC_2(ic,jc) = HFC_2(ic,jc)-20.0_wp*(real(conjg(DIPSOm(ic,Iss,Jss))*DIPSOf(jc,Iss,Jss)))/dlt_E

                  HFC_3(ic,jc) = HFC_3(ic,jc)-(real(conjg(DIPSOm(ic,Iss,Jss))*DIPSOf(jc,Iss,Jss))+ &
                                               real(DIPSOm(ic,Iss,Jss)*conjg(DIPSOf(jc,Iss,Jss))))/dlt_E

                  !CurieT(ic,jc) = CurieT(ic,jc)-Zero*(real(conjg(DIPSOm(ic,Iss,Jss))*DIPSOf(jc,Iss,Jss))+ &
                  !                                    real(DIPSOm(ic,Iss,Jss)*conjg(DIPSOf(jc,Iss,Jss))))/dlt_E

                  !DiamT(ic,jc) = DiamT(ic,jc)-Zero*(real(conjg(DIPSOm(ic,Iss,Jss))*DIPSOf(jc,Iss,Jss))+ &
                  !                                  real(DIPSOm(ic,Iss,Jss)*conjg(DIPSOf(jc,Iss,Jss))))/dlt_E

                end if
              end do
            end do
            HFC_2(:,:) = HFC_2(:,:)+1.0e-6_wp*real(DIMSO(:,:,Iss,Jss))/(ALPHA2*auTocm)
            DiamT(:,:) = DiamT(:,:)+1.0e-6_wp*real(DIMSO(:,:,Iss,Jss))/(ALPHA2*auTocm)
          end do !Jss
          PNMRT(iT,:,:) = PNMRT(iT,:,:)+p_Boltz*HFC_2(:,:)

          PNMR(iT,:,:) = PNMR(iT,:,:)+p_Boltz*HFC_3(:,:)

          PNMRC(iT,:,:) = PNMRC(iT,:,:)+p_Boltz*CurieT(:,:)

          PNMRCPS(iT,Iss,:,:) = PNMRCPS(iT,Iss,:,:)+CurieT(:,:)

          PNMRD(iT,:,:) = PNMRD(iT,:,:)+p_Boltz*DiamT(:,:)
        end do !Iss
        PNMRT(iT,:,:) = 1.0e6_wp*auTocm*ALPHA2*(PNMRT(iT,:,:)/Zstat)
        PNMR(iT,:,:) = 1.0e6_wp*auTocm*ALPHA2*(PNMR(iT,:,:)/Zstat)
        PNMRC(iT,:,:) = 1.0e6_wp*auTocm*ALPHA2*(PNMRC(iT,:,:)/Zstat)
        PNMRD(iT,:,:) = 1.0e6_wp*auTocm*ALPHA2*(PNMRD(iT,:,:)/Zstat)
        PNMRCPS(iT,:,:,:) = 1.0e6_wp*auTocm*ALPHA2*PNMRCPS(iT,:,:,:)

      end do ! iT

      !PNMRT(:,:,:) = 1.0e6_wp*PNMRT(:,:,:)

      write(u6,'(/)')
      write(u6,'(A)') repeat('-',120)
      write(u6,'(30X,A)') ' TOTAL SONCINI PNMR TENSOR in (ppm)'
      write(u6,'(A)') repeat('-',120)
      write(u6,*)
      !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(u6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
      write(u6,*)
      do iT=1,NTP
        write(u6,'(4X,F6.1,3X,11(F15.4,2X),F8.4)') TMPf(iT),((PNMRT(iT,ic,jc),jc=1,3),ic=1,3)
      end do
      !write(u6,*) '*******'
      !do iT=1,NTP
      !  do ic=1,3
      !    do jc=1,3
      !      write(u6,*) PNMRT_tens(iT,ic,jc)
      !    end do
      !  end do
      !end do
      !write(u6,*) '*******'

      write(u6,'(/)')
      write(u6,'(A)') repeat('-',120)
      write(u6,'(30X,A)') 'paramagnetic SONCINI PNMR TENSOR in (ppm)'
      write(u6,'(A)') repeat('-',120)
      write(u6,*)
      !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(u6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
      write(u6,*)
      do iT=1,NTP
        write(u6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPf(iT),((PNMR(iT,ic,jc),jc=1,3),ic=1,3)
      end do

      write(u6,'(/)')
      write(u6,'(A)') repeat('-',120)
      write(u6,'(30X,A)') 'Diamagnetic  SONCINI PNMR TENSOR in (ppm)'
      write(u6,'(A)') repeat('-',120)
      write(u6,*)
      !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(u6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
      write(u6,*)
      do iT=1,NTP
        write(u6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPf(iT),((PNMRD(iT,ic,jc),jc=1,3),ic=1,3)
      end do

      write(u6,'(/)')
      write(u6,'(A)') repeat('-',120)
      write(u6,'(30X,A)') 'Curie term contrib. to SONCINI PNMR TENSOR in (ppm)'
      write(u6,'(A)') repeat('-',120)
      write(u6,*)
      !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
      write(u6,'(6X,A,9X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
      write(u6,*)
      do iT=1,NTP
        write(u6,'(4X,F6.1,3X,11(F15.4,2X),F8.4)') TMPf(iT),((PNMRC(iT,ic,jc),jc=1,3),ic=1,3)
        !if (.true.) then
        write(u6,'(/)')
        write(u6,'(A)') repeat('-',120)
        write(u6,'(30X,A)') 'Curie per state in (ppm)'
        write(u6,'(A)') repeat('-',120)
        write(u6,*)
        !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
        write(u6,'(6X,A,8X,9(A,9X))') 'NSS','xx','xy','xz','yx','yy','yz','zx','zy','zz'
        write(u6,*)
        !write(u6,*)
        !write(u6,*) PNMRC(iT,3,3)
        !write(u6,*)
        do Iss=1,NSS
          write(u6,'(4X,I5,3X,11(F10.4,2X),F8.4)') Iss,((PNMRCPS(iT,Iss,ic,jc),jc=1,3),ic=1,3)
        end do

        ! end if
        write(u6,'(/)')
        write(u6,'(A)') repeat('=',130)
        write(u6,'(/)')
      end do
      write(u6,*)
      write(u6,*)
      write(u6,*) '  A-Matrix'
      write(u6,*) '  =========='
    end if !IFSONCINI

    ISS = 1
    do while ((ISS <= NSS) .and. (ENSOR(min(ISS,NSS))-ENSOR(1) <= EPRATHR))

      GTENS(:,:) = Zero

      KDGN = 1
      do JSS=ISS+1,NSS
        EDIFF = ENSOR(JSS)-ENSOR(ISS)
        if (IFATCALSA .and. IFGTSHSA) then
          KDGN = MULTIP
          !write(u6,*) 'KDGN=',KDGN
        else if (abs(EDIFF) < 1.0e-6_wp) then
          KDGN = KDGN+1
        end if
      end do

      write(u6,*)
      do I=1,KDGN
        write(u6,'(3x,A9,I4,3x,A4,F18.8)') 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
      end do
      write(u6,'(3x,A46)') '----------------------------------------------'
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (IFATCALSA) then
        if (ISS == 1) IFUNCT = 0
        call SINANI(KDGN,IFUNCT,NSS,DIPSOf,SPNSFS)
        IFUNCT = IFUNCT+KDGN
      end if
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (KDGN /= 2) then
        write(u6,*) 'no twofold degeneracy'
      else

        if ((ISS == 1) .and. (KDGN == 2) .and. (IPGLOB >= 3)) then
          write(u6,*) 'Experimental: SFS contributions to G=gg+'
          write(u6,*)
          write(u6,'(a6,9(5x,a2,5x))') 'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
          write(u6,*)
          do ISTATE=1,NSTATE
            write(u6,'(2x,I2,2x,9(F12.6))') ISTATE,real(GCONT(:,:,ISTATE))
          end do

          write(u6,*)
          write(u6,'(A6,9(F12.6))') 'total ',GTOTAL(:,:)
        end if

        JSS = ISS+1

        do IXYZ=1,3
          do JXYZ=1,3
            GTIJ = Zero
            do ISO=ISS,JSS
              GTIJ = GTIJ+sum(ZXYZR(ISO,ISS:JSS,IXYZ)*ZXYZR(ISS:JSS,ISO,JXYZ)-ZXYZI(ISO,ISS:JSS,IXYZ)*ZXYZI(ISS:JSS,ISO,JXYZ))
            end do
            GTENS(IXYZ,JXYZ) = Two*GTIJ
          end do
        end do

        !if (IPGLOB > 3) then
        write(u6,*) 'G tensor = gg+'
        write(u6,*)
        write(u6,'(6x,3(6x,a2,4x))') (xyzchr(IXYZ),IXYZ=1,3)
        do IXYZ=1,3
          write(u6,'(2x,a2,2x,3(1x,f18.8,1x))') xyzchr(IXYZ),GTENS(IXYZ,:)
        end do
        !end if

        EVR(:) = Zero
        EVI(:) = Zero
        TMPMAT(:,:) = GTENS(:,:)
        call unitmat(TMPVEC,3)

        call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

        ! construct g_s matrix from G by back-transormation of the
        ! square root of the G eigenvalues
        do IXYZ=1,3
          do JXYZ=1,3
            GTENS(IXYZ,JXYZ) = sum(TMPVEC(IXYZ,:)*sqrt(EVR(:))*TMPVEC(JXYZ,:))
          end do
        end do

        write(u6,'(6x,3(5x,a2,5x),4x,4x,2x,8x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('a_',IXYZ,IXYZ=1,3)
        write(u6,*)
        do IXYZ=1,3
          write(u6,'(2x,a2,2x,3(1x,f10.6,1x),4x,a2,i1,a1,2x,f12.9,2x,a2,2x,3(1x,f8.4,1x))') xyzchr(IXYZ),GTENS(IXYZ,:),'a_',IXYZ, &
                                                                                            ':',sqrt(EVR(IXYZ)),xyzchr(IXYZ), &
                                                                                            TMPVEC(IXYZ,:)
        end do

        call ADD_INFO('ATENS',GTENS,9,5)
        call ADD_INFO('ATENS2',EVR,3,5)

      end if

      ISS = ISS+KDGN

    end do

    call Deallocate_Z()

    if (IFACALFCON) then

      write(u6,*)
      write(u6,*) '  ========================================='
      write(u6,*) '  A (FC)-Matrix for center:',ICEN
      write(u6,*) '  ========================================='

      call Allocate_and_Load_PSOP()

      call mma_allocate(MXR,NSS,NSS,Label='MXR')
      call mma_allocate(MXI,NSS,NSS,Label='MXI')
      MXR(:,:) = Zero
      MXI(:,:) = Zero
      call mma_allocate(MYR,NSS,NSS,Label='MYR')
      call mma_allocate(MYI,NSS,NSS,Label='MYI')
      MYR(:,:) = Zero
      MYI(:,:) = Zero
      call mma_allocate(MZR,NSS,NSS,Label='MZR')
      call mma_allocate(MZI,NSS,NSS,Label='MZI')
      MZR(:,:) = Zero
      MZI(:,:) = Zero
      call SMMAT(PROP,MXR,NSS,0,1)
      call SMMAT(PROP,MYI,NSS,0,2)
      call SMMAT(PROP,MZR,NSS,0,3)

      MXR(:,:) = FEGVAL*MXR(:,:)
      MZR(:,:) = FEGVAL*MZR(:,:)

      MXI(:,:) = MXI(:,:)+LXI(:,:)
      MYI(:,:) = FEGVAL*MYI(:,:)+LYI(:,:)
      MZI(:,:) = MZI(:,:)+LZI(:,:)

      call Deallocate_PSOP()

      call ZTRNSF(NSS,USOR,USOI,MXR,MXI)
      call ZTRNSF(NSS,USOR,USOI,MYR,MYI)
      call ZTRNSF(NSS,USOR,USOI,MZR,MZI)

      call mma_deallocate(MXR)
      call mma_deallocate(MXI)
      call mma_deallocate(MYR)
      call mma_deallocate(MYI)
      call mma_deallocate(MZR)
      call mma_deallocate(MZI)

      IAMFI1 = 0
      IAMFI2 = 0
      IAMFI3 = 0
      IAMFI4 = 0
      IAMFI5 = 0
      IAMFI6 = 0
      write(SDPROP,'(a4,i4)') 'ASD ',ICEN
      write(u6,*) 'Looking for ',SDPROP
      do KPROP=1,NPROP
        if ((PNAME(KPROP)(1:3) == SDPROP(1:3)) .and. (PNAME(KPROP)(5:8) == SDPROP(5:8))) then
          if (ICOMP(KPROP) == 1) IAMFI1 = KPROP
          if (ICOMP(KPROP) == 2) IAMFI2 = KPROP
          if (ICOMP(KPROP) == 3) IAMFI3 = KPROP
          if (ICOMP(KPROP) == 4) IAMFI4 = KPROP
          if (ICOMP(KPROP) == 5) IAMFI5 = KPROP
          if (ICOMP(KPROP) == 6) IAMFI6 = KPROP
        end if
      end do

      call Allocate_Z()

      do ISS=1,NSS
        ISTATE = MAPST(ISS)
        MPLET1 = MAPSP(ISS)
        MSPROJ1 = MAPMS(ISS)
        S1 = Half*real(MPLET1-1,kind=wp)
        SM1 = Half*real(MSPROJ1,kind=wp)
        do JSS=1,NSS
          JSTATE = MAPST(JSS)
          MPLET2 = MAPSP(JSS)
          MSPROJ2 = MAPMS(JSS)
          S2 = Half*real(MPLET2-1,kind=wp)
          SM2 = Half*real(MSPROJ2,kind=wp)
          AMFI1 = Zero
          AMFI2 = Zero
          AMFI3 = Zero
          AMFI4 = Zero
          AMFI5 = Zero
          AMFI6 = Zero
          ! SD
          if (IAMFI1 /= 0) AMFI1 = PROP(ISTATE,JSTATE,IAMFI1)
          if (IAMFI2 /= 0) AMFI2 = PROP(ISTATE,JSTATE,IAMFI2)
          if (IAMFI3 /= 0) AMFI3 = PROP(ISTATE,JSTATE,IAMFI3)
          if (IAMFI4 /= 0) AMFI4 = PROP(ISTATE,JSTATE,IAMFI4)
          if (IAMFI5 /= 0) AMFI5 = PROP(ISTATE,JSTATE,IAMFI5)
          if (IAMFI6 /= 0) AMFI6 = PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = Zero

          !if (IFACALFC) then
          ACNT = AMFI6/OneHalf
          !else if ((ISS == 1) .and. (JSS == 1)) then
          !  write(u6,*) '********************'
          !  write(u6,*) '* Skipping FC Part *'
          !  write(u6,*) '********************'
          !end if
          !AMFI6 = -AMFI1-AMFI4

          !if (IFACALSD) then
          !  AMFI1 = -AMFI1
          !  AMFI2 = -AMFI2
          !  AMFI3 = -AMFI3
          !  AMFI4 = -AMFI4
          !  AMFI5 = -AMFI5
          !  AMFI6 = -AMFI6
          !else
          !  if ((ISS == 1) .and. (JSS == 1)) then
          !    write(u6,*) '********************'
          !    write(u6,*) '* Skipping SD Part *'
          !    write(u6,*) '********************'
          !  end if
          AMFI1 = Zero
          AMFI2 = Zero
          AMFI3 = Zero
          AMFI4 = Zero
          AMFI5 = Zero
          AMFI6 = Zero
          !end if

          AMFI1 = AMFI1+ACNT
          AMFI4 = AMFI4+ACNT
          AMFI6 = AMFI6+ACNT

          ! WIGNER-ECKART THEOREM:
          FACT = One/sqrt(real(MPLET1,kind=wp))
          if (MPLET1 == MPLET2-2) FACT = -FACT
          CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
          CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
          CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
          CGX = sqrt(Half)*(CGM-CGP)
          CGY = sqrt(Half)*(CGM+CGP)

          ZXYZR(ISS,JSS,1) = CGX*AMFI1+CG0*AMFI3
          ZXYZI(ISS,JSS,1) = CGY*AMFI2
          ZXYZR(ISS,JSS,2) = CGX*AMFI2+CG0*AMFI5
          ZXYZI(ISS,JSS,2) = CGY*AMFI4
          ZXYZR(ISS,JSS,3) = CGX*AMFI3+CG0*AMFI6
          ZXYZI(ISS,JSS,3) = CGY*AMFI5
        end do
      end do

      GSENERGY = ENERGY(1)
      do ISTATE=2,NSTATE
        if (ENERGY(ISTATE) < GSENERGY) GSENERGY = ENERGY(ISTATE)
      end do

      IMLTPL = 0
      do ISTATE=1,NSTATE
        ISGS(IMLTPL+1:IMLTPL+MLTPLT(JBNUM(ISTATE))) = abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp
        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do
      ZEKL(:,:,:,:) = cZero

      IMLTPL = 0
      do ISTATE=1,NSTATE

        ISTART = IMLTPL+1
        IFINAL = IMLTPL+MLTPLT(JBNUM(ISTATE))

        if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

          ! Contribution of the GS spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                  !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
                end if
              end do
            end do
          end do

        else

          ! Contributions of the ES spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                end if

                !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS, SS,ZEKL(:,:,IXYZ,ISTATE)
              end do
            end do
          end do

        end if

        !do IXYZ=1,3
        !  write(u6,720) 'ZEKL',IXYZ,ISTATE,ZEKL(:,:,IXYZ,ISTATE)
        !end do

        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do

      GCONT(:,:,:) = cZero
      GTOTAL(:,:) = Zero

      do ISTATE=1,NSTATE
        do IXYZ=1,3
          do JXYZ=1,3

            if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

              ! Contributions for the GS's
              do JSTATE=1,NSTATE
                if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) then
                  do I=1,2
                    do J=1,2
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    end do
                  end do
                end if
              end do

            else

              ! Contributions for the ES's
              do JSTATE=1,NSTATE
                do I=1,2
                  do J=1,2
                    GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                       ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) &
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                  end do
                end do
              end do

            end if

          end do
        end do

        GTOTAL(:,:) = GTOTAL(:,:)+real(GCONT(:,:,ISTATE))
      end do

      ! Continue original calculation of G tensor (=gg^*)

      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,1),ZXYZI(:,:,1))
      call MULMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1),eex,Z)
      DIPSOfc(1,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,2),ZXYZI(:,:,2))
      call MULMAT(NSS,ZXYZR(:,:,2),ZXYZI(:,:,2),eey,Z)
      DIPSOfc(2,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,3),ZXYZI(:,:,3))
      call MULMAT(NSS,ZXYZR(:,:,3),ZXYZI(:,:,3),eez,Z)
      DIPSOfc(3,:,:) = Z(:,:)

      ISS = 1
      do while ((ENSOR(min(ISS,NSS))-ENSOR(1) <= EPRATHR) .and. (ISS <= NSS))

        GTENS(:,:) = Zero

        KDGN = 1
        do JSS=ISS+1,NSS
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (IFATCALSA .and. IFGTSHSA) then
            KDGN = MULTIP
            !write(u6,*) 'KDGN=',KDGN
          else if (abs(EDIFF) < 1.0e-6_wp) then
            KDGN = KDGN+1
          end if
        end do

        write(u6,*)
        do I=1,KDGN
          write(u6,'(3x,A9,I4,3x,A4,F18.8)') 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
        end do
        write(u6,'(3x,A46)') '----------------------------------------------'
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (IFATCALSA) then
          if (ISS == 1) IFUNCT = 0
          call SINANI(KDGN,IFUNCT,NSS,DIPSOfc,SPNSFS)
          IFUNCT = IFUNCT+KDGN
        end if
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (KDGN /= 2) then
          write(u6,*) 'no twofold degeneracy'
        else

          if ((ISS == 1) .and. (KDGN == 2) .and. (IPGLOB >= 3)) then
            write(u6,*) 'Experimental: SFS contributions to G=gg+'
            write(u6,*)
            write(u6,'(a6,9(5x,a2,5x))') 'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
            write(u6,*)
            do ISTATE=1,NSTATE
              write(u6,'(2x,I2,2x,9(F12.6))') ISTATE,real(GCONT(:,:,ISTATE))
            end do

            write(u6,*)
            write(u6,'(A6,9(F12.6))') 'total ',GTOTAL(:,:)
          end if

          JSS = ISS+1

          do IXYZ=1,3
            do JXYZ=1,3
              GTIJ = Zero
              do ISO=ISS,JSS
                GTIJ = GTIJ+sum(ZXYZR(ISO,ISS:JSS,IXYZ)*ZXYZR(ISS:JSS,ISO,JXYZ)-ZXYZI(ISO,ISS:JSS,IXYZ)*ZXYZI(ISS:JSS,ISO,JXYZ))
              end do
              GTENS(IXYZ,JXYZ) = Two*GTIJ
            end do
          end do

          !if (IPGLOB > 3) then
          write(u6,*) 'G tensor = gg+'
          write(u6,*)
          write(u6,'(6x,3(6x,a2,4x))') (xyzchr(IXYZ),IXYZ=1,3)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f18.8,1x))') xyzchr(IXYZ),GTENS(IXYZ,:)
          end do
          !end if

          EVR(:) = Zero
          EVI(:) = Zero
          TMPMAT(:,:) = GTENS(:,:)
          call unitmat(TMPVEC,3)

          call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

          do IXYZ=1,3
            do JXYZ=1,3
              do KXYZ=1,3
                GTENS(IXYZ,JXYZ) = sum(TMPVEC(IXYZ,:)*sqrt(EVR(:))*TMPVEC(JXYZ,:))
              end do
            end do
          end do

          write(u6,'(6x,3(5x,a2,5x),4x,4x,2x,8x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('a_',IXYZ,IXYZ=1,3)
          write(u6,*)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f10.6,1x),4x,a2,i1,a1,2x,f12.9,2x,a2,2x,3(1x,f8.4,1x))') xyzchr(IXYZ),GTENS(IXYZ,:),'a_', &
                                                                                              IXYZ,':',sqrt(EVR(IXYZ)), &
                                                                                              xyzchr(IXYZ),TMPVEC(IXYZ,:)
          end do

        end if

        ISS = ISS+KDGN

      end do

      call Deallocate_Z()

    end if

    if (IFACALSDON) then

      write(u6,*)
      write(u6,*) '  ========================================='
      write(u6,*) '  A (SD)-Matrix for center:',ICEN
      write(u6,*) '  ========================================='

      IAMFI1 = 0
      IAMFI2 = 0
      IAMFI3 = 0
      IAMFI4 = 0
      IAMFI5 = 0
      IAMFI6 = 0
      write(SDPROP,'(a4,i4)') 'ASD ',ICEN
      write(u6,*) 'Looking for ',SDPROP
      do KPROP=1,NPROP
        if ((PNAME(KPROP)(1:3) == SDPROP(1:3)) .and. (PNAME(KPROP)(5:8) == SDPROP(5:8))) then
          if (ICOMP(KPROP) == 1) IAMFI1 = KPROP
          if (ICOMP(KPROP) == 2) IAMFI2 = KPROP
          if (ICOMP(KPROP) == 3) IAMFI3 = KPROP
          if (ICOMP(KPROP) == 4) IAMFI4 = KPROP
          if (ICOMP(KPROP) == 5) IAMFI5 = KPROP
          if (ICOMP(KPROP) == 6) IAMFI6 = KPROP
        end if
      end do

      call Allocate_Z()

      do ISS=1,NSS
        ISTATE = MAPST(ISS)
        MPLET1 = MAPSP(ISS)
        MSPROJ1 = MAPMS(ISS)
        S1 = Half*real(MPLET1-1,kind=wp)
        SM1 = Half*real(MSPROJ1,kind=wp)
        do JSS=1,NSS
          JSTATE = MAPST(JSS)
          MPLET2 = MAPSP(JSS)
          MSPROJ2 = MAPMS(JSS)
          S2 = Half*real(MPLET2-1,kind=wp)
          SM2 = Half*real(MSPROJ2,kind=wp)
          AMFI1 = Zero
          AMFI2 = Zero
          AMFI3 = Zero
          AMFI4 = Zero
          AMFI5 = Zero
          AMFI6 = Zero

          ! SD
          if (IAMFI1 /= 0) AMFI1 = PROP(ISTATE,JSTATE,IAMFI1)
          if (IAMFI2 /= 0) AMFI2 = PROP(ISTATE,JSTATE,IAMFI2)
          if (IAMFI3 /= 0) AMFI3 = PROP(ISTATE,JSTATE,IAMFI3)
          if (IAMFI4 /= 0) AMFI4 = PROP(ISTATE,JSTATE,IAMFI4)
          if (IAMFI5 /= 0) AMFI5 = PROP(ISTATE,JSTATE,IAMFI5)
          if (IAMFI6 /= 0) AMFI6 = PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = Zero

          !if (IFACALFC) then
          !  ACNT = AMFI6/OneHalf
          !else if ((ISS == 1) .and. (JSS == 1)) then
          !  write(u6,*) '********************'
          !  write(u6,*) '* Skipping FC Part *'
          !  write(u6,*) '********************'
          !end if
          AMFI6 = -AMFI1-AMFI4

          !if (IFACALSD) then
          AMFI1 = -AMFI1
          AMFI2 = -AMFI2
          AMFI3 = -AMFI3
          AMFI4 = -AMFI4
          AMFI5 = -AMFI5
          AMFI6 = -AMFI6
          !else
          !  if ((ISS == 1) .and. (JSS == 1)) then
          !    write(u6,*) '********************'
          !    write(u6,*) '* Skipping SD Part *'
          !    write(u6,*) '********************'
          !  end if
          !  AMFI1 = Zero
          !  AMFI2 = Zero
          !  AMFI3 = Zero
          !  AMFI4 = Zero
          !  AMFI5 = Zero
          !  AMFI6 = Zero
          !end if

          AMFI1 = AMFI1+ACNT
          AMFI4 = AMFI4+ACNT
          AMFI6 = AMFI6+ACNT

          ! WIGNER-ECKART THEOREM:
          FACT = One/sqrt(real(MPLET1,kind=wp))
          if (MPLET1 == MPLET2-2) FACT = -FACT
          CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
          CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
          CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
          CGX = sqrt(Half)*(CGM-CGP)
          CGY = sqrt(Half)*(CGM+CGP)

          ZXYZR(ISS,JSS,1) = CGX*AMFI1+CG0*AMFI3
          ZXYZI(ISS,JSS,1) = CGY*AMFI2
          ZXYZR(ISS,JSS,2) = CGX*AMFI2+CG0*AMFI5
          ZXYZI(ISS,JSS,2) = CGY*AMFI4
          ZXYZR(ISS,JSS,3) = CGX*AMFI3+CG0*AMFI6
          ZXYZI(ISS,JSS,3) = CGY*AMFI5
        end do
      end do

      GSENERGY = ENERGY(1)
      do ISTATE=2,NSTATE
        if (ENERGY(ISTATE) < GSENERGY) GSENERGY = ENERGY(ISTATE)
      end do

      IMLTPL = 0
      do ISTATE=1,NSTATE
        ISGS(IMLTPL+1:IMLTPL+MLTPLT(JBNUM(ISTATE))) = abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp
        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do

      ZEKL(:,:,:,:) = cZero

      IMLTPL = 0
      do ISTATE=1,NSTATE

        ISTART = IMLTPL+1
        IFINAL = IMLTPL+MLTPLT(JBNUM(ISTATE))

        if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

          ! Contribution of the GS spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                  !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
                end if
              end do
            end do
          end do

        else

          ! Contributions of the ES spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                end if

                !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
              end do
            end do
          end do

        end if

        !do IXYZ=1,3
        !  write(u6,720) 'ZEKL',IXYZ,ISTATE,ZEKL(:,:,IXYZ,ISTATE)
        !end do

        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do

      GCONT(:,:,:) = cZero
      GTOTAL(:,:) = Zero

      do ISTATE=1,NSTATE
        do IXYZ=1,3
          do JXYZ=1,3

            if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

              ! Contributions for the GS's
              do JSTATE=1,NSTATE
                if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) then
                  do I=1,2
                    do J=1,2
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                 ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    end do
                  end do
                end if
              end do

            else

              ! Contributions for the ES's
              do JSTATE=1,NSTATE
                do I=1,2
                  do J=1,2
                    GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                       ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) &
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                  end do
                end do
              end do

            end if

          end do
        end do

        GTOTAL(:,:) = GTOTAL(:,:)+real(GCONT(:,:,ISTATE))
      end do

      ! Continue original calculation of G tensor (=gg^*)

      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,1),ZXYZI(:,:,1))
      call MULMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1),eex,Z)
      DIPSOfsd(1,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,2),ZXYZI(:,:,2))
      call MULMAT(NSS,ZXYZR(:,:,2),ZXYZI(:,:,2),eey,Z)
      DIPSOfsd(2,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,3),ZXYZI(:,:,3))
      call MULMAT(NSS,ZXYZR(:,:,3),ZXYZI(:,:,3),eez,Z)
      DIPSOfsd(3,:,:) = Z(:,:)

      ISS = 1
      do while ((ENSOR(min(ISS,NSS))-ENSOR(1) <= EPRATHR) .and. (ISS <= NSS))

        GTENS(:,:) = Zero

        KDGN = 1
        do JSS=ISS+1,NSS
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (IFATCALSA .and. IFGTSHSA) then
            KDGN = MULTIP
          else if (abs(EDIFF) < 1.0e-6_wp) then
            KDGN = KDGN+1
          end if
        end do

        write(u6,*)
        do I=1,KDGN
          write(u6,'(3x,A9,I4,3x,A4,F18.8)') 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
        end do
        write(u6,'(3x,A46)') '----------------------------------------------'

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (IFATCALSA) then
          if (ISS == 1) IFUNCT = 0
          call SINANI(KDGN,IFUNCT,NSS,DIPSOfsd,SPNSFS)
          IFUNCT = IFUNCT+KDGN
        end if
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (KDGN /= 2) then
          write(u6,*) 'no twofold degeneracy'
        else

          if ((ISS == 1) .and. (KDGN == 2) .and. (IPGLOB >= 3)) then
            write(u6,*) 'Experimental: SFS contributions to G=gg+'
            write(u6,*)
            write(u6,'(a6,9(5x,a2,5x))') 'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
            write(u6,*)
            do ISTATE=1,NSTATE
              write(u6,'(2x,I2,2x,9(F12.6))') ISTATE,real(GCONT(:,:,ISTATE))
            end do

            write(u6,*)
            write(u6,'(A6,9(F12.6))') 'total ',GTOTAL(:,:)
          end if

          JSS = ISS+1

          do IXYZ=1,3
            do JXYZ=1,3
              GTIJ = Zero
              do ISO=ISS,JSS
                GTIJ = GTIJ+sum(ZXYZR(ISO,ISS:JSS,IXYZ)*ZXYZR(ISS:JSS,ISO,JXYZ)-ZXYZI(ISO,ISS:JSS,IXYZ)*ZXYZI(ISS:JSS,ISO,JXYZ))
              end do
              GTENS(IXYZ,JXYZ) = Two*GTIJ
            end do
          end do

          !if (IPGLOB > 3) then
          write(u6,*) 'G tensor = gg+'
          write(u6,*)
          write(u6,'(6x,3(6x,a2,4x))') (xyzchr(IXYZ),IXYZ=1,3)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f18.8,1x))') xyzchr(IXYZ),GTENS(IXYZ,:)
          end do
          !end if

          EVR(:) = Zero
          EVI(:) = Zero
          TMPMAT(:,:) = GTENS(:,:)
          call unitmat(TMPVEC,3)

          call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

          do IXYZ=1,3
            do JXYZ=1,3
              GTENS(IXYZ,JXYZ) = sum(TMPVEC(IXYZ,:)*sqrt(EVR(:))*TMPVEC(JXYZ,:))
            end do
          end do

          !do IXYZ=1,3
          !  write(u6,*) 'IXYZ',IXYZ
          !  write(u6,*)
          !  do JXYZ=1,3
          !    write(u6,*) 'JXYZ',JXYZ
          !    write(u6,*)
          !    GTENS(IXYZ,JXYZ) = Zero
          !    dO KXYZ=1,3
          !      write(u6,*) 'KXYZ',KXYZ
          !      write(u6,*)
          !      write(u6,*) 'GTENS(IXYZ,JXYZ)',GTENS(IXYZ,JXYZ),IXYZ,JXYZ
          !      write(u6,*)
          !      write(u6,*) 'TMPVEC(IXYZ,KXYZ)',TMPVEC(IXYZ,KXYZ),IXYZ,KXYZ
          !      write(u6,*)
          !      write(u6,*) 'SQRT(EVR(KXYZ))',SQRT(EVR(KXYZ)),KXYZ
          !      write(u6,*)
          !      write(u6,*) 'TMPVEC(JXYZ,KXYZ)',TMPVEC(JXYZ,KXYZ),JXYZ,KXYZ
          !      write(u6,*) '*************************'
          !      GTENS(IXYZ,JXYZ) = GTENS(IXYZ,JXYZ)+TMPVEC(IXYZ,KXYZ)*SQRT(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
          !      write(u6,*) 'GTENS(IXYZ,JXYZ)',GTENS(IXYZ,JXYZ),IXYZ,JXYZ
          !      write(u6,*) '*************************'
          !    end do
          !  end do
          !end do

          write(u6,'(6x,3(5x,a2,5x),4x,4x,2x,8x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('a_',IXYZ,IXYZ=1,3)
          write(u6,*)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f10.6,1x),4x,a2,i1,a1,2x,f12.9,2x,a2,2x,3(1x,f8.4,1x))') xyzchr(IXYZ),GTENS(IXYZ,:),'a_', &
                                                                                              IXYZ,':',sqrt(EVR(IXYZ)), &
                                                                                              xyzchr(IXYZ),TMPVEC(IXYZ,:)
          end do

        end if

        ISS = ISS+KDGN

      end do

      call Deallocate_Z()

    end if

    ! Skip if not a hyperfine calculation
    if (IFACALFCSDON) then

      write(u6,*)
      write(u6,*) '  ========================================='
      write(u6,*) '  A (FC+SD)-Matrix for center:',ICEN
      write(u6,*) '  ========================================='

      IAMFI1 = 0
      IAMFI2 = 0
      IAMFI3 = 0
      IAMFI4 = 0
      IAMFI5 = 0
      IAMFI6 = 0
      write(SDPROP,'(a4,i4)') 'ASD ',ICEN
      write(u6,*) 'Looking for ',SDPROP
      do KPROP=1,NPROP
        if ((PNAME(KPROP)(1:3) == SDPROP(1:3)) .and. (PNAME(KPROP)(5:8) == SDPROP(5:8))) then
          if (ICOMP(KPROP) == 1) IAMFI1 = KPROP
          if (ICOMP(KPROP) == 2) IAMFI2 = KPROP
          if (ICOMP(KPROP) == 3) IAMFI3 = KPROP
          if (ICOMP(KPROP) == 4) IAMFI4 = KPROP
          if (ICOMP(KPROP) == 5) IAMFI5 = KPROP
          if (ICOMP(KPROP) == 6) IAMFI6 = KPROP
        end if
      end do

      call Allocate_Z()

      do ISS=1,NSS
        ISTATE = MAPST(ISS)
        MPLET1 = MAPSP(ISS)
        MSPROJ1 = MAPMS(ISS)
        S1 = Half*real(MPLET1-1,kind=wp)
        SM1 = Half*real(MSPROJ1,kind=wp)
        do JSS=1,NSS
          JSTATE = MAPST(JSS)
          MPLET2 = MAPSP(JSS)
          MSPROJ2 = MAPMS(JSS)
          S2 = Half*real(MPLET2-1,kind=wp)
          SM2 = Half*real(MSPROJ2,kind=wp)
          AMFI1 = Zero
          AMFI2 = Zero
          AMFI3 = Zero
          AMFI4 = Zero
          AMFI5 = Zero
          AMFI6 = Zero

          ! SD
          if (IAMFI1 /= 0) AMFI1 = PROP(ISTATE,JSTATE,IAMFI1)
          if (IAMFI2 /= 0) AMFI2 = PROP(ISTATE,JSTATE,IAMFI2)
          if (IAMFI3 /= 0) AMFI3 = PROP(ISTATE,JSTATE,IAMFI3)
          if (IAMFI4 /= 0) AMFI4 = PROP(ISTATE,JSTATE,IAMFI4)
          if (IAMFI5 /= 0) AMFI5 = PROP(ISTATE,JSTATE,IAMFI5)
          if (IAMFI6 /= 0) AMFI6 = PROP(ISTATE,JSTATE,IAMFI6)

          ACNT = Zero

          !if (IFACALFC) then
          ACNT = AMFI6/OneHalf
          !else if ((ISS == 1) .and. (JSS == 1)) then
          !  write(u6,*) '********************'
          !  write(u6,*) '* Skipping FC Part *'
          !  write(u6,*) '********************'
          !end if

          AMFI6 = -AMFI1-AMFI4

          !if (IFACALSD) then
          AMFI1 = -AMFI1
          AMFI2 = -AMFI2
          AMFI3 = -AMFI3
          AMFI4 = -AMFI4
          AMFI5 = -AMFI5
          AMFI6 = -AMFI6
          !else
          !  if ((ISS == 1) .and. (JSS == 1)) then
          !    write(u6,*) '********************'
          !    write(u6,*) '* Skipping SD Part *'
          !    write(u6,*) '********************'
          !  end if
          !  AMFI1 = Zero
          !  AMFI2 = Zero
          !  AMFI3 = Zero
          !  AMFI4 = Zero
          !  AMFI5 = Zero
          !  AMFI6 = Zero
          !end if

          AMFI1 = AMFI1+ACNT
          AMFI4 = AMFI4+ACNT
          AMFI6 = AMFI6+ACNT

          ! WIGNER-ECKART THEOREM:
          FACT = One/sqrt(real(MPLET1,kind=wp))
          if (MPLET1 == MPLET2-2) FACT = -FACT
          CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
          CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
          CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
          CGX = sqrt(Half)*(CGM-CGP)
          CGY = sqrt(Half)*(CGM+CGP)

          ZXYZR(ISS,JSS,1) = CGX*AMFI1+CG0*AMFI3
          ZXYZI(ISS,JSS,1) = CGY*AMFI2
          ZXYZR(ISS,JSS,2) = CGX*AMFI2+CG0*AMFI5
          ZXYZI(ISS,JSS,2) = CGY*AMFI4
          ZXYZR(ISS,JSS,3) = CGX*AMFI3+CG0*AMFI6
          ZXYZI(ISS,JSS,3) = CGY*AMFI5

        end do
      end do

      GSENERGY = ENERGY(1)
      do ISTATE=2,NSTATE
        if (ENERGY(ISTATE) < GSENERGY) GSENERGY = ENERGY(ISTATE)
      end do

      IMLTPL = 0
      do ISTATE=1,NSTATE
        ISGS(IMLTPL+1:IMLTPL+MLTPLT(JBNUM(ISTATE))) = abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp
        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do

      ZEKL(:,:,:,:) = cZero

      IMLTPL = 0
      do ISTATE=1,NSTATE

        ISTART = IMLTPL+1
        IFINAL = IMLTPL+MLTPLT(JBNUM(ISTATE))

        if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

          ! Contribution of the GS spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                  !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
                end if
              end do
            end do
          end do

        else

          ! Contributions of the ES spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                end if
                !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
              end do
            end do
          end do

        end if

        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do
      GCONT(:,:,:) = cZero
      GTOTAL(:,:) = Zero

      do ISTATE=1,NSTATE
        do IXYZ=1,3
          do JXYZ=1,3

            if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

              ! Contributions for the GS's
              do JSTATE=1,NSTATE
                if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) then
                  do I=1,2
                    do J=1,2
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    end do
                  end do
                end if
              end do

            else
              ! Contributions for the ES's
              do JSTATE=1,NSTATE
                do I=1,2
                  do J=1,2
                    GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                       ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) &
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                  end do
                end do
              end do

            end if

          end do
        end do

        GTOTAL(:,:) = GTOTAL(:,:)+real(GCONT(:,:,ISTATE))
      end do

      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,1),ZXYZI(:,:,1))
      call MULMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1),eex,Z)
      DIPSOfcsd(1,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,2),ZXYZI(:,:,2))
      call MULMAT(NSS,ZXYZR(:,:,2),ZXYZI(:,:,2),eey,Z)
      DIPSOfcsd(2,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,3),ZXYZI(:,:,3))
      call MULMAT(NSS,ZXYZR(:,:,3),ZXYZI(:,:,3),eez,Z)
      DIPSOfcsd(3,:,:) = Z(:,:)

      ISS = 1
      do while ((ENSOR(min(ISS,NSS))-ENSOR(1) <= EPRATHR) .and. (ISS <= NSS))

        GTENS(:,:) = Zero

        KDGN = 1
        do JSS=ISS+1,NSS
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (IFATCALSA .and. IFGTSHSA) then
            KDGN = MULTIP
            !write(u6,*) 'KDGN=',KDGN
          else if (abs(EDIFF) < 1.0e-6_wp) then
            KDGN = KDGN+1
          end if
        end do

        write(u6,*)
        do I=1,KDGN
          write(u6,'(3x,A9,I4,3x,A4,F18.8)') 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
        end do
        write(u6,'(3x,A46)') '----------------------------------------------'
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (IFATCALSA) then
          if (ISS == 1) IFUNCT = 0
          call SINANI(KDGN,IFUNCT,NSS,DIPSOfcsd,SPNSFS)
          IFUNCT = IFUNCT+KDGN
        end if
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (KDGN /= 2) then
          write(u6,*) 'no twofold degeneracy'
        else

          if ((ISS == 1) .and. (KDGN == 2) .and. (IPGLOB >= 3)) then
            write(u6,*) 'Experimental: SFS contributions to G=gg+'
            write(u6,*)
            write(u6,'(a6,9(5x,a2,5x))') 'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
            write(u6,*)
            do ISTATE=1,NSTATE
              write(u6,'(2x,I2,2x,9(F12.6))') ISTATE,real(GCONT(:,:,ISTATE))
            end do

            write(u6,*)
            write(u6,'(A6,9(F12.6))') 'total ',GTOTAL(:,:)
          end if
          JSS = ISS+1

          do IXYZ=1,3
            do JXYZ=1,3
              GTIJ = Zero
              do ISO=ISS,JSS
                GTIJ = GTIJ+sum(ZXYZR(ISO,ISS:JSS,IXYZ)*ZXYZR(ISS:JSS,ISO,JXYZ)-ZXYZI(ISO,ISS:JSS,IXYZ)*ZXYZI(ISS:JSS,ISO,JXYZ))
              end do
              GTENS(IXYZ,JXYZ) = Two*GTIJ
            end do
          end do

          !if (IPGLOB > 3) then
          write(u6,*) 'G tensor = gg+'
          write(u6,*)
          write(u6,'(6x,3(6x,a2,4x))') (xyzchr(IXYZ),IXYZ=1,3)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f18.8,1x))') xyzchr(IXYZ),GTENS(IXYZ,:)
          end do
          !end if

          EVR(:) = Zero
          EVI(:) = Zero
          TMPMAT(:,:) = GTENS(:,:)
          call unitmat(TMPVEC,3)

          call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

          ! construct g_s matrix from G by back-transormation of the
          ! square root of the G eigenvalues
          do IXYZ=1,3
            do JXYZ=1,3
              GTENS(IXYZ,JXYZ) = sum(TMPVEC(IXYZ,:)*sqrt(EVR(:))*TMPVEC(JXYZ,:))
            end do
          end do

          write(u6,'(6x,3(5x,a2,5x),4x,4x,2x,8x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('a_',IXYZ,IXYZ=1,3)
          write(u6,*)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f10.6,1x),4x,a2,i1,a1,2x,f12.9,2x,a2,2x,3(1x,f8.4,1x))') xyzchr(IXYZ),GTENS(IXYZ,:),'a_', &
                                                                                              IXYZ,':',sqrt(EVR(IXYZ)), &
                                                                                              xyzchr(IXYZ),TMPVEC(IXYZ,:)
          end do

        end if

        ISS = ISS+KDGN

      end do

      call Deallocate_Z()

    end if

    ! Skip if not a hyperfine calculation
    if (IFACALPSO) then

      write(u6,*)
      write(u6,*) '  ========================================='
      write(u6,*) '  A (PSO)-Matrix for center:',ICEN
      write(u6,*) '  ========================================='

      write(PSOPROP,'(a4,i4)') 'PSOP',ICEN
      write(u6,*) 'Looking for ',PSOPROP

      call Allocate_and_Load_PSOP()

      call Allocate_Z()

      ZXYZI(:,:,1) = ZXYZI(:,:,1)+LXI(:,:)
      ZXYZI(:,:,2) = ZXYZI(:,:,2)+LYI(:,:)
      ZXYZI(:,:,3) = ZXYZI(:,:,3)+LZI(:,:)

      call Deallocate_PSOP()

      GSENERGY = ENERGY(1)
      do ISTATE=2,NSTATE
        if (ENERGY(ISTATE) < GSENERGY) GSENERGY = ENERGY(ISTATE)
      end do

      IMLTPL = 0
      do ISTATE=1,NSTATE
        ISGS(IMLTPL+1:IMLTPL+MLTPLT(JBNUM(ISTATE))) = abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp
        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do
      !*****************************************************
      ZEKL(:,:,:,:) = cZero

      IMLTPL = 0
      do ISTATE=1,NSTATE

        ISTART = IMLTPL+1
        IFINAL = IMLTPL+MLTPLT(JBNUM(ISTATE))

        if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

          ! Contribution of the GS spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                  !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
                end if
              end do
            end do
          end do

        else

          ! Contributions of the ES spin components
          do IXYZ=1,3
            do ISS=ISTART,IFINAL
              do JSS=1,NSS
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                if (ISGS(JSS)) then
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
                  call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
                end if
                !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
              end do
            end do
          end do

        end if

        !do IXYZ=1,3
        !  write(u6,720) 'ZEKL',IXYZ,ISTATE,ZEKL(:,:,IXYZ,ISTATE)
        !end do

        IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
      end do

      GCONT(:,:,:) = cZero
      GTOTAL(:,:) = Zero

      do ISTATE=1,NSTATE
        do IXYZ=1,3
          do JXYZ=1,3

            if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

              ! Contributions for the GS's
              do JSTATE=1,NSTATE
                if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) then
                  do I=1,2
                    do J=1,2
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    end do
                  end do
                end if
              end do

            else

              ! Contributions for the ES's
              do JSTATE=1,NSTATE
                do I=1,2
                  do J=1,2
                    GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                       ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                    if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) &
                      GCONT(JXYZ,IXYZ,ISTATE) = GCONT(JXYZ,IXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE)+ &
                                                                         ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
                  end do
                end do
              end do

            end if

          end do
        end do

        GTOTAL(:,:) = GTOTAL(:,:)+real(GCONT(:,:,ISTATE))
      end do

      ! Continue original calculation of G tensor (=gg^*)

      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,1),ZXYZI(:,:,1))
      call MULMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1),eex,Z)
      DIPSOfpso(1,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,2),ZXYZI(:,:,2))
      call MULMAT(NSS,ZXYZR(:,:,2),ZXYZI(:,:,2),eey,Z)
      DIPSOfpso(2,:,:) = Z(:,:)
      call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,3),ZXYZI(:,:,3))
      call MULMAT(NSS,ZXYZR(:,:,3),ZXYZI(:,:,3),eez,Z)
      DIPSOfpso(3,:,:) = Z(:,:)

      ISS = 1
      do while ((ENSOR(min(ISS,NSS))-ENSOR(1) <= EPRATHR) .and. (ISS <= NSS))

        GTENS(:,:) = Zero

        KDGN = 1
        do JSS=ISS+1,NSS
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (IFATCALSA .and. IFGTSHSA) then
            KDGN = MULTIP
            !write(u6,*) 'KDGN=',KDGN
          else if (abs(EDIFF) < 1.0e-6_wp) then
            KDGN = KDGN+1
          end if
        end do

        write(u6,*)
        do I=1,KDGN
          write(u6,'(3x,A9,I4,3x,A4,F18.8)') 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
        end do
        write(u6,'(3x,A46)') '----------------------------------------------'
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (IFATCALSA) then
          if (ISS == 1) IFUNCT = 0
          call SINANI(KDGN,IFUNCT,NSS,DIPSOfpso,SPNSFS)
          IFUNCT = IFUNCT+KDGN
        end if
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (KDGN /= 2) then
          write(u6,*) 'no twofold degeneracy'
        else

          if ((ISS == 1) .and. (KDGN == 2) .and. (IPGLOB >= 3)) then
            write(u6,*) 'Experimental: SFS contributions to G=gg+'
            write(u6,*)
            write(u6,'(a6,9(5x,a2,5x))') 'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
            write(u6,*)
            do ISTATE=1,NSTATE
              write(u6,'(2x,I2,2x,9(F12.6))') ISTATE,real(GCONT(:,:,ISTATE))
            end do

            write(u6,*)
            write(u6,'(A6,9(F12.6))') 'total ',GTOTAL(:,:)
          end if

          JSS = ISS+1

          do IXYZ=1,3
            do JXYZ=1,3
              GTIJ = Zero
              do ISO=ISS,JSS
                GTIJ = GTIJ+sum(ZXYZR(ISO,ISS:JSS,IXYZ)*ZXYZR(ISS:JSS,ISO,JXYZ)-ZXYZI(ISO,ISS:JSS,IXYZ)*ZXYZI(ISS:JSS,ISO,JXYZ))
              end do
              GTENS(IXYZ,JXYZ) = Two*GTIJ
            end do
          end do

          !if (IPGLOB > 3) then
          write(u6,*) 'G tensor = gg+'
          write(u6,*)
          write(u6,'(6x,3(6x,a2,4x))') (xyzchr(IXYZ),IXYZ=1,3)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f18.8,1x))') xyzchr(IXYZ),GTENS(IXYZ,:)
          end do
          !end if

          EVR(:) = Zero
          EVI(:) = Zero
          TMPMAT(:,:) = GTENS(:,:)
          call unitmat(TMPVEC,3)

          call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

          ! construct g_s matrix from G by back-transormation of the
          ! square root of the G eigenvalues

          write(u6,*)

          ICOUNT = 0
          do IXYZ=1,3
            do JXYZ=1,3
              GTENS(IXYZ,JXYZ) = Zero
              do KXYZ=1,3
                if (EVR(KXYZ) < THRSH) then ! EVR(KXYZ) should probably be ZERO
                  if (EVR(KXYZ) < -THRSH) then
                    write(u6,*) 'WARNNING : negative real eigenvalue of G'
                    write(u6,*) 'WARNNING : this may cause numerical problem'
                    write(u6,*) 'WARNNING : and may give "NAN" in G tensor'
                  else
                    EVR(KXYZ) = ZERO
                    if (ICOUNT == 0) then
                      write(u6,*) ' The eigenvalue "EVR(KXYZ)" will be set to ZERO'
                      write(u6,*) ' each time it is very small'
                      write(u6,*) ' to avoid the roots of negative eigenvalues'
                    end if
                    ICOUNT = 1
                  end if
                end if
              end do
              GTENS(IXYZ,JXYZ) = sum(TMPVEC(IXYZ,:)*sqrt(EVR(:))*TMPVEC(JXYZ,:))
            end do
          end do
          write(u6,*)
          write(u6,'(6x,3(5x,a2,5x),4x,4x,2x,8x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('a_',IXYZ,IXYZ=1,3)
          write(u6,*)
          do IXYZ=1,3
            write(u6,'(2x,a2,2x,3(1x,f10.6,1x),4x,a2,i1,a1,2x,f12.9,2x,a2,2x,3(1x,f8.4,1x))') xyzchr(IXYZ),GTENS(IXYZ,:),'a_', &
                                                                                              IXYZ,':',sqrt(EVR(IXYZ)), &
                                                                                              xyzchr(IXYZ),TMPVEC(IXYZ,:)
          end do

          call ADD_INFO('ATENS_PSO',GTENS,9,5)
          call ADD_INFO('ATENS_PSOEVR',EVR,3,5)

          write(u6,*)
        end if

        ISS = ISS+KDGN

      end do

      call Deallocate_Z()

    end if

  end if
  ! End loop over CNT properties
end do

call mma_deallocate(MAPST)
call mma_deallocate(MAPSP)
call mma_deallocate(MAPMS)
call mma_deallocate(ISGS)
call mma_deallocate(PNMRT)
call mma_deallocate(PNMR)
call mma_deallocate(PNMRC)
call mma_deallocate(PNMRD)
call mma_deallocate(TMPf)
call mma_deallocate(DIPSOf)
call mma_deallocate(DIPSOfc)
call mma_deallocate(DIPSOfcsd)
call mma_deallocate(DIPSOfsd)
call mma_deallocate(DIPSOfpso)
call mma_deallocate(DIMSO)
call mma_deallocate(SPNSFS)
call mma_deallocate(Z)
call mma_deallocate(ZEKL)
call mma_deallocate(GCONT)

!710  FORMAT(A4,4I4,4(2X,'('F12.8','F12.8')'))
!720  FORMAT(A4,2I4,4(2X,'('F12.8','F12.8')'))

contains

subroutine Allocate_and_Load_PSOP()

  integer KPROP
  integer IAMX, IAMY, IAMZ

  IAMX = 0
  IAMY = 0
  IAMZ = 0
  do KPROP=1,NPROP
    if (PNAME(KPROP) == PSOPROP) then
      if (ICOMP(KPROP) == 1) IAMX = KPROP
      if (ICOMP(KPROP) == 2) IAMY = KPROP
      if (ICOMP(KPROP) == 3) IAMZ = KPROP
    end if
  end do
  call mma_allocate(LXI,NSS,NSS,Label='LXI')
  LXI(:,:) = Zero
  call mma_allocate(LYI,NSS,NSS,Label='LYI')
  LYI(:,:) = Zero
  call mma_allocate(LZI,NSS,NSS,Label='LZI')
  LZI(:,:) = Zero
  if (IAMX > 0) call SMMAT(PROP,LXI,NSS,IAMX,1)
  if (IAMY > 0) call SMMAT(PROP,LYI,NSS,IAMY,2)
  if (IAMZ > 0) call SMMAT(PROP,LZI,NSS,IAMZ,3)

end subroutine Allocate_and_Load_PSOP

subroutine Deallocate_PSOP()

  call mma_deallocate(LXI)
  call mma_deallocate(LYI)
  call mma_deallocate(LZI)

end subroutine Deallocate_PSOP

subroutine Allocate_Z()

  call mma_allocate(ZXYZR,NSS,NSS,3,Label='ZXYZR')
  call mma_allocate(ZXYZI,NSS,NSS,3,Label='ZXYZI')
  ZXYZR(:,:,:) = Zero
  ZXYZI(:,:,:) = Zero

end subroutine Allocate_Z

subroutine Deallocate_Z()

  call mma_deallocate(ZXYZR)
  call mma_deallocate(ZXYZI)

end subroutine Deallocate_Z

end subroutine HFCTS
