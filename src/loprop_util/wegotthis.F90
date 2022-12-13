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

subroutine WeGotThis(nAt,nB,MP,nij,EC,lMax,iPrint,Pot_Expo,Pot_Point,Pot_Fac,Diffed)

use Index_Functions, only: nTri3_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAt, nB, nij, lMax, iPrint
real(kind=wp), intent(in) :: MP(nij,*), EC(3,nij), Pot_Expo(nij*2), Pot_Point(nij), Pot_Fac(nij*4)
logical(kind=iwp), intent(in) :: Diffed(nij*2)
integer(kind=iwp) :: iA, iComp, iDC, iOpt, iPP, irc, iSmLbl, jA, k, kaunt, kauntA, kComp, l, nDens, nEPP, nImprove, nShitty
real(kind=wp) :: chP, CorrCoeff, DeNom, Dif1, Dif2, ElPot_APP, ElPot_MP, ElPot_REF, ErrAv1, ErrAv2, ErrCorr, ErrDe1, ErrDe2, &
                 ErrMax1, ErrMax2, ErrRe1, ErrRe2, ErrVar1, ErrVar2, Expo(4), PImp, PP, PShi, r, rinv, x, y, z
logical(kind=iwp) :: D1, D2, Found, Que
character(len=10) :: OneFile, Label
real(kind=wp), allocatable :: D1ao(:), dMullig(:), ElP(:), EPCo(:,:)
character(len=10), parameter :: DistType(2) = ['Monopole  ','Dipole    ']
real(kind=wp) :: Ddot_, ElPot
interface
  subroutine Diff_Aux1(nEPotPoints,EPCo,nB,OneFile)
    import :: wp, iwp
    integer(kind=iwp), intent(out) :: nEPotPoints
    real(kind=wp), allocatable, intent(out) :: EPCo(:,:)
    integer(kind=iwp), intent(in) :: nB
    character(len=10), intent(in) :: OneFile
  end subroutine Diff_Aux1
end interface

! Print exponents and factors.

write(u6,*)
write(u6,'(A)') ' **********************************************************************************'
write(u6,'(A)') ' *                                                                                *'
write(u6,'(A)') ' *                  Result for the Diffuse Distribution                           *'
write(u6,'(A)') ' *                                                                                *'
write(u6,'(A)') ' **********************************************************************************'
write(u6,*)
write(u6,401) 'Centre','Coordinate','Multipole','Factor','Exponent','Point-charge'
write(u6,402) '|-------------------------------------------------------------------------------------------------|'
kauntA = 0
do iA=1,nAt
  do jA=1,iA
    kauntA = kauntA+1
    do iDC=1,2
      PP = -One
      if (iDC == 1) PP = Pot_Point(kauntA)
      if (Diffed(2*(kauntA-1)+iDC)) then
        if (iDC == 1) then
          write(u6,403) kauntA,(EC(k,kauntA),k=1,3),DistType(iDC),Pot_Fac(4*(kauntA-1)+iDC),Two*Pot_Expo(2*(kauntA-1)+iDC),PP
        else if (iDC == 2) then
          write(u6,404) kauntA,(EC(k,kauntA),k=1,3),DistType(iDC),(Pot_Fac(4*(kauntA-1)+iDC+k),k=0,2),Two*Pot_Expo(2*(kauntA-1)+iDC)
        end if
      else
        if (iDC == 1) then
          write(u6,405) kauntA,(EC(k,kauntA),k=1,3),DistType(iDC),Pot_Fac(4*(kauntA-1)+iDC),'Point(inf)',PP
        else if (iDC == 2) then
          write(u6,406) kauntA,(EC(k,kauntA),k=1,3),DistType(iDC),(Pot_Fac(4*(kauntA-1)+iDC+k),k=0,2),'Point(inf)'
        end if
      end if
    end do
  end do
end do
401 format(' ',A,'     ',A,'             ',A,'             ',A,'          ',A,'     ',A)
402 format(A)
403 format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,'          ',F7.3,'          ',F7.3,'      ',F7.3)
404 format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,' (',2(F7.3,','),F7.3,') ',F7.3)
405 format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,'          ',F7.3,'            ',A,'      ',F7.3)
406 format(' ',I3,' (',2(F7.3,','),F7.3,')      ',A,' (',2(F7.3,','),F7.3,')   ',A)

! If the extra ONEINT-file exist, then do an error analysis.

Que = .false.
call F_Inquire('ONEINTP',Que)
if (Que) then
  write(u6,*)
  write(u6,'(A)') ' Found Test-grid for error analysis.'
  write(u6,*)
  ErrAv1 = Zero
  ErrAv2 = Zero
  ErrDe1 = Zero
  ErrDe2 = Zero
  ErrVar1 = Zero
  ErrVar2 = Zero
  ErrMax1 = Zero
  ErrMax2 = Zero
  ErrCorr = Zero
  nImprove = 0
  nShitty = 0
  DeNom = Zero
  write(OneFile,'(A)') 'ONEINTP'
  call Diff_Aux1(nEPP,EPCo,nB,OneFile)
  call Qpg_dArray('D1ao',Found,nDens)
  if (Found .and. (nDens /= 0)) then
    call mma_allocate(D1ao,nDens,Label='D1ao')
  else
    write(u6,*) 'WeGotThis: do not think so!'
    call Abend()
  end if
  call Get_dArray_chk('D1ao',D1ao,nDens)
  call mma_allocate(ElP,nDens+4,label='ElPot')
  if (iPrint >= 2) then
    write(u6,*)
    write(u6,'(A)') ' Electric Potential'
    write(u6,'(A)') '  Reference   Approximate MP-expanded'
  end if

  ! Loop over all points where the electric potential has been sampled.

  call mma_allocate(dMullig,nTri3_Elem1(lMax),label='dMullig')

  do iPP=1,nEPP

    ! First, get the true electric potential, the reference.

    write(Label,'(A3,I5)') 'EF0',iPP
    irc = -1
    iOpt = 0
    iSmLbl = 0
    iComp = 1
    call RdOne(irc,iOpt,Label,iComp,ElP,iSmLbl)
    ElPot_REF = ElP(nDens+4)
    ElPot_REF = ElPot_REF-Ddot_(nDens,D1ao,1,ElP,1)

    ! Second, get the approximate electric potential and also the
    ! completely multipole expanded potential.

    ElPot_APP = Zero
    ElPot_MP = Zero
    kauntA = 1
    !rMin = 1.0e10_wp
    do iA=1,nAt
      do jA=1,iA
        x = EPCo(1,iPP)-EC(1,kauntA)
        y = EPCo(2,iPP)-EC(2,kauntA)
        z = EPCo(3,iPP)-EC(3,kauntA)
        r = sqrt(x**2+y**2+z**2)
        rinv = One/r
        D1 = Diffed(2*(kauntA-1)+1)
        D2 = Diffed(2*(kauntA-1)+2)
        Expo(1) = Pot_Expo(2*(kauntA-1)+1)
        Expo(2) = Pot_Expo(2*(kauntA-1)+2)
        chP = Pot_Point(kauntA)
        kaunt = 1
        !rmin = min(r,rmin)
        do l=0,lMax
          kComp = (l+1)*(l+2)/2
          do k=1,kComp
            dMullig(kaunt) = MP(kauntA,kaunt)
            kaunt = kaunt+1
          end do
        end do
        ElPot_APP = ElPot_APP+ElPot(r,rinv,x,y,z,dMullig,lMax,Expo,chP,D1,D2)
        ElPot_MP = ElPot_MP+ElPot(r,rinv,x,y,z,dMullig,lMax,Expo,chP,.false.,.false.)
        kauntA = kauntA+1
      end do
    end do
    !write(u6,*) 'Minimum Dist:',rMin

    ! Print if requested.

    if (iPrint >= 2) then
      write(u6,441) ElPot_REF,ElPot_APP,ElPot_MP,EPCo(:,iPP)
    end if

    ! Third, accumulate to error analysis.

    ! The difference
    Dif1 = ElPot_APP-ElPot_REF
    Dif2 = ElPot_MP-ElPot_REF
    ! Accumulate to average error
    ErrAv1 = ErrAv1+Dif1
    ErrAv2 = ErrAv2+Dif2
    ! Accumulate to variance of error
    ErrVar1 = ErrVar1+Dif1**2
    ErrVar2 = ErrVar2+Dif2**2
    ! Accumulate to deviation
    ErrDe1 = ErrDe1+abs(Dif1)
    ErrDe2 = ErrDe2+abs(Dif2)
    ! Maximal error
    if (ErrMax1 < abs(Dif1)) ErrMax1 = abs(Dif1)
    if (ErrMax2 < abs(Dif2)) ErrMax2 = abs(Dif2)
    ! Accumulate to covariance of errors
    ErrCorr = ErrCorr+Dif1*Dif2
    ! Better or worse
    if (abs(Dif1) <= abs(Dif2)) nImprove = nImprove+1
    if (abs(Dif2) < abs(Dif1)) nShitty = nShitty+1
    ! Accumulate to denominator in relative error
    DeNom = DeNom+abs(ElPot_REF)
  end do
  ErrRe1 = ErrDe1/DeNom
  ErrRe2 = ErrDe2/DeNom
  ErrAv1 = ErrAv1/real(nEPP,kind=wp)
  ErrAv2 = ErrAv2/real(nEPP,kind=wp)
  ErrDe1 = ErrDe1/real(nEPP,kind=wp)
  ErrDe2 = ErrDe2/real(nEPP,kind=wp)
  ErrVar1 = ErrVar1/real(nEPP,kind=wp)
  ErrVar2 = ErrVar2/real(nEPP,kind=wp)
  ErrCorr = ErrCorr/real(nEPP,kind=wp)
  CorrCoeff = (ErrCorr-ErrAv1*ErrAv2)/sqrt((ErrVar1-ErrAv1**2)*(ErrVar2-ErrAv2**2))
  PImp = 100.0_wp*real(nImprove,kind=wp)/real(nEPP,kind=wp)
  PShi = 100.0_wp*real(nShitty,kind=wp)/real(nEPP,kind=wp)

  ! Four, print the analysis.

  write(u6,*)
  write(u6,'(A)') '   Error Analysis'
  write(u6,'(A)') ' |----------------------------------------------------|'
  write(u6,'(A)') '                                Diffuse     MP-expanded'
  write(u6,442) '  Average absolute error:      ',ErrAv1,ErrAv2
  write(u6,442) '  Average absolute deviation:  ',ErrDe1,ErrDe2
  write(u6,442) '  Average relative error:      ',ErrRe1,ErrRe2
  write(u6,442) '  Maximal deviation:           ',ErrMax1,ErrMax2
  write(u6,*)
  write(u6,443) '  Error correlation:           ',CorrCoeff
  write(u6,444) '  Smaller error than MP:       ',PImp,'%'
  write(u6,444) '  Larger error than MP:        ',PShi,'%'
  write(u6,'(A)') ' |----------------------------------------------------|'
  write(u6,*)

  ! Deallocate

  call mma_deallocate(ElP)
  call mma_deallocate(D1ao)
  call mma_deallocate(EPCo)
  call mma_deallocate(dMullig)
  irc = -1
  call ClsOne(irc,0)
end if

441 format(3(F12.8),'    In: ',3(F8.4))
442 format(A,2(F12.8))
443 format(A,F12.8)
444 format(A,F12.8,A)

return

end subroutine WeGotThis
