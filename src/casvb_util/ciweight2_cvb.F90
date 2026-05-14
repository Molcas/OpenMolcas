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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine ciweight2_cvb(civec,civbs,civb,citmp,civec5,orbs,sorbs,orbinv,owrk,ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)

use casvb_global, only: form2AD, formAD, gjorb, gjorb2, gjorb3, iciweights, iprec, nalf, nbet, nda, ndet, nel, norb, npcf, svb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: civec(0:ndet), civbs(0:ndet), civb(0:ndet), citmp(0:ndet), civec5(0:ndet)
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), intent(out) :: sorbs(norb,norb), orbinv(norb,norb), owrk(norb,norb)
integer(kind=iwp), intent(in) :: ionmin, ionmax, mxrem, mxsng, mxasg, ncnfcas, mxdetcas
integer(kind=iwp) :: i, ia, iaorb, ib, ibeg, ibegt, iborb, ic, idet, ilin, inda, indab, indasg, indb, indion, indsng, ion, iorb, &
                     ix1, lenfld, mp, mrem, nalfsng, nbetsng, nc, nindasg, nprint, nsing, rc
real(kind=wp) :: c1, c2, c3, c4, cnrm, cprint(6), fac, fac1, fac2, s11, s12, s22, sm1, sm2, sum1, sum2, total1, total2, total3, &
                 total4, total5, total6
character(len=240) :: line
integer(kind=iwp), allocatable :: iaocc(:), ibocc(:), indavec(:), indbvec(:), locasg(:), locion(:), locsng(:), lunasg(:), &
                                  lunion(:), lunsng(:), maxgasg(:), maxgion(:), maxgrph(:), maxgsng(:), mingasg(:), mingion(:), &
                                  mingrph(:), mingsng(:), nkasg(:), nkion(:), nksng(:), xalf(:,:), xasg(:), xbet(:,:), xion(:), &
                                  xsng(:)
real(kind=wp), allocatable :: gal1(:), gal2(:), wghtion1(:), wghtion2(:), wghtion3(:), wghtion4(:), wghtion5(:), wghtion6(:)
integer(kind=iwp), external :: indget_cvb

call cidot_cvb(civb,civbs,cnrm)
fac = svb/sqrt(cnrm)

call cicopy_cvb(civec,citmp)
orbinv(:,:) = orbs(:,:)
call mxinv_cvb(orbinv,norb)
call gaussj_cvb(orbinv,gjorb)
call applyt_cvb(civec,gjorb)
! Chirgwin-Coulson weights
if (btest(iciweights,0)) then
  call trnsps(norb,norb,orbs,owrk)
  call gaussj_cvb(owrk,gjorb2)
  call applyt_cvb(citmp,gjorb2)
  do idet=1,ndet
    civbs(idet) = (citmp(idet)-fac*civbs(idet))*(civec(idet)-fac*civb(idet))
  end do
  do idet=1,ndet
    citmp(idet) = citmp(idet)*civec(idet)
  end do
end if
do idet=1,ndet
  civb(idet) = civec(idet)-fac*civb(idet)
end do

call mma_allocate(mingrph,[0,norb],label='mingrph')
call mma_allocate(maxgrph,[0,norb],label='maxgrph')
call mma_allocate(xalf,[0,nel],[0,nalf],label='xalf')
call mma_allocate(xbet,[0,nel],[0,nbet],label='xbet')
call mma_allocate(iaocc,norb,label='iaocc')
call mma_allocate(ibocc,norb,label='ibocc')
call mma_allocate(mingion,[0,norb],label='mingion')
call mma_allocate(maxgion,[0,norb],label='maxgion')
call mma_allocate(nkion,[0,norb],label='nkion')
call mma_allocate(xion,(norb+1)*(ionmax+1),label='xion')
call mma_allocate(locion,norb,label='locion')
call mma_allocate(lunion,norb,label='lunion')
call mma_allocate(mingsng,[0,norb],label='mingsng')
call mma_allocate(maxgsng,[0,norb],label='maxgsng')
call mma_allocate(nksng,[0,norb],label='nksng')
call mma_allocate(xsng,(mxrem+1)*(mxsng+1),label='xsng')
call mma_allocate(locsng,norb,label='locsng')
call mma_allocate(lunsng,norb,label='lunsng')
call mma_allocate(mingasg,[0,norb],label='mingasg')
call mma_allocate(maxgasg,[0,norb],label='maxgasg')
call mma_allocate(nkasg,[0,norb],label='nkasg')
call mma_allocate(xasg,(mxsng+1)*(mxasg+1),label='xasg')
call mma_allocate(locasg,norb,label='locasg')
call mma_allocate(lunasg,norb,label='lunasg')

! Inverse-overlap weights
if (btest(iciweights,2)) then
  call mma_allocate(gal1,ncnfcas,label='gal1')
  call mma_allocate(gal2,ncnfcas,label='gal2')
  call mma_allocate(indavec,mxdetcas,label='indavec')
  call mma_allocate(indbvec,mxdetcas,label='indbvec')
  call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
  orbinv(:,:) = sorbs(:,:)
  call mxinv_cvb(orbinv,norb)
  call gaussj_cvb(orbinv,gjorb)
  ! Alpha weight array:
  do iorb=0,norb
    mingrph(iorb) = max(iorb-norb+nalf,0)
    maxgrph(iorb) = min(iorb,nalf)
  end do
  call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
  ! Beta weight array:
  do iorb=0,norb
    mingrph(iorb) = max(iorb-norb+nbet,0)
    maxgrph(iorb) = min(iorb,nbet)
  end do
  call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

  nc = 0
  do ion=ionmin,ionmax
    nsing = nel-2*ion
    mrem = norb-ion
    nalfsng = nalf-ion
    nbetsng = nbet-ion
    ! Initialise loop for ionic orbitals
    do iorb=0,norb
      mingion(iorb) = max(iorb-norb+ion,0)
      maxgion(iorb) = min(iorb,ion)
    end do
    call weight_cvb(xion,mingion,maxgion,ion,norb)
    nkion(:) = maxgion(:)
    call occupy_cvb(nkion,norb,locion,lunion)
    ! Initialise loop for singly occupied orbitals
    do iorb=0,mrem
      mingsng(iorb) = max(iorb-mrem+nsing,0)
      maxgsng(iorb) = min(iorb,nsing)
    end do
    call weight_cvb(xsng,mingsng,maxgsng,nsing,mrem)
    nksng(:) = maxgsng(:)
    call occupy_cvb(nksng,mrem,locsng,lunsng)
    ! Initialise loop for singly occupied alpha orbitals
    do iorb=0,nsing
      mingasg(iorb) = max(iorb-nsing+nalfsng,0)
      maxgasg(iorb) = min(iorb,nalfsng)
    end do
    call weight_cvb(xasg,mingasg,maxgasg,nalfsng,nsing)
    nkasg(:) = maxgasg(:)
    call occupy_cvb(nkasg,nsing,locasg,lunasg)

    ! Loop ionic
    indion = 1
    do
      ! Loop singly occupied
      indsng = 1
      do
        civec5(1:) = Zero
        s11 = Zero
        s22 = Zero
        s12 = Zero
        ! Loop singly occupied alpha
        indasg = 1

        do
          iaocc(:) = 0
          ibocc(:) = 0
          do i=1,ion
            iaocc(locion(i)) = 1
            ibocc(locion(i)) = 1
          end do
          do ia=1,nalfsng
            iaorb = lunion(locsng(locasg(ia)))
            iaocc(iaorb) = 1
          end do
          do ib=1,nbetsng
            iborb = lunion(locsng(lunasg(ib)))
            ibocc(iborb) = 1
          end do
          inda = indget_cvb(iaocc,nalf,norb,xalf)
          indb = indget_cvb(ibocc,nbet,norb,xbet)
          indavec(indasg) = inda
          indbvec(indasg) = indb

          indab = (indb-1)*nda+inda
          civec5(indab) = civec(indab)
          s11 = s11+civec(indab)*civec(indab)
          s22 = s22+civb(indab)*civb(indab)
          s12 = s12+civec(indab)*civb(indab)

          call loind_cvb(nsing,nalfsng,nkasg,mingasg,maxgasg,locasg,lunasg,indasg,xasg,rc)
          if (rc == 0) exit
        end do

        call applyt_cvb(civec5,gjorb)

        call icomb_cvb(nsing,nalfsng,nindasg)
        sm1 = Zero
        do indasg=1,nindasg
          inda = indavec(indasg)
          indb = indbvec(indasg)
          indab = (indb-1)*nda+inda
          sm1 = sm1+civec(indab)*civec5(indab)
        end do

        if (abs(sm1) > 1.0e-20_wp) then
          sm1 = s11*s11/sm1
        else if (abs(sm1) <= 1.0e-20_wp) then
          sm1 = Zero
        end if

        civec5(1:) = Zero
        do indasg=1,nindasg
          inda = indavec(indasg)
          indb = indbvec(indasg)
          indab = (indb-1)*nda+inda
          civec5(indab) = civb(indab)
        end do

        call applyt_cvb(civec5,gjorb)

        sm2 = Zero
        do indasg=1,nindasg
          inda = indavec(indasg)
          indb = indbvec(indasg)
          indab = (indb-1)*nda+inda
          sm2 = sm2+civb(indab)*civec5(indab)
        end do

        if (abs(sm2) > 1.0e-20_wp) then
          sm2 = s22*s22/sm2
        else if (abs(sm2) <= 1.0e-20_wp) then
          sm2 = Zero
        end if

        nc = nc+1
        gal1(nc) = sm1
        gal2(nc) = sm2

        call loind_cvb(mrem,nsing,nksng,mingsng,maxgsng,locsng,lunsng,indsng,xsng,rc)
        if (rc == 0) exit
      end do
      call loind_cvb(norb,ion,nkion,mingion,maxgion,locion,lunion,indion,xion,rc)
      if (rc == 0) exit
    end do
  end do
  sum1 = Zero
  sum2 = Zero
  do ic=1,ncnfcas
    sum1 = sum1+gal1(ic)
    sum2 = sum2+gal2(ic)
  end do
  fac1 = One/sum1
  if ((abs(One-svb*svb) < 1.0e-20_wp) .and. (abs(sum2) < 1.0e-20_wp)) then
    fac2 = One
  else
    fac2 = (One-svb*svb)/sum2
  end if
  do ic=1,ncnfcas
    gal1(ic) = fac1*gal1(ic)
    gal2(ic) = fac2*gal2(ic)
  end do
end if

! Weights of Lowdin orthonormalized structures
if (btest(iciweights,1)) then
  call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
  call mxsqrt_cvb(sorbs,norb,1)
  call gaussj_cvb(sorbs,gjorb3)
  call applyt_cvb(civec,gjorb3)
  call applyt_cvb(civb,gjorb3)
  call cidot_cvb(civb,civb,cnrm)
  fac = svb/sqrt(cnrm)
  civb(1:) = civb(1:)**2
  civec(1:) = civec(1:)**2
end if

write(u6,'(/,2a)') ' Weights of CASSCF configurations in VB basis (c_res=c_cas-Svb*c_vb) :'
write(u6,'(2a)') ' ---------------------------------------------------------------------'
if (btest(iciweights,2)) then
  write(u6,'(a)') ' Sum of inverse-overlap weights :'
  write(u6,form2AD) ' c_cas :',sum1,' expected :',One
  write(u6,form2AD) ' c_res :',sum2,' expected :',One-svb*svb
  write(u6,'(a)') ' '
end if
lenfld = 8+iprec
ix1 = max(0,min(3,2*lenfld-16))
line = ''
if ((npcf > 0) .or. (npcf == -1)) then
  line(1:21) = '  Conf. =>  Orbitals '
  ibeg = max(14+3*nel,22)
  ibegt = ibeg
  if (btest(iciweights,0)) then
    line(ibegt+ix1:ibegt+2*lenfld-1) = 'Chirgwin-Coulson'
    ibegt = ibegt+2*lenfld
  end if
  if (btest(iciweights,1)) then
    line(ibegt+ix1:ibegt+2*lenfld-1) = 'Lowdin'
    ibegt = ibegt+2*lenfld
  end if
  if (btest(iciweights,2)) then
    line(ibegt+ix1:ibegt+2*lenfld-1) = 'Inverse'
    ibegt = ibegt+2*lenfld
  end if
  write(u6,'(a)') trim(line)
  line = ''
  if (btest(iciweights,0)) then
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
    ibeg = ibeg+lenfld
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
    ibeg = ibeg+lenfld
  end if
  if (btest(iciweights,1)) then
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
    ibeg = ibeg+lenfld
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
    ibeg = ibeg+lenfld
  end if
  if (btest(iciweights,2)) then
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
    ibeg = ibeg+lenfld
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
    ibeg = ibeg+lenfld
  end if
  write(u6,'(a)') trim(line)
end if
line = ''

! Alpha weight array:
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nalf,0)
  maxgrph(iorb) = min(iorb,nalf)
end do
call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
! Beta weight array:
do iorb=0,norb
  mingrph(iorb) = max(iorb-norb+nbet,0)
  maxgrph(iorb) = min(iorb,nbet)
end do
call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)

call mma_deallocate(mingrph)
call mma_deallocate(maxgrph)

nc = 0
call mma_allocate(wghtion1,[ionmin,ionmax],label='wghtion1')
call mma_allocate(wghtion2,[ionmin,ionmax],label='wghtion2')
call mma_allocate(wghtion3,[ionmin,ionmax],label='wghtion3')
call mma_allocate(wghtion4,[ionmin,ionmax],label='wghtion4')
call mma_allocate(wghtion5,[ionmin,ionmax],label='wghtion5')
call mma_allocate(wghtion6,[ionmin,ionmax],label='wghtion6')
wghtion1(:) = Zero
wghtion2(:) = Zero
wghtion3(:) = Zero
wghtion4(:) = Zero
wghtion5(:) = Zero
wghtion6(:) = Zero
do ion=ionmin,ionmax
  nsing = nel-2*ion
  mrem = norb-ion
  nalfsng = nalf-ion
  nbetsng = nbet-ion
  ! Initialise loop for ionic orbitals
  do iorb=0,norb
    mingion(iorb) = max(iorb-norb+ion,0)
    maxgion(iorb) = min(iorb,ion)
  end do
  call weight_cvb(xion,mingion,maxgion,ion,norb)
  nkion(:) = maxgion(:)
  call occupy_cvb(nkion,norb,locion,lunion)
  ! Initialise loop for singly occupied orbitals
  do iorb=0,mrem
    mingsng(iorb) = max(iorb-mrem+nsing,0)
    maxgsng(iorb) = min(iorb,nsing)
  end do
  call weight_cvb(xsng,mingsng,maxgsng,nsing,mrem)
  nksng(:) = maxgsng(:)
  call occupy_cvb(nksng,mrem,locsng,lunsng)
  ! Initialise loop for singly occupied alpha orbitals
  do iorb=0,nsing
    mingasg(iorb) = max(iorb-nsing+nalfsng,0)
    maxgasg(iorb) = min(iorb,nalfsng)
  end do
  call weight_cvb(xasg,mingasg,maxgasg,nalfsng,nsing)
  nkasg(:) = maxgasg(:)
  call occupy_cvb(nkasg,nsing,locasg,lunasg)

  ! Loop ionic
  indion = 1
  do
    ! Loop singly occupied
    indsng = 1
    do
      c1 = Zero
      c2 = Zero
      c3 = Zero
      c4 = Zero
      ! Loop singly occupied alpha
      indasg = 1

      do
        iaocc(:) = 0
        ibocc(:) = 0
        do i=1,ion
          iaocc(locion(i)) = 1
          ibocc(locion(i)) = 1
        end do
        do ia=1,nalfsng
          iaorb = lunion(locsng(locasg(ia)))
          iaocc(iaorb) = 1
        end do
        do ib=1,nbetsng
          iborb = lunion(locsng(lunasg(ib)))
          ibocc(iborb) = 1
        end do
        inda = indget_cvb(iaocc,nalf,norb,xalf)
        indb = indget_cvb(ibocc,nbet,norb,xbet)

        indab = (indb-1)*nda+inda
        c1 = c1+citmp(indab)
        c2 = c2+civbs(indab)
        c3 = c3+civec(indab)
        c4 = c4+civb(indab)

        call loind_cvb(nsing,nalfsng,nkasg,mingasg,maxgasg,locasg,lunasg,indasg,xasg,rc)
        if (rc == 0) exit
      end do
      nc = nc+1
      wghtion1(ion) = wghtion1(ion)+c1
      wghtion2(ion) = wghtion2(ion)+c2
      wghtion3(ion) = wghtion3(ion)+c3
      wghtion4(ion) = wghtion4(ion)+c4
      wghtion5(ion) = wghtion5(ion)+gal1(nc)
      wghtion6(ion) = wghtion6(ion)+gal2(nc)
      if ((nc <= npcf) .or. (npcf == -1)) then
        call int2char_cvb(line,nc,7)
        line(9:10) = '=>'
        ilin = 11
        do i=1,ion
          call int2char_cvb(line(ilin:ilin+2),locion(i),3)
          ilin = ilin+3
          call int2char_cvb(line(ilin:ilin+2),locion(i),3)
          ilin = ilin+3
        end do
        do i=1,nsing
          call int2char_cvb(line(ilin:ilin+2),lunion(locsng(i)),3)
          ilin = ilin+3
        end do
        ilin = max(ilin,19)
        nprint = 0
        if (btest(iciweights,0)) then
          cprint(nprint+1) = c1
          cprint(nprint+2) = c2
          nprint = nprint+2
        end if
        if (btest(iciweights,1)) then
          cprint(nprint+1) = c3
          cprint(nprint+2) = c4
          nprint = nprint+2
        end if
        if (btest(iciweights,2)) then
          cprint(nprint+1) = gal1(nc)
          cprint(nprint+2) = gal2(nc)
          nprint = nprint+2
        end if
        write(u6,formAD) line(1:ilin),(cprint(mp),mp=1,nprint)
      end if
      call loind_cvb(mrem,nsing,nksng,mingsng,maxgsng,locsng,lunsng,indsng,xsng,rc)
      if (rc == 0) exit
    end do
    call loind_cvb(norb,ion,nkion,mingion,maxgion,locion,lunion,indion,xion,rc)
    if (rc == 0) exit
  end do
end do

call mma_deallocate(xalf)
call mma_deallocate(xbet)
call mma_deallocate(iaocc)
call mma_deallocate(ibocc)
call mma_deallocate(mingion)
call mma_deallocate(maxgion)
call mma_deallocate(nkion)
call mma_deallocate(xion)
call mma_deallocate(locion)
call mma_deallocate(lunion)
call mma_deallocate(mingsng)
call mma_deallocate(maxgsng)
call mma_deallocate(nksng)
call mma_deallocate(xsng)
call mma_deallocate(locsng)
call mma_deallocate(lunsng)
call mma_deallocate(mingasg)
call mma_deallocate(maxgasg)
call mma_deallocate(nkasg)
call mma_deallocate(xasg)
call mma_deallocate(locasg)
call mma_deallocate(lunasg)
call mma_deallocate(gal1)
call mma_deallocate(gal2)
call mma_deallocate(indavec)
call mma_deallocate(indbvec)

line = ' Accumulated weights:'
ibeg = 22
if (btest(iciweights,0)) then
  line(ibeg+ix1:ibeg+2*lenfld-1) = 'Chirgwin-Coulson'
  ibeg = ibeg+2*lenfld
end if
if (btest(iciweights,1)) then
  line(ibeg+ix1:ibeg+2*lenfld-1) = 'Lowdin'
  ibeg = ibeg+2*lenfld
end if
if (btest(iciweights,2)) then
  line(ibeg+ix1:ibeg+2*lenfld-1) = 'Inverse'
  ibeg = ibeg+2*lenfld
end if
write(u6,'(/,a)') trim(line)
line = ''
ibeg = 22
if (btest(iciweights,0)) then
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
  ibeg = ibeg+lenfld
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
  ibeg = ibeg+lenfld
end if
if (btest(iciweights,1)) then
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
  ibeg = ibeg+lenfld
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
  ibeg = ibeg+lenfld
end if
if (btest(iciweights,2)) then
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
  ibeg = ibeg+lenfld
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
  ibeg = ibeg+lenfld
end if
write(u6,'(a)') trim(line)
line = ''
total1 = Zero
total2 = Zero
total3 = Zero
total4 = Zero
total5 = Zero
total6 = Zero
do ion=ionmin,ionmax
  total1 = total1+wghtion1(ion)
  total2 = total2+wghtion2(ion)
  total3 = total3+wghtion3(ion)
  total4 = total4+wghtion4(ion)
  total5 = total5+wghtion5(ion)
  total6 = total6+wghtion6(ion)
  line(1:19) = ' For ionicity    : '
  call int2char_cvb(line(14:16),ion,3)
  nprint = 0
  if (btest(iciweights,0)) then
    cprint(nprint+1) = wghtion1(ion)
    cprint(nprint+2) = wghtion2(ion)
    nprint = nprint+2
  end if
  if (btest(iciweights,1)) then
    cprint(nprint+1) = wghtion3(ion)
    cprint(nprint+2) = wghtion4(ion)
    nprint = nprint+2
  end if
  if (btest(iciweights,2)) then
    cprint(nprint+1) = wghtion5(ion)
    cprint(nprint+2) = wghtion6(ion)
    nprint = nprint+2
  end if
  write(u6,formAD) line(1:19),(cprint(mp),mp=1,nprint)
end do
nprint = 0
if (btest(iciweights,0)) then
  cprint(nprint+1) = total1
  cprint(nprint+2) = total2
  nprint = nprint+2
end if
if (btest(iciweights,1)) then
  cprint(nprint+1) = total3
  cprint(nprint+2) = total4
  nprint = nprint+2
end if
if (btest(iciweights,2)) then
  cprint(nprint+1) = total5
  cprint(nprint+2) = total6
  nprint = nprint+2
end if
write(u6,formAD) ' Total all       : ',(cprint(mp),mp=1,nprint)
write(u6,'(a)') ' '

call mma_deallocate(wghtion1)
call mma_deallocate(wghtion2)
call mma_deallocate(wghtion3)
call mma_deallocate(wghtion4)
call mma_deallocate(wghtion5)
call mma_deallocate(wghtion6)

return

end subroutine ciweight2_cvb
