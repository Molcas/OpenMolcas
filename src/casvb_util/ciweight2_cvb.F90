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
subroutine ciweight2_cvb(civec,civbs,civb,citmp,civec5,orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3,vec1,vec2,vec3,vec4,vec5, &
                         wghtion1,wghtion2,wghtion3,wghtion4,wghtion5,wghtion6,mingrph,maxgrph,xalf,xbet,iaocc,ibocc,mingion, &
                         maxgion,nkion,xion,locion,lunion,mingsng,maxgsng,nksng,xsng,locsng,lunsng,mingasg,maxgasg,nkasg,xasg, &
                         locasg,lunasg,gal1,gal2,indavec,indbvec,ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)

implicit real*8(a-h,o-w,y-z),integer(x)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "formats_cvb.fh"
character*240 line
dimension civec(ndet), civbs(ndet), civb(ndet), citmp(ndet), civec5(ndet)
dimension orbs(norb,norb), sorbs(norb,norb)
dimension orbinv(norb,norb), owrk(norb,norb)
dimension gjorb(*), gjorb2(*), gjorb3(*)
dimension vec1(ndet), vec2(ndet), vec3(ndet), vec4(ndet), vec5(ndet)
dimension wghtion1(ionmin:ionmax), wghtion2(ionmin:ionmax)
dimension wghtion3(ionmin:ionmax), wghtion4(ionmin:ionmax)
dimension wghtion5(ionmin:ionmax), wghtion6(ionmin:ionmax)
dimension mingrph(0:norb), maxgrph(0:norb)
dimension xalf(0:norb,0:nalf), xbet(0:norb,0:nbet)
dimension iaocc(norb), ibocc(norb)
dimension mingion(0:norb), maxgion(0:norb), nkion(0:norb)
dimension xion((norb+1)*(ionmax+1)), locion(norb), lunion(norb)
dimension mingsng(0:norb), maxgsng(0:norb), nksng(0:norb)
dimension xsng((mxrem+1)*(mxsng+1)), locsng(norb), lunsng(norb)
dimension mingasg(0:norb), maxgasg(0:norb), nkasg(0:norb)
dimension xasg((mxsng+1)*(mxasg+1)), locasg(norb), lunasg(norb)
dimension gal1(ncnfcas), gal2(ncnfcas)
dimension indavec(mxdetcas), indbvec(mxdetcas)
dimension cprint(6)
integer rc

call cidot_cvb(civb,civbs,cnrm)
fac = svb/sqrt(cnrm)

call cicopy_cvb(civec,citmp)
call fmove_cvb(orbs,orbinv,norb*norb)
call mxinv_cvb(orbinv,norb)
call gaussj_cvb(orbinv,gjorb)
call applyt_cvb(civec,gjorb)
! Chirgwin-Coulson weights
if (mod(iciweights,2) == 1) then
  call transp_cvb(orbs,owrk,norb,norb)
  call gaussj_cvb(owrk,gjorb2)
  call applyt_cvb(citmp,gjorb2)
  do idet=1,ndet
    vec2(idet) = (vec1(idet)-fac*vec2(idet))*(vec3(idet)-fac*vec4(idet))
  end do
  do idet=1,ndet
    vec1(idet) = vec1(idet)*vec3(idet)
  end do
end if
do idet=1,ndet
  vec4(idet) = vec3(idet)-fac*vec4(idet)
end do

! Inverse-overlap weights
if (mod(iciweights,8) > 3) then
  call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
  call fmove_cvb(sorbs,orbinv,norb*norb)
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
    call imove_cvb(maxgion,nkion,norb+1)
    call occupy_cvb(nkion,norb,locion,lunion)
    ! Initialise loop for singly occupied orbitals
    do iorb=0,mrem
      mingsng(iorb) = max(iorb-mrem+nsing,0)
      maxgsng(iorb) = min(iorb,nsing)
    end do
    call weight_cvb(xsng,mingsng,maxgsng,nsing,mrem)
    call imove_cvb(maxgsng,nksng,mrem+1)
    call occupy_cvb(nksng,mrem,locsng,lunsng)
    ! Initialise loop for singly occupied alpha orbitals
    do iorb=0,nsing
      mingasg(iorb) = max(iorb-nsing+nalfsng,0)
      maxgasg(iorb) = min(iorb,nalfsng)
    end do
    call weight_cvb(xasg,mingasg,maxgasg,nalfsng,nsing)
    call imove_cvb(maxgasg,nkasg,nsing+1)
    call occupy_cvb(nkasg,nsing,locasg,lunasg)

    ! Loop ionic
    indion = 1
    do
      ! Loop singly occupied
      indsng = 1
      do
        call fzero(vec5,ndet)
        s11 = zero
        s22 = zero
        s12 = zero
        ! Loop singly occupied alpha
        indasg = 1

        do
          call izero(iaocc,norb)
          call izero(ibocc,norb)
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
          vec5(indab) = vec3(indab)
          s11 = s11+vec3(indab)*vec3(indab)
          s22 = s22+vec4(indab)*vec4(indab)
          s12 = s12+vec3(indab)*vec4(indab)

          call loind_cvb(nsing,nalfsng,nkasg,mingasg,maxgasg,locasg,lunasg,indasg,xasg,rc)
          if (rc == 0) exit
        end do

        call applyt_cvb(civec5,gjorb)

        call icomb_cvb(nsing,nalfsng,nindasg)
        sm1 = zero
        do indasg=1,nindasg
          inda = indavec(indasg)
          indb = indbvec(indasg)
          indab = (indb-1)*nda+inda
          sm1 = sm1+vec3(indab)*vec5(indab)
        end do

        if (abs(sm1) > 1.d-20) then
          sm1 = s11*s11/sm1
        else if (abs(sm1) <= 1.d-20) then
          sm1 = zero
        end if

        call fzero(vec5,ndet)
        do indasg=1,nindasg
          inda = indavec(indasg)
          indb = indbvec(indasg)
          indab = (indb-1)*nda+inda
          vec5(indab) = vec4(indab)
        end do

        call applyt_cvb(civec5,gjorb)

        sm2 = zero
        do indasg=1,nindasg
          inda = indavec(indasg)
          indb = indbvec(indasg)
          indab = (indb-1)*nda+inda
          sm2 = sm2+vec4(indab)*vec5(indab)
        end do

        if (abs(sm2) > 1.d-20) then
          sm2 = s22*s22/sm2
        else if (abs(sm2) <= 1.d-20) then
          sm2 = zero
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
  sum1 = zero
  sum2 = zero
  do ic=1,ncnfcas
    sum1 = sum1+gal1(ic)
    sum2 = sum2+gal2(ic)
  end do
  fac1 = one/sum1
  if ((abs(one-svb*svb) < 1.d-20) .and. (abs(sum2) < 1.d-20)) then
    fac2 = one
  else
    fac2 = (one-svb*svb)/sum2
  end if
  do ic=1,ncnfcas
    gal1(ic) = fac1*gal1(ic)
    gal2(ic) = fac2*gal2(ic)
  end do
end if

! Weights of Lowdin orthonormalized structures
if (mod(iciweights,4) > 1) then
  call mxattb_cvb(orbs,orbs,norb,norb,norb,sorbs)
  call mxsqrt_cvb(sorbs,norb,1)
  call gaussj_cvb(sorbs,gjorb3)
  call applyt_cvb(civec,gjorb3)
  call applyt_cvb(civb,gjorb3)
  call cidot_cvb(civb,civb,cnrm)
  fac = svb/sqrt(cnrm)
  do idet=1,ndet
    vec4(idet) = vec4(idet)*vec4(idet)
    vec3(idet) = vec3(idet)*vec3(idet)
  end do
end if

write(6,'(/,2a)') ' Weights of CASSCF configurations in VB basis (c_res=c_cas-Svb*c_vb) :'
write(6,'(2a)') ' ---------------------------------------------------------------------'
if (mod(iciweights,8) > 3) then
  write(6,'(a)') ' Sum of inverse-overlap weights :'
  write(6,form2AD) ' c_cas :',sum1,' expected :',one
  write(6,form2AD) ' c_res :',sum2,' expected :',one-svb*svb
  write(6,'(a)') ' '
end if
lenfld = 8+iprec
ix1 = max(0,min(3,2*lenfld-16))
call cblank_cvb(line,240)
if ((npcf > 0) .or. (npcf == -1)) then
  line(1:21) = '  Conf. =>  Orbitals '
  ibeg = max(14+3*nel,22)
  ibegt = ibeg
  if (mod(iciweights,2) == 1) then
    line(ibegt+ix1:ibegt+2*lenfld-1) = 'Chirgwin-Coulson'
    ibegt = ibegt+2*lenfld
  end if
  if (mod(iciweights,4) > 1) then
    line(ibegt+ix1:ibegt+2*lenfld-1) = 'Lowdin'
    ibegt = ibegt+2*lenfld
  end if
  if (mod(iciweights,8) > 3) then
    line(ibegt+ix1:ibegt+2*lenfld-1) = 'Inverse'
    ibegt = ibegt+2*lenfld
  end if
  write(6,'(a)') line(1:len_trim_cvb(line))
  call cblank_cvb(line,240)
  if (mod(iciweights,2) == 1) then
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
    ibeg = ibeg+lenfld
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
    ibeg = ibeg+lenfld
  end if
  if (mod(iciweights,4) > 1) then
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
    ibeg = ibeg+lenfld
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
    ibeg = ibeg+lenfld
  end if
  if (mod(iciweights,8) > 3) then
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
    ibeg = ibeg+lenfld
    line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
    ibeg = ibeg+lenfld
  end if
  write(6,'(a)') line(1:len_trim_cvb(line))
end if
call cblank_cvb(line,240)

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
call fzero(wghtion1,ionmax-ionmin+1)
call fzero(wghtion2,ionmax-ionmin+1)
call fzero(wghtion3,ionmax-ionmin+1)
call fzero(wghtion4,ionmax-ionmin+1)
call fzero(wghtion5,ionmax-ionmin+1)
call fzero(wghtion6,ionmax-ionmin+1)
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
  call imove_cvb(maxgion,nkion,norb+1)
  call occupy_cvb(nkion,norb,locion,lunion)
  ! Initialise loop for singly occupied orbitals
  do iorb=0,mrem
    mingsng(iorb) = max(iorb-mrem+nsing,0)
    maxgsng(iorb) = min(iorb,nsing)
  end do
  call weight_cvb(xsng,mingsng,maxgsng,nsing,mrem)
  call imove_cvb(maxgsng,nksng,mrem+1)
  call occupy_cvb(nksng,mrem,locsng,lunsng)
  ! Initialise loop for singly occupied alpha orbitals
  do iorb=0,nsing
    mingasg(iorb) = max(iorb-nsing+nalfsng,0)
    maxgasg(iorb) = min(iorb,nalfsng)
  end do
  call weight_cvb(xasg,mingasg,maxgasg,nalfsng,nsing)
  call imove_cvb(maxgasg,nkasg,nsing+1)
  call occupy_cvb(nkasg,nsing,locasg,lunasg)

  ! Loop ionic
  indion = 1
  do
    ! Loop singly occupied
    indsng = 1
    do
      c1 = zero
      c2 = zero
      c3 = zero
      c4 = zero
      ! Loop singly occupied alpha
      indasg = 1

      do
        call izero(iaocc,norb)
        call izero(ibocc,norb)
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
        c1 = c1+vec1(indab)
        c2 = c2+vec2(indab)
        c3 = c3+vec3(indab)
        c4 = c4+vec4(indab)

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
        if (mod(iciweights,2) == 1) then
          cprint(nprint+1) = c1
          cprint(nprint+2) = c2
          nprint = nprint+2
        end if
        if (mod(iciweights,4) > 1) then
          cprint(nprint+1) = c3
          cprint(nprint+2) = c4
          nprint = nprint+2
        end if
        if (mod(iciweights,8) > 3) then
          cprint(nprint+1) = gal1(nc)
          cprint(nprint+2) = gal2(nc)
          nprint = nprint+2
        end if
        write(6,formAD) line(1:ilin),(cprint(mp),mp=1,nprint)
      end if
      call loind_cvb(mrem,nsing,nksng,mingsng,maxgsng,locsng,lunsng,indsng,xsng,rc)
      if (rc == 0) exit
    end do
    call loind_cvb(norb,ion,nkion,mingion,maxgion,locion,lunion,indion,xion,rc)
    if (rc == 0) exit
  end do
end do

call cblank_cvb(line,240)
line(1:21) = ' Accumulated weights:'
ibeg = 22
if (mod(iciweights,2) == 1) then
  line(ibeg+ix1:ibeg+2*lenfld-1) = 'Chirgwin-Coulson'
  ibeg = ibeg+2*lenfld
end if
if (mod(iciweights,4) > 1) then
  line(ibeg+ix1:ibeg+2*lenfld-1) = 'Lowdin'
  ibeg = ibeg+2*lenfld
end if
if (mod(iciweights,8) > 3) then
  line(ibeg+ix1:ibeg+2*lenfld-1) = 'Inverse'
  ibeg = ibeg+2*lenfld
end if
write(6,'(/,a)') line(1:len_trim_cvb(line))
call cblank_cvb(line,240)
ibeg = 22
if (mod(iciweights,2) == 1) then
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
  ibeg = ibeg+lenfld
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
  ibeg = ibeg+lenfld
end if
if (mod(iciweights,4) > 1) then
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
  ibeg = ibeg+lenfld
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
  ibeg = ibeg+lenfld
end if
if (mod(iciweights,8) > 3) then
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_cas'
  ibeg = ibeg+lenfld
  line(ibeg+ix1:ibeg+lenfld-1) = 'c_res'
  ibeg = ibeg+lenfld
end if
write(6,'(a)') line(1:len_trim_cvb(line))
call cblank_cvb(line,240)
total1 = zero
total2 = zero
total3 = zero
total4 = zero
total5 = zero
total6 = zero
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
  if (mod(iciweights,2) == 1) then
    cprint(nprint+1) = wghtion1(ion)
    cprint(nprint+2) = wghtion2(ion)
    nprint = nprint+2
  end if
  if (mod(iciweights,4) > 1) then
    cprint(nprint+1) = wghtion3(ion)
    cprint(nprint+2) = wghtion4(ion)
    nprint = nprint+2
  end if
  if (mod(iciweights,8) > 3) then
    cprint(nprint+1) = wghtion5(ion)
    cprint(nprint+2) = wghtion6(ion)
    nprint = nprint+2
  end if
  write(6,formAD) line(1:19),(cprint(mp),mp=1,nprint)
end do
nprint = 0
if (mod(iciweights,2) == 1) then
  cprint(nprint+1) = total1
  cprint(nprint+2) = total2
  nprint = nprint+2
end if
if (mod(iciweights,4) > 1) then
  cprint(nprint+1) = total3
  cprint(nprint+2) = total4
  nprint = nprint+2
end if
if (mod(iciweights,8) > 3) then
  cprint(nprint+1) = total5
  cprint(nprint+2) = total6
  nprint = nprint+2
end if
write(6,formAD) ' Total all       : ',(cprint(mp),mp=1,nprint)
write(6,'(a)') ' '

return

end subroutine ciweight2_cvb
