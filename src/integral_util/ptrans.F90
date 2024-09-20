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

!#define _DEBUGPRINT_
subroutine ptrans(npam,ipam,nxpam,PSOPam,nPSOPam,Cred,nC,Scr1,nS1,Scr2,nS2)
! -------------------------------------------------------------------
! The following subroutine calculates a batch of 2-electron
! density matrix elements in SO basis.
! Need to know:
!   nish(),nash(),nbas()
!   npam(indpos,isym)= Nr of SO indices at index position 1..4,
!       with symmetry label isym=0..7.
!   ipam()= A consecutive list of SO indices.
! Also:
!   mxpam=Similar, ipam array
!   mxSO=Largest batch of SO indices in one single symmetry
! Returns:
!   PSOPam()=a four-index array containing the selected matrix
!         elements.
! -------------------------------------------------------------------

use Index_Functions, only: iTri, nTri_Elem
use etwas, only: CoulFac, mBas, mIrrep, nAsh, nIsh, npSOp
use pso_stuff, only: CMO, D0, G1, G2
use Constants, only: Zero, One, Two, Quart
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: npam(4,0:*), nxpam, nPSOPam, nC, nS1, nS2
real(kind=wp), intent(in) :: ipam(nxpam)
real(kind=wp), intent(out) :: PSOPam(nPSOPam), Cred(nC), Scr1(nS1), Scr2(nS2)
integer(kind=iwp) :: i, iEnd, ijSym, Ind, iOCMOI, iOCMOJ, iOCMOK, iOCMOL, iOCMOT, iOCMOU, iOCMOV, iOCMOX, ioDQ, ioDR, iods, iOff1, &
                     iOff2, ioPam1, ioPam2, ioPam3, ioPam4, ip, ipq, ipr, ips, ipSO, iq, ir, irq, irs, is, iScr, isq, iSta, iSym, &
                     it, itEnd, itSta, itu, ituvx, itv, itx, iu, iuEnd, iuSta, iv, iVEnd, ivSta, ivu, ivx, ix, ixEnd, ixSta, ixu, &
                     j, jEnd, jklOf1, jklOff, jSta, jSym, k, kEnd, klOf1, klOff, klSym, kSta, kSym, l, lEnd, lOf1, lOff, lSta, &
                     lSym, nbi, nbj, nbk, nbl, nCopy, ni, nijkl, nj, njkl, nk, nkl, nKLT, nl, nLTU, nnPam1, nnPam2, nnPam3, &
                     nnPam4, nSkip1, nSkip2, nt, nTUV, nu, nv, nx, nxv, nxvu, nxvut
real(kind=wp) :: Fact, t14

#ifdef _DEBUGPRINT_
write(u6,*) ' iPam',iPam
write(u6,*) ' nPam',(nPam(1,i),i=0,mIrrep-1)
write(u6,*) ' nPam',(nPam(2,i),i=0,mIrrep-1)
write(u6,*) ' nPam',(nPam(3,i),i=0,mIrrep-1)
write(u6,*) ' nPam',(nPam(4,i),i=0,mIrrep-1)
#endif
t14 = Quart
! Offsets into the ipam array:
nnpam1 = sum(npam(1,0:mirrep-1))
nnpam2 = sum(npam(2,0:mirrep-1))
nnpam3 = sum(npam(3,0:mirrep-1))
nnpam4 = sum(npam(4,0:mirrep-1))
nPSOP = nnpam1*nnpam2*nnpam3*nnpam4
PSOPam(1:nPSOP) = Zero
iopam1 = 0
iopam2 = nnpam1
iopam3 = iopam2+nnpam2
iopam4 = iopam3+nnpam3
! Loop over all possible symmetry combinations.
lend = 0
ixend = 0
iocmol = 0
ioDs = 0
do lsym=0,mirrep-1
  nl = npam(4,lsym)
  lsta = lend+1
  lend = lend+nl
  nx = nash(lsym)
  ixsta = ixend+1
  ixend = ixend+nx
  iocmox = iocmol+nish(lsym)*mbas(lsym)
  kend = 0
  ivend = 0
  iocmok = 0
  ioDr = 0
  do ksym=0,mirrep-1
    klsym = ieor(ksym,lsym)
    nk = npam(3,ksym)
    ksta = kend+1
    kend = kend+nk
    nkl = nk*nl
    nv = nash(ksym)
    ivsta = ivend+1
    ivend = ivend+nv
    nxv = nx*nv
    iocmov = iocmok+nish(ksym)*mbas(ksym)
    jend = 0
    iuend = 0
    iocmoj = 0
    ioDq = 0
    do jsym=0,mirrep-1
      nj = npam(2,jsym)
      jsta = jend+1
      jend = jend+nj
      njkl = nj*nkl
      nu = nash(jsym)
      iusta = iuend+1
      iuend = iuend+nu
      nxvu = nxv*nu
      iocmou = iocmoj+nish(jsym)*mbas(jsym)
      iend = 0
      itend = 0
      iocmoi = 0
      do isym=0,mirrep-1
        ni = npam(1,isym)
        ista = iend+1
        iend = iend+ni
        nijkl = ni*njkl
        nt = nash(isym)
        itsta = itend+1
        itend = itend+nt
        nxvut = nxvu*nt
        iocmot = iocmoi+nish(isym)*mbas(isym)
        ! Break loop if not acceptable symmetry combination.
        ! Break loop if no such symmetry block:
        ijsym = ieor(isym,jsym)
        if ((klsym == ijsym) .and. (nijkl /= 0)) then
#         ifdef _DEBUGPRINT_
          write(u6,*) ' i,j,k,lsym=',iSym,jSym,kSym,lSym
#         endif
          ! Bypass transformation if no active orbitals:
          if (nxvut /= 0) then
            ! Pick up matrix elements and put in a full symmetry block:
            ind = 0
            do ix=ixsta,ixend
              do iv=ivsta,ivend
                ivx = iTri(iv,ix)
                do iu=iusta,iuend
                  do it=itsta,itend
                    itu = iTri(it,iu)
                    ituvx = iTri(itu,ivx)
                    ind = ind+1
                    scr1(ind) = G2(ituvx,1)
                    if (isym == jsym) then
                      fact = One
                      if ((itu >= ivx) .and. (iv == ix)) fact = Two
                      if ((itu < ivx) .and. (it == iu)) fact = Two
                      scr1(ind) = fact*scr1(ind)
                      !hjw multiplying the G1 product with coulfac gives wrong result
                      scr1(ind) = scr1(ind)-G1(itu,1)*G1(ivx,1)
                    end if
                    if (isym == lsym) then
                      itx = iTri(it,ix)
                      ivu = iTri(iv,iu)
                      !hjw t14 includes exfac, why not coulfac above? What are these terms?
                      scr1(ind) = scr1(ind)+t14*G1(itx,1)*G1(ivu,1)
                    end if
                    if (isym == ksym) then
                      itv = iTri(it,iv)
                      ixu = iTri(ix,iu)
                      scr1(ind) = scr1(ind)+t14*G1(itv,1)*G1(ixu,1)
                    end if
                  end do
                end do
              end do
            end do
#           ifdef _DEBUGPRINT_
            call RecPrt('G2-G1G1(MO)',' ',Scr1,nash(iSym)*nash(jSym),nash(kSym)*nash(lsym))
#           endif

            ! Transform:
            !  scr2(l,tuv)= sum cmo(sl,x)*scr1(tuv,x)
            ncopy = nash(lsym)
            nskip1 = mbas(lsym)
            ioff2 = 0
            nskip2 = npam(4,lsym)
            ntuv = nash(isym)*nash(jsym)*nash(ksym)
            do l=lsta,lend
              ioff1 = iocmox+int(ipam(iopam4+l))
              ioff2 = ioff2+1
              call dcopy_(ncopy,CMO(ioff1,1),nskip1,Cred(ioff2),nskip2)
            end do
            call DGEMM_('N','T',nskip2,ntuv,ncopy,One,Cred,nskip2,Scr1,ntuv,Zero,Scr2,nskip2)
            ! Transform:
            !  scr3(k,ltu)= sum cmo(rk,v)*scr2(ltu,v)
            ncopy = nash(ksym)
            nskip1 = mbas(ksym)
            ioff2 = 0
            nskip2 = npam(3,ksym)
            nltu = nash(isym)*nash(jsym)*npam(4,lsym)
            do k=ksta,kend
              ioff1 = iocmov+int(ipam(iopam3+k))
              ioff2 = ioff2+1
              call dcopy_(ncopy,CMO(ioff1,1),nskip1,Cred(ioff2),nskip2)
            end do
            call DGEMM_('N','T',nskip2,nltu,ncopy,One,Cred,nskip2,Scr2,nltu,Zero,Scr1,nskip2)
            ! Transform:
            !  scr4(j,klt)= sum cmo(qj,u)*scr3(klt,u)
            ncopy = nash(jsym)
            nskip1 = mbas(jsym)
            ioff2 = 0
            nskip2 = npam(2,jsym)
            nklt = nash(isym)*npam(3,ksym)*npam(4,lsym)
            do j=jsta,jend
              ioff1 = iocmou+int(ipam(iopam2+j))
              ioff2 = ioff2+1
              call dcopy_(ncopy,CMO(ioff1,1),nskip1,Cred(ioff2),nskip2)
            end do
            call DGEMM_('N','T',nskip2,nklt,ncopy,One,Cred,nskip2,Scr1,nklt,Zero,Scr2,nskip2)
            ! Transform:
            !  scr5(i,jkl)= sum cmo(pi,t)*scr4(jkl,t)
            ncopy = nash(isym)
            nskip1 = mbas(isym)
            ioff2 = 0
            nskip2 = npam(1,isym)
            njkl = npam(2,jsym)*npam(3,ksym)*npam(4,lsym)
            do i=ista,iend
              ioff1 = iocmot+int(ipam(iopam1+i))
              ioff2 = ioff2+1
              call dcopy_(ncopy,CMO(ioff1,1),nskip1,Cred(ioff2),nskip2)
            end do
            call DGEMM_('N','T',nskip2,njkl,ncopy,One,Cred,nskip2,Scr2,njkl,Zero,Scr1,nskip2)
#           ifdef _DEBUGPRINT_
            call RecPrt('G2-G1G1(SO)',' ',Scr1,nPam(1,iSym)*nPam(2,jSym),nPam(3,kSym)*nPam(4,lSym))
#           endif

            ! Put results into correct positions in PSOPam:
            do l=lsta,lend
              loff = nnpam3*(l-1)
              lof1 = nPam(3,ksym)*(l-lsta)
              do k=ksta,kend
                kloff = nnpam2*(k-1+loff)
                klof1 = nPam(2,jsym)*(k-ksta+lof1)
                do j=jsta,jend
                  jkloff = nnpam1*(j-1+kloff)
                  jklof1 = nPam(1,isym)*(j-jsta+klof1)
                  do i=ista,iend
                    ipso = i+jkloff
                    iscr = 1+i-ista+jklof1
                    PSOPam(ipso) = Scr1(iscr)
                  end do
                end do
              end do
            end do
          end if
          ! Add contributions from 1-el density matrix:
          do l=lsta,lend
            is = int(ipam(iopam4+l))
            loff = nnpam3*(l-1)
            do k=ksta,kend
              ir = int(ipam(iopam3+k))
              irs = iTri(ir,is)
              kloff = nnpam2*(k-1+loff)
              do j=jsta,jend
                iq = int(ipam(iopam2+j))
                jkloff = nnpam1*(j-1+kloff)
                do i=ista,iend
                  ip = int(ipam(iopam1+i))
                  ipq = iTri(ip,iq)
                  ipso = i+jkloff
                  if (isym == lsym) then
                    ips = iTri(ip,is)
                    irq = iTri(ir,iq)
                    PSOPam(ipso) = PSOPam(ipso)-t14*D0(ioDs+ips,1)*D0(ioDr+irq,1)
                  end if
                  if (isym == ksym) then
                    ipr = iTri(ip,ir)
                    isq = iTri(is,iq)
                    PSOPam(ipso) = PSOPam(ipso)-t14*D0(ioDr+ipr,1)*D0(ioDs+isq,1)
                  end if
                  if (isym == jsym) then
                    PSOPam(ipso) = PSOPam(ipso)+D0(ioDq+ipq,1)*D0(ioDs+irs,1)*coulfac
                  end if
                end do
              end do
            end do
          end do
        end if
        ! End of loop over symmetry labels.
        nbi = mbas(isym)
        iocmoi = iocmoi+nbi**2
      end do
      nbj = mbas(jsym)
      iocmoj = iocmoj+nbj**2
      ioDq = ioDq+nTri_Elem(nbj)
    end do
    nbk = mbas(ksym)
    iocmok = iocmok+nbk**2
    ioDr = ioDr+nTri_Elem(nbk)
  end do
  nbl = mbas(lsym)
  iocmol = iocmol+nbl**2
  ioDs = ioDs+nTri_Elem(nbl)
end do
#ifdef _DEBUGPRINT_
call RecPrt('PSOPam',' ',PSOPam,nnPam1*nnPam2,nnPam3*nnPam4)
#endif

return

end subroutine ptrans
