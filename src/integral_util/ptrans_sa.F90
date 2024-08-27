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
subroutine ptrans_sa(npam,ipam,nxpam,PSOPam,nPSOPam,Cred,nC,Scr1,nS1,Scr2,nS2,ScrP,nsp)
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

use pso_stuff, only: nSSDM, SSDM, DSO => D0, CMO, G2
use Constants, only: Zero, Quart
use etwas, only: mIrrep, npSOp, mBas, nAsh, nIsh

implicit none
integer nxpam, nPSOPam, nC, nS1, nS2, nSP
integer npam(4,0:*), indi(4)
real*8 ipam(nxpam)
real*8 PSOPam(nPSOPam), Cred(*), Scr1(nS1), Scr2(nS2), ScrP(nsP)
integer i, j, i3adr
integer nnPam1, nnPam2, nnPam3, nnPam4, iSym, jSym, kSym, lSym, ioPam1, ioPam2, ioPam3, ioPam4, iEnd, jEnd, kEnd, lEnd, ni, nj, &
        nk, nl, ip, iq, is, it, ir, ipq, irs, ipr, ips, irq, isq, nbi, nbj, nbk, nbl, nx, nkl, nv, nxv, njkl, nu, nxvu, nt, &
        nxvut, ix, iv, iu, itu, ituvx, iScr, ixEnd, iOCMOL, iods, lSta, ixSta, iOCMOX, iVEnd, iOCMOK, ioDR, klSym, kSta, ivSta, &
        iOCMOV, iuEnd, iOCMOJ, ioDQ, jSta, iuSta, iOCMOU, itEnd, iOCMOI, iSta, nijkl, iOCMOT, ijSym, Ind, ivx, nCopy, nSkip1, &
        iOff2, nSkip2, nTUV, l, nLTU, k, nKLT, lOff, lOf1, klOff, klOf1, jklOff, jklOf1, itSta, iOff1, ipSO, ioIT, isSDM
real*8 Fact
! Triangular addressing without symmetry:
i3adr(i,j) = ((max(i,j))*((max(i,j))-1))/2+min(i,j)

! Offsets into the ipam array:
nnpam1 = 0
nnpam2 = 0
nnpam3 = 0
nnpam4 = 0
do isym=0,mirrep-1
  nnpam1 = nnpam1+npam(1,isym)
  nnpam2 = nnpam2+npam(2,isym)
  nnpam3 = nnpam3+npam(3,isym)
  nnpam4 = nnpam4+npam(4,isym)
end do
nPSOP = nnpam1*nnpam2*nnpam3*nnpam4
call dcopy_(nPSOP,[Zero],0,PSOPam,1)
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
        ijsym = ieor(isym,jsym)
        if (klsym /= ijsym) goto 1005
        ! Break loop if no such symmetry block:
        if (nijkl == 0) goto 1005
#       ifdef _DEBUGPRINT_
        write(6,*) ' i,j,k,lsym=',iSym,jSym,kSym,lSym
#       endif
        ! Bypass transformation if no active orbitals:
        if (nxvut == 0) goto 300
        ! Pick up matrix elements and put in a full symmetry block:
        ! State-average 2-DM
        ind = 0
        do ix=ixsta,ixend
          do iv=ivsta,ivend
            ivx = i3adr(iv,ix)
            do iu=iusta,iuend
              do it=itsta,itend
                itu = i3adr(it,iu)
                ituvx = i3adr(itu,ivx)
                ind = ind+1
                scrP(ind) = G2(ituvx,2)
                if (isym == jsym) then
                  fact = 1.0d00
                  if ((itu >= ivx) .and. (iv == ix)) fact = 2.0d00
                  if ((itu < ivx) .and. (it == iu)) fact = 2.0d00
                  scrP(ind) = fact*scrP(ind)
                end if
              end do
            end do
          end do
        end do

        ! Transform:
        !  scr2(l,tuv)= sum cmo(sl,x)*scr1(tuv,x)
        do ioit=1,4
          call icopy(4,[1],0,indi,1)
          indi(ioit) = 2
          ncopy = nash(lsym)
          nskip1 = mbas(lsym)
          ioff2 = 0
          nskip2 = npam(4,lsym)
          ntuv = nash(isym)*nash(jsym)*nash(ksym)
          do l=lsta,lend
            ioff1 = iocmox+int(ipam(iopam4+l))
            ioff2 = ioff2+1
            call dcopy_(ncopy,CMO(ioff1,indi(1)),nskip1,Cred(ioff2),nskip2)
          end do
          call DGEMM_('N','T',nskip2,ntuv,ncopy,1.0d0,Cred,nskip2,ScrP,ntuv,0.0d0,Scr2,nskip2)
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
            call dcopy_(ncopy,CMO(ioff1,indi(2)),nskip1,Cred(ioff2),nskip2)
          end do
          call DGEMM_('N','T',nskip2,nltu,ncopy,1.0d0,Cred,nskip2,Scr2,nltu,0.0d0,Scr1,nskip2)
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
            call dcopy_(ncopy,CMO(ioff1,indi(3)),nskip1,Cred(ioff2),nskip2)
          end do

          call DGEMM_('N','T',nskip2,nklt,ncopy,1.0d0,Cred,nskip2,Scr1,nklt,0.0d0,Scr2,nskip2)
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
            call dcopy_(ncopy,CMO(ioff1,indi(4)),nskip1,Cred(ioff2),nskip2)
          end do
          call DGEMM_('N','T',nskip2,njkl,ncopy,1.0d0,Cred,nskip2,Scr2,njkl,0.0d0,Scr1,nskip2)
#         ifdef _DEBUGPRINT_
          call RecPrt('G2(SO1)',' ',Scr1,nPam(1,iSym)*nPam(2,jSym),nPam(3,kSym)*nPam(4,lSym))
#         endif

          !============================================================*
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
                  PSOPam(ipso) = Scr1(iscr)+PSOPAM(ipSO)
                end do
              end do
            end do
          end do
        end do
        !==============================================================*
        ! State-specific 2-DM
        ind = 0
        do ix=ixsta,ixend
          do iv=ivsta,ivend
            ivx = i3adr(iv,ix)

            do iu=iusta,iuend
              do it=itsta,itend
                itu = i3adr(it,iu)
                ituvx = i3adr(itu,ivx)
                ind = ind+1
                scr1(ind) = G2(ituvx,1)
                if (isym == jsym) then
                  fact = 1.0d00
                  if ((itu >= ivx) .and. (iv == ix)) fact = 2.0d00
                  if ((itu < ivx) .and. (it == iu)) fact = 2.0d00
                  scr1(ind) = fact*scr1(ind)
                end if
              end do
            end do
          end do
        end do
#       ifdef _DEBUGPRINT_
        call RecPrt('G2(MO)',' ',Scr1,nash(iSym)*nash(jSym),nash(kSym)*nash(lsym))
#       endif

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
        call DGEMM_('N','T',nskip2,ntuv,ncopy,1.0d0,Cred,nskip2,Scr1,ntuv,0.0d0,Scr2,nskip2)
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
        call DGEMM_('N','T',nskip2,nltu,ncopy,1.0d0,Cred,nskip2,Scr2,nltu,0.0d0,Scr1,nskip2)
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

        call DGEMM_('N','T',nskip2,nklt,ncopy,1.0d0,Cred,nskip2,Scr1,nklt,0.0d0,Scr2,nskip2)
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
        call DGEMM_('N','T',nskip2,njkl,ncopy,1.0d0,Cred,nskip2,Scr2,njkl,0.0d0,Scr1,nskip2)
#       ifdef _DEBUGPRINT_
        call RecPrt('G2(SO2)',' ',Scr1,nPam(1,iSym)*nPam(2,jSym),nPam(3,kSym)*nPam(4,lSym))
        call RecPrt('PSOPam 0',' ',PSOPam,nnPam1*nnPam2,nnPam3*nnPam4)
#       endif

        !==============================================================*
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
                PSOPam(ipso) = Scr1(iscr)+PSOPAM(ipSO)
              end do
            end do
          end do
        end do
#       ifdef _DEBUGPRINT_
        call RecPrt('PSOPam no 1-el',' ',PSOPam,nnPam1*nnPam2,nnPam3*nnPam4)
#       endif

        !==============================================================*

        ! Add contributions from 1-el density matrix:

        !write(6,*) 'information before introducing 1-el density matrix' ! yma
        !write(6,*) 'lsta,lend',lsta,lend
        !write(6,*) 'ksta,kend',ksta,kend
        !write(6,*) 'jsta,jend',jsta,jend
        !write(6,*) 'ista,iend',ista,iend

        !do i=1,nDSO
        !  write(6,*) i,'DSO-1',DSO(i,1)
        !  write(6,*) i,'DSO-2',DSO(i,2)
        !  write(6,*) i,'DSO-3',DSO(i,3)
        !  write(6,*) i,'DSO-4',DSO(i,4)
        !  write(6,*) i,'DSO-5',DSO(i,5)
        !end do

300     continue
        do l=lsta,lend
          is = int(ipam(iopam4+l))
          loff = nnpam3*(l-1)
          do k=ksta,kend
            ir = int(ipam(iopam3+k))
            irs = i3adr(ir,is)
            kloff = nnpam2*(k-1+loff)
            do j=jsta,jend
              iq = int(ipam(iopam2+j))
              jkloff = nnpam1*(j-1+kloff)
              do i=ista,iend
                ip = int(ipam(iopam1+i))
                ipq = i3adr(ip,iq)
                ipso = i+jkloff

                ! ANDERS HAS ADDED LAGRANGIAN HERE
                ! BY USING A VECTOR OF DENSITIES, THIS CAN EASILY BE
                ! GENERALIZED TO CALCULATE FOR EXAMPLE THE E_2 CONTRIBUTION
                ! FOR RAMAN SPECTRA

                if (isym == lsym) then
                  ips = i3adr(ip,is)
                  irq = i3adr(ir,iq)
                  PSOPam(ipso) = PSOPam(ipso)-Quart*DSO(ioDs+ips,1)*DSO(ioDr+irq,2)-Quart*DSO(ioDs+ips,2)*DSO(ioDr+irq,1)- &
                                 Quart*DSO(ioDs+ips,3)*DSO(ioDr+irq,4)-Quart*DSO(ioDs+ips,4)*DSO(ioDr+irq,3)- &
                                 Quart*DSO(ioDs+ips,1)*DSO(ioDr+irq,6)-Quart*DSO(ioDs+ips,6)*DSO(ioDr+irq,1)
                  !ANDREW - uncomment
                  !-Quart*DSO(ioDs+ips,1)*DSO(ioDr+irq,5)-Quart*DSO(ioDs+ips,5)*DSO(ioDr+irq,1)
                  !END ANDREW
                  if (nSSDM /= 0) then
                    ! The last four lines subtract unnecessary contributions
                    do iSSDM=1,nSSDM
                      PSOPam(ipso) = PSOPam(ipso)-Quart*SSDM(ioDs+ips,1,iSSDM)*SSDM(ioDr+irq,2,iSSDM)- &
                                     Quart*SSDM(ioDs+ips,2,iSSDM)*SSDM(ioDr+irq,1,iSSDM)
                    end do
                  end if
                end if
                if (isym == ksym) then
                  ipr = i3adr(ip,ir)
                  isq = i3adr(is,iq)
                  PSOPam(ipso) = PSOPam(ipso)-Quart*DSO(ioDr+ipr,1)*DSO(ioDs+isq,2)-Quart*DSO(ioDr+ipr,2)*DSO(ioDs+isq,1)- &
                                 Quart*DSO(ioDr+ipr,3)*DSO(ioDs+isq,4)-Quart*DSO(ioDr+ipr,4)*DSO(ioDs+isq,3)- &
                                 Quart*DSO(ioDr+ipr,1)*DSO(ioDs+isq,6)-Quart*DSO(ioDr+ipr,6)*DSO(ioDs+isq,1)
                  !ANDREW - uncomment
                  !-Quart*DSO(ioDr+ipr,1)*DSO(ioDs+isq,5)-Quart*DSO(ioDr+ipr,5)*DSO(ioDs+isq,1)
                  !END ANDREW
                  if (nSSDM /= 0) then
                    do iSSDM=1,nSSDM
                      PSOPam(ipso) = PSOPam(ipso)-Quart*SSDM(ioDr+ipr,1,iSSDM)*SSDM(ioDs+isq,2,iSSDM)- &
                                     Quart*SSDM(ioDr+ipr,2,iSSDM)*SSDM(ioDs+isq,1,iSSDM)
                    end do
                  end if
                end if
                if (isym == jsym) then
                  PSOPam(ipso) = PSOPam(ipso)+DSO(ioDq+ipq,1)*DSO(ioDs+irs,2)+DSO(ioDq+ipq,2)*DSO(ioDs+irs,1)+ &
                                 DSO(ioDq+ipq,3)*DSO(ioDs+irs,4)+DSO(ioDq+ipq,4)*DSO(ioDs+irs,3)+DSO(ioDq+ipq,1)*DSO(ioDs+irs,5)+ &
                                 DSO(ioDq+ipq,5)*DSO(ioDs+irs,1)
                  if (nSSDM /= 0) then
                    issdm = 1
                    do iSSDM=1,nSSDM
                      PSOPam(ipso) = PSOPam(ipso)+SSDM(ioDq+ipq,1,iSSDM)*SSDM(ioDs+irs,2,iSSDM)+ &
                                     SSDM(ioDq+ipq,2,iSSDM)*SSDM(ioDs+irs,1,iSSDM)
                    end do
                  end if
                end if
                ! IT STOPS HERE
              end do
            end do
          end do
        end do
        ! End of loop over symmetry labels.
1005    continue
        nbi = mbas(isym)
        iocmoi = iocmoi+nbi**2
      end do
      nbj = mbas(jsym)
      iocmoj = iocmoj+nbj**2
      ioDq = ioDq+(nbj*(nbj+1))/2
    end do
    nbk = mbas(ksym)
    iocmok = iocmok+nbk**2
    ioDr = ioDr+(nbk*(nbk+1))/2
  end do
  nbl = mbas(lsym)
  iocmol = iocmol+nbl**2
  ioDs = ioDs+(nbl*(nbl+1))/2
end do
#ifdef _DEBUGPRINT_
call RecPrt('PSOPam',' ',PSOPam,nnPam1*nnPam2,nnPam3*nnPam4)
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nC)

end subroutine ptrans_sa
