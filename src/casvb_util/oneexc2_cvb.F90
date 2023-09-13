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

subroutine oneexc2_cvb(cfrom,cto,vij,i1alf,i1bet,iato,ibto,phato,phbto,iapr,ixapr,ibpr,ixbpr,npvb,nda,ndb,n1a,n1b,nam1,nbm1,norb, &
                       commut,sc,absym,diag,idens,iPvb)
! Calculates Cto = Pvb Eij Cfrom

implicit real*8(a-h,o-z)
logical commut, sc, absym, diag
dimension cfrom(nda,ndb), cto(nda,ndb)
dimension vij(*)
dimension i1alf(n1a,norb), i1bet(n1b,norb)
dimension iato(norb,0:nam1), ibto(norb,0:nbm1)
dimension phato(norb,nam1), phbto(norb,nbm1)
dimension iapr(npvb), ixapr(nda+1), ibpr(npvb), ixbpr(ndb+1)
save two, half, thresh
data two/2d0/,half/.5d0/,thresh/1.d-10/

if (diag) then
  nvij = norb*norb
else if (idens == 1) then
  nvij = norb*(norb-1)
end if
if (absym) then
  if (idens == 0) then
    call dscal_(nda*ndb,half,cto,1)
  else
    call dscal_(nvij,half,vij,1)
  end if
end if
iprm = 0
do jorb=1,norb
  do iorb=1,norb
    if ((iorb == jorb) .and. (.not. diag)) goto 1200
    iprm = iprm+1
    if ((idens == 0) .and. (abs(vij(iprm)) < thresh)) goto 1200
    if (.not. sc) then
      ! a) Alpha excitation
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          if (idens == 0) then
            tcof = vij(iprm)*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
            if (iPvb == 0) then
              call daxpy_(ndb,tcof,cfrom(jax,1),nda,cto(iax,1),nda)
            else if (iPvb == 1) then
              do ixa=ixapr(jax),ixapr(jax+1)-1
                ibx = iapr(ixa)
                cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(jax,ibx)
              end do
            else if (iPvb == 2) then
              do ixa=ixapr(iax),ixapr(iax+1)-1
                ibx = iapr(ixa)
                cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(jax,ibx)
              end do
            end if
          else
            tcof = phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
            if (iPvb == 0) then
              vij(iprm) = vij(iprm)+tcof*ddot_(ndb,cto(iax,1),nda,cfrom(jax,1),nda)
            else if (iPvb == 1) then
              do ixa=ixapr(jax),ixapr(jax+1)-1
                ibx = iapr(ixa)
                vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
              end do
            else if (iPvb == 2) then
              do ixa=ixapr(iax),ixapr(iax+1)-1
                ibx = iapr(ixa)
                vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
              end do
            end if
          end if
        end if
      end do

      if (.not. absym) then
        ! c) Beta excitation
        do ib=1,n1b
          ibxtmp = i1bet(ib,iorb)
          jbx = ibto(jorb,ibxtmp)
          if (jbx /= 0) then
            ibx = ibto(iorb,ibxtmp)
            if (idens == 0) then
              tcof = vij(iprm)*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
              if (iPvb == 0) then
                call daxpy_(nda,tcof,cfrom(1,jbx),1,cto(1,ibx),1)
              else if (iPvb == 1) then
                do ixb=ixbpr(jbx),ixbpr(jbx+1)-1
                  iax = ibpr(ixb)
                  cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(iax,jbx)
                end do
              else if (iPvb == 2) then
                do ixb=ixbpr(ibx),ixbpr(ibx+1)-1
                  iax = ibpr(ixb)
                  cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(iax,jbx)
                end do
              end if
            else
              tcof = phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
              if (iPvb == 0) then
                vij(iprm) = vij(iprm)+tcof*ddot_(nda,cto(1,ibx),1,cfrom(1,jbx),1)
              else if (iPvb == 1) then
                do ixb=ixbpr(jbx),ixbpr(jbx+1)-1
                  iax = ibpr(ixb)
                  vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
                end do
              else if (iPvb == 2) then
                do ixb=ixbpr(ibx),ixbpr(ibx+1)-1
                  iax = ibpr(ixb)
                  vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
                end do
              end if
            end if
          end if
        end do
      end if
    else if (sc) then
      ! a) Alpha excitation
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          if (idens == 0) then
            tcof = vij(iprm)*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
            if (iPvb == 0) then
              call daxpy_(ndb,tcof,cfrom(jax,1),nda,cto(iax,1),nda)
            else if (iPvb == 1) then
              ibx = ndb-jax+1
              cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(jax,ibx)
            else if (iPvb == 2) then
              ibx = ndb-iax+1
              cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(jax,ibx)
            end if
          else
            tcof = phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
            if (iPvb == 0) then
              vij(iprm) = vij(iprm)+tcof*ddot_(ndb,cto(iax,1),nda,cfrom(jax,1),nda)
            else if (iPvb == 1) then
              ibx = ndb-jax+1
              vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
            else if (iPvb == 2) then
              ibx = ndb-iax+1
              vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(jax,ibx)
            end if
          end if
        end if
      end do

      if (.not. absym) then
        ! c) Beta excitation
        do ib=1,n1b
          ibxtmp = i1bet(ib,iorb)
          jbx = ibto(jorb,ibxtmp)
          if (jbx /= 0) then
            ibx = ibto(iorb,ibxtmp)
            if (idens == 0) then
              tcof = vij(iprm)*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
              if (iPvb == 0) then
                call daxpy_(nda,tcof,cfrom(1,jbx),1,cto(1,ibx),1)
              else if (iPvb == 1) then
                iax = nda-jbx+1
                cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(iax,jbx)
              else if (iPvb == 2) then
                iax = nda-ibx+1
                cto(iax,ibx) = cto(iax,ibx)+tcof*cfrom(iax,jbx)
              end if
            else
              tcof = phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
              if (iPvb == 0) then
                vij(iprm) = vij(iprm)+tcof*ddot_(nda,cto(1,ibx),1,cfrom(1,jbx),1)
              else if (iPvb == 1) then
                iax = nda-jbx+1
                vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
              else if (iPvb == 2) then
                iax = nda-ibx+1
                vij(iprm) = vij(iprm)+tcof*cto(iax,ibx)*cfrom(iax,jbx)
              end if
            end if
          end if
        end do
      end if
    end if
1200 continue
  end do
end do
if (absym) then
  if (idens == 0) then
    do ia=1,nda
      do ib=ia,nda
        cto(ia,ib) = cto(ia,ib)+cto(ib,ia)
        cto(ib,ia) = cto(ia,ib)
      end do
    end do
  else
    call dscal_(nvij,two,vij,1)
  end if
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_logical(commut)

end subroutine oneexc2_cvb
