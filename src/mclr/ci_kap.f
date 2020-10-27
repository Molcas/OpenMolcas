************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine CI_KAP(ipcid,fock,fockOut,isym)
*     use ipPage, only: W
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "dmrginfo_mclr.fh"
      Real*8 Fock(*),FockOut(*)

! Added for DMRG calculation
      real*8,allocatable::tmpDe(:,:),tmpP(:),tmpDeM(:,:),tmpPM(:,:,:,:)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      irc=ipnout(-1)
      Call Getmem('De','allo','real',ipde,ntash**2)
      Call Getmem('P','allo','real',ipp,n2dens)

      Call CIDens_SA(.true.,ipCI,ipCid,State_sym,State_Sym,
     &                Work(ipP),Work(ipde))

! ====================================================================
      If(doDMRG)then
        call dmrg_dim_change_mclr(LRras2(1:8),ntash,0)
        call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)

        ! will be changed to getmem/mma_allocate
        allocate(tmpDe(ndim,ndim))
        allocate(tmpP(ndim**2*(ndim**2+1)/2))
        allocate(tmpDeM(ntash,ntash))
        allocate(tmpPM(ntash,ntash,ntash,ntash))
        tmpDe=0.0d0
        tmpP=0.0d0
        tmpDeM=0.0d0
        tmpPM=0.0d0

        ij=0
        do i=1,ntash
          do j=1,ntash
            ij=ij+1
            if(abs(Work(ipDe+ij-1)).lt.1.0e-12)then
              Work(ipDe+ij-1)=0.0d0
            end if
            tmpDeM(i,j)=Work(ipDe+ij-1)
          end do
        end do

        ij=0
        do i=1,ndim
          do j=1,ndim
            ij=ij+1
            if(i.gt.ntash.or.j.gt.ntash)then
              tmpDe(i,j)=0.0d0
            else
              tmpDe(i,j)=tmpDeM(i,j)
            end if
          end do
        end do

        Do i=1,ntash
          Do j=1,ntash
            Do k=1,ntash
              Do l=1,ntash
                ij1=ntash*(i-1)+j
                ij2=ntash*(j-1)+i
                kl1=ntash*(k-1)+l
                kl2=ntash*(l-1)+k
                if(ij1.ge.kl1)then
                  if(abs(Work(ipP+itri(ij1,kl1)-1)).lt.1.0e-12)then
                    Work(ipP+itri(ij1,kl1)-1)=0.0d0
                  end if
                  tmpPM(i,j,k,l)=Work(ipP+itri(ij1,kl1)-1)
                end if
              End Do
            End Do
          End Do
        End Do

        do i=1,ndim
          do j=1,ndim
            do k=1,ndim
              do l=1,ndim
                ij1=ndim*(i-1)+j
                ij2=ndim*(j-1)+i
                kl1=ndim*(k-1)+l
                kl2=ndim*(l-1)+k
                if(ij1.ge.kl1)then
            if(i.gt.ntash.or.j.gt.ntash.or.k.gt.ntash.or.l.gt.ntash)then
                  tmpP(itri(ij1,kl1))=0.0d0
                else
                  tmpP(itri(ij1,kl1))=tmpPM(i,j,k,l)
                end if
                end if
              end do
            end do
          end do
        end do

! ====================================================================
        call dmrg_dim_change_mclr(RGras2(1:8),ntash,0)
        call dmrg_dim_change_mclr(RGras2(1:8),nna,0)

*       irc=ipin(ipCID)
*       irc=ipin(ipci)
*       call projecter(W(ipCID)%Vec,W(ipci)%Vec,
*    &                 Work(ipDe),Work(ipp))
        call dcopy_(ndens2,[0.0d0],0,Fock,1)
        call dcopy_(ndens2,[0.0d0],0,FockOut,1)
        d0=0.0d0

        Call FockGen(d0,tmpDe,tmpP,Fock,FockOut,isym) ! yma modified

      else

*
*       irc=ipin(ipCID)
*       irc=ipin(ipci)
*       call projecter(W(ipCID)%Vec,W(ipci)%Vec,
*    &                 Work(ipDe),Work(ipp))
        call dcopy_(ndens2,[0.0d0],0,Fock,1)
        call dcopy_(ndens2,[0.0d0],0,FockOut,1)
        d0=0.0d0
        Call FockGen(d0,Work(ipDe),Work(ipP),Fock,FockOut,isym)
      end if
      Call Getmem('De','free','real',ipde,ntash**2)
      Call Getmem('P','free','real',ipp,n2dens)
      Call Getmem('P','chec','real',ipp,n2dens)

! ===================================================================
      if(doDMRG)then  !yma
        call dmrg_dim_change_mclr(LRras2(1:8),nna,0)
        call dmrg_dim_change_mclr(LRras2(1:8),ntash,0)
        ! mma_deallocate later
        deallocate(tmpDe)
        deallocate(tmpDeM)
        deallocate(tmpP)
        deallocate(tmpPM)
      end if
! ===================================================================

      return
      end

#ifdef _NOTUSED_
      Subroutine Projecter(CID,CI,D,P)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 CI(*),CID(*),P(*),D(*)
      Call Getmem('De','allo','real',ipde,ntash**2)
      Call Getmem('P','allo','real',ipp,n2dens)
      Do i=0,nroots-1
       Do j=0,nroots-1
        r=ddot_(nconf1,ci(1+nconf1*i),1,
     &                   cid(1+nconf1*j),1)
        Call CIDENS2(ci(1+nconf1*i),
     &               ci(1+nconf1*j),State_sym,State_sym,
     *     work(ipp),work(ipde))
       call daxpy_(ntash**2,-r*weight(1+i),work(ipde),1,d,1)
       call daxpy_(n2dens,-r*weight(1+i),work(ipp),1,p,1)
       End Do
      End Do
      Call Getmem('De','free','real',ipde,ntash**2)
      Call Getmem('P','free','real',ipp,n2dens)
      Return
      End
#endif
