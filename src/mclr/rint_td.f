************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Jonna Stalring                                         *
************************************************************************
      SubRoutine RInt_td(ekappa,mkappa,isym)
c
c
*******************************************************
*            [2]                                      *
* Modifies  E    *kappa by subtracting Omega*S*kappa  *
* to make timedep PT possible.                        *
*                                                     *
*******************************************************
c
c     mkappa       is d/dx(kappa) in matrix form (derivative of orb rot mat w r t PT)
c     ekappa       is E*d/dx(kappa) in matrix form, REPLACED BY
c                  the entire timedep expression
c     isym         Rinttd is called once for each symmtry
c
c Local variables
c
c     dens The density matrix
c     wDKt   Omega*(density matrix)*(kappa transposed)
c     wKtD   As above but different order
c
c
      Implicit Real*8 (a-h,o-z)
c
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "Input.fh"
#include "real.fh"
#include "Pointers.fh"
      Real*8 ekappa(ndens2),mkappa(ndens2)
      Real*8, Allocatable:: Dens(:), wDKt(:), wKtD(:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
c
c---------------------------------------------------------
c
c  Allocate memory
c
      lDens = 0
      Do iS=1,nSym
         lDens = lDens + nBas(iS)**2
      End Do
      Call mma_allocate(Dens,lDens,Label='Dens')
      Call mma_allocate(wDKt,ndens2,Label='wDKt')
      Call mma_allocate(wKtD,ndens2,Label='wKtD')
c
c---------------------------------------------------------
c***************************************
c Construct the density matrix
c***************************************
c
c For HF
c
      Dens(:)=Zero
      ip3 = 1
      Do iS=1,nSym
          inc = nBas(iS) + 1
          call dcopy_(nIsh(iS),[2.0d0],0,Dens(ip3),inc)
          ip3 = ip3 + nBas(iS)*nBas(iS)
      End Do
*
*         do is=1,nsym
*            size=nbas(is)*nbas(is)      ! The size of the sym i and j martix
*            do k=-1,(size-1)             ! Loop over the entire submatrix
*               residual=DMOD((k+1),(nbas(is)+1))
*               if ((residual.eq.zero).and.(k.lt.(nbas(is)*nIsh
*     &      (is))).and.(nIsh(is).ne.0)) then
*                     Dens(1+ipCM(is)+k)=two
*               else
*                     Dens(1+ipCM(is)+k)=zero
*               end if
*            end do
*          end do
c
c For a CASSCF wavefunc. From Anders subrut r2elint
c Add the active active dens
c
          Do iS=1,nSym
             Do iB=1,nAsh(iS)
                Do jB=1,nAsh(iS)
                   ip=ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
                   iA=nA(is)+ib
                   jA=nA(is)+jb
                   ip2=itri(iA,jA)
                   Dens(ip)=Work(ipG1t+ip2-1)
                End Do
             End Do
          End Do
c*******************************************
c Multiply D and mkappa wDKt and wKtD
c*******************************************
*      Call RecPrt('dens ',' ',dens,ldens,1)
*      Call RecPrt('ekappa ',' ',ekappa,ndens2,1)
      do is=1,nsym
           js=iEOr(is-1,isym-1)+1
C      wDKt
           If ( (nBas(iS).gt.0).and.(nBas(jS).gt.0) ) Then
             call DGEMM_('n','n',nbas(is),nbas(js),nbas(is),
     &                   Two*Omega,Dens(ipCM(is)),nbas(is),
     &                             mkappa(ipmat(is,js)),nbas(is),
     &                   zero,wDKt(ipmat(is,js)),nbas(is))
C      wKtD
             call DGEMM_('n','n',nbas(is),nbas(js),nbas(js),
     &                   Two*Omega,mkappa(ipmat(is,js)),nbas(is),
     &                             Dens(ipCM(js)),nbas(js),
     &                   zero,wKtD(ipmat(is,js)),nbas(is))

c*****************************************************
c            Replace ekappa ekappa=ekappa-wDKt+wKtD
c*****************************************************
             incx=1
             incy=incx
             lenght=nbas(is)*nbas(js)
             call daxpy_(lenght,1.0d0,wDKt(ipmat(is,js)),
     &            incx,ekappa(ipmat(is,js)),incy)
             call daxpy_(lenght,-1.0d0,wKtD(ipmat(is,js)),
     &            incx,ekappa(ipmat(is,js)),incy)
          End If
      end do
*      Call RecPrt('wDKt ',' ',wDKt,ndens2,1)
*      Call RecPrt('wKtD ',' ',wKtD,ndens2,1)
*      Call RecPrt('ekappa ',' ',ekappa,ndens2,1)
c
c  Free memory
c
       Call mma_deallocate(Dens)
       Call mma_deallocate(wDKt)
       Call mma_deallocate(wKtD)
      return
      end
