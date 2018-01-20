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
      SUBROUTINE PRCEVC(NSS,FRAC,SOENE,MAPST,MAPSP,MAPMS,UMATR,UMATI)
      IMPLICIT NONE
#include "Molcas.fh"
#include "cntrl.fh"
#include "prgm.fh"
#include "stdalloc.fh"
#include "rassiwfn.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='PRCEVC')

      INTEGER NSS
      REAL*8 FRAC,SOENE(NSS),UMATR(NSS,NSS),UMATI(NSS,NSS)
      INTEGER MAPST(NSS),MAPSP(NSS),MAPMS(NSS)

      real*8, allocatable :: weight(:),sstate(:)
      real*8 :: wmax(5),smax(5)
      integer :: nmax(5)

      INTEGER I
      INTEGER ISFS,IMAXSTATE,ISTATE,ISS,JSS,JSTA,JEND,NW
      REAL*8 CFFLIM,S,SSMAX,SZ,TST,WGTMAX,XMAX

      call mma_allocate(weight,nss)
      call mma_allocate(sstate,nss)

C Write out the complex eigenvectors of the spin-orbit states.
C Four states at a time are written out.
C They are written as complex numbers, four on each line.
C Coefficients are written out if any of the four complex coeffs
C on the same line have absolute value at least as large as
C a certain fraction (FRAC) of the largest such value for the
C foursome of states.

      CALL QENTER(ROUTINE)

      If(IPGLOB.ge.verbose) then
      DO JSTA=1,NSS,4
       JEND=MIN(NSS,JSTA+3)
       WRITE(6,*)
       WRITE(6,'(1X,A16,F16.8,3(2X,F16.8))')
     &            '    Energy (au) ',(SOENE(JSS),JSS=JSTA,JEND)
       WRITE(6,'(1X,A16,6X,I4,3(14X,I4))')
     &            ' SFS  S     Ms  ',(JSS,JSS=JSTA,JEND)
* Scan coefficients to pick out the largest:
       WGTMAX=0.0D0
       DO ISS=1,NSS
        DO JSS=JSTA,JEND
         WGTMAX=MAX(WGTMAX,UMATR(ISS,JSS)**2+UMATI(ISS,JSS)**2)
        END DO
       END DO

* Scan coefficients, write if large enough:
       CFFLIM=FRAC*SQRT(WGTMAX)
       DO ISS=1,NSS
        ISTATE=MAPST(ISS)
        S=DBLE(MAPSP(ISS)-1)*0.5D0
        SZ=DBLE(MAPMS(ISS))*0.5D0
        TST=0.0D0
        DO JSS=JSTA,JEND
         TST=MAX(TST,UMATR(ISS,JSS)**2+UMATI(ISS,JSS)**2)
        END DO
        IF(TST.GE.CFFLIM**2) THEN
         WRITE(6,
     &        '(I4,1X,F4.1,1X,F5.1,3X,4(A1,F7.4,A1,F7.4,A1,1x))')
     &         ISTATE,S,SZ,('(',UMATR(ISS,JSS),','
     &                            ,UMATI(ISS,JSS),')',JSS=JSTA,JEND)
        END IF
       END DO
      END DO
      else If(IPGLOB.ge.usual) then
*     Write out the weights of the five most important spin-free states
*     For each spin-orbit state (BOR in Krapperup 070226)
*
      Write(6,*)
      Write(6,*)
      Write(6,*)'Weights of the five most important'//
     &         ' spin-orbit-free states for each spin-orbit state.'
      Write(6,*)
      Write(6,*) 'SO State  Total energy (au)           Spin-free '//
     & 'states, spin, and weights'
      Write(6,*)'----------------------------------------------------'//
     &          '---------------------------------------------------'
      Do iss=1,nss
       Do isfs=1,nstate
        weight(isfs)=0.d0
       Enddo
       Do jss=1,nss
        istate=mapst(jss)
        sstate(istate)=dble(mapsp(jss)-1)*0.5D0
        weight(istate)=weight(istate)+
     &  umatr(jss,iss)**2+umati(jss,iss)**2
       Enddo
*      Sort the weights
       imaxstate=0
  100  Continue
        xmax=0.d0
        ssmax=0.d0
        nw=0
       Do istate=1,nstate
        if(weight(istate).ge.xmax) then
         nw=istate
         xmax=weight(istate)
         ssmax=sstate(istate)
        Endif
       Enddo
       weight(nw)=-1.d0
       imaxstate=imaxstate+1
       nmax(imaxstate)=nw
       wmax(imaxstate)=xmax
       smax(imaxstate)=ssmax
       if(imaxstate.lt.min(nstate,5)) go to 100
*
       write(6,'(i5,1x,f16.6,3x,5(i5,f4.1,f8.4))')
     &  iss,soene(iss),(nmax(i),smax(i),wmax(i),i=1,(min(nstate,5)))

      Enddo
      Write(6,*)'----------------------------------------------------'//
     &          '---------------------------------------------------'
      Endif

C Added by Ungur Liviu on 04.11.2009
C Addition of UMATR and UMATI on RunFile

      CALL Put_dArray('UMATR_SINGLE',UMATR,NSS**2)
      CALL Put_dArray('UMATI_SINGLE',UMATI,NSS**2)

      call mma_deallocate(weight)
      call mma_deallocate(sstate)

      CALL QEXIT(ROUTINE)

      RETURN
      END
