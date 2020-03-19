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
* Copyright (C) 2007, Bjorn O. Roos                                    *
************************************************************************
      SUBROUTINE AFREEZ(NSYM,NBAS,NFRO,NISH,NASH,NSSH,NDEL,NAME,
     &           NAMFRO,LNFRO,DPQ,THRFR,THRDE,IFQCAN,CMO,NCMO)
*****************************************************************************
*                                                                           *
* Purpose: to select orbitals, which will be frozen in the CASPT2           *
* calculations based on a selection of atoms controlled by the input        *
* keyword ATOMS.                                                            *
* Each inactive orbital is checked for the fraction of electrons located    *
* on the selected atoms. If smaller than a given threshold, the orbital     *
* will be frozen.                                                           *
* Called by READIN_CASPT2                                                   *
* Author: B. O. Roos in July 2007 for MOLCAS-7                              *
*     Calling parameters:                                                   *
*     NSYM   : Number of symmetries                                         *
*     NFRO   : Number of frozen orbitals (modified by the program)          *
*     NISH   : Number of inactive orbitals                                  *
*     Name   : Center and function type label per basis function            *
*     Namfro : names of atoms to be selected (length lnfro)                 *
*     Labfro : labels for orbitals to be frozen                             *
*     CMO    : Orbital coefficients                                         *
*     OccN   : Orbital occupations                                          *
*     SMat   : Overlap matrix                                               *
*     DPQ    : The charge matrix for a given orbital                        *
*     THRFR : Threshold for freezing orbitals                               *
*     THRDE : Threshold for deleting orbitals                               *
*                                                                           *
*****************************************************************************
      USE REFWFN
      IMPLICIT REAL*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
*
      CHARACTER(4) NAME(2,*),NAMFRO(*)
      DIMENSION NBAS(NSYM),NFRO(NSYM),NISH(NSYM),NASH(NSYM),NSSH(NSYM),
     &          NDEL(NSYM)
      DIMENSION LABFRO(mxbas),DPQ(*)
      REAL*8, ALLOCATABLE :: SMAT(:)
      REAL*8 CMO(*)
*
*
*----------------------------------------------------------------------*
*     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
*----------------------------------------------------------------------*
*
*      Write(6,*) 'Entering AFreez'
      NBAST=0
      ntri=0
      Do I=1,NSYM
        NBAST=NBAST+NBAS(I)
        ntri=(nbas(i)+nbas(i)**2)/2+ntri
      End Do
      IF(NBAST.GT.MXBAS) then
       Write(6,'(/6X,A)')
     & 'The number of basis functions exceeds the present limit'
       Call Abend
      Endif
*
*----------------------------------------------------------------------*
*     Read the overlap matrix                                          *
*----------------------------------------------------------------------*
      NSMAT=NTRI+6
      CALL MMA_ALLOCATE(SMAT,NSMAT)
      isymlbl=1
      Call RdOne(irc,6,'Mltpl  0',1,SMAT,isymlbl)
*
*----------------------------------------------------------------------*
*      write(6,*)'molecular orbitals before localization'
*      imo=0
*      do isym=1,nsym
*       nbi=nbas(isym)
*       do ib=1,nbi
*        write(6,*) 'orbital', isym, ib
*        write(6,'(4E18.12)') (CMO(imo+i),i=1,nbi)
*       imo=imo+nbi
*       enddo
*      enddo
*----------------------------------------------------------------------*
*     Localize the inactive and virtual orbitals                       *
*----------------------------------------------------------------------*
      Thrs=1.d-06
      Call Cho_x_Loc(irc,Thrs,nSym,nBas,nFro,nIsh,nAsh,nSsh,CMO)
      If(irc.ne.0) then
       write(6,*) 'Localization failed. The AFRE option cannot be used'
       Call Abend
      Endif
*      write(6,*)'molecular orbitals after localization'
*      imo=0
*      do isym=1,nsym
*       nbi=nbas(isym)
*       do ib=1,nbi
*        write(6,*) 'orbital', isym, ib
*        write(6,'(4E18.12)') (CMO(imo+i),i=1,nbi)
*       imo=imo+nbi
*       enddo
*      enddo
*----------------------------------------------------------------------*
*     Compute Mulliken atomic charges for each center and              *
*     each orbital.                                                    *
*----------------------------------------------------------------------*
*
      nb2=0
      Do isym=1,nsym
       nb2=nb2+nbas(isym)*(nbas(isym)+1)/2
      Enddo
*      write(6,*) 'Starting the calculation',nb2
      Do i=1,nb2
       DPQ(i)=0.0d0
      Enddo
      ib=0
      imo0=0
      ipq0=0
      Do isym=1,nsym
       nbi=nbas(isym)
       nfi=nfro(isym)
       nin=nish(isym)
       imo=imo0+nbi*nfi
       If(nin.ne.0) then
        Do i=1,nin
         labfro(i)=0
        Enddo
         Do ni=1,nin
*         write(6,*) 'loop over sym and inactive orbitals',isym,ni
          ipq=ipq0
          ipq1=0
          Do np=1,nbi
           Do nq=1,np
            ipq=ipq+1
            ipq1=ipq1+1
            DPQ(ipq1)=
     &      CMO(imo+np)*CMO(imo+nq)*SMAT(ipq)
           Enddo
          Enddo
*         DPQ is the charge matrix for orbital ni in symmetry isym
*         Now add non-diagonal elements to the diagonal
          ipq1=0
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
           iqq=0
           Do nq=1,np
            iqq=iqq+nq
            ipq1=ipq1+1
            If(np.ne.nq) then
             DPQ(ipp)=DPQ(ipp)+DPQ(ipq1)
             DPQ(iqq)=DPQ(iqq)+DPQ(ipq1)
            Endif
           Enddo
          Enddo
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
*          write(6,*) 'diagonal element',ipp,DPQ(ipp)
          Enddo

*         The diagonal now contains the charges for each basis function
*         Add charges for basis functions centered on the selected atoms
*         First check that the sum is equal to one
          chksum=0.0d0
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
           chksum=chksum+DPQ(ipp)
          Enddo
          If(abs(chksum-1.d0).gt.1.d-08) then
           Write(6,*) 'Error on Checksum in Afreez.',
     &     'Value is not equal to 1:', isym, ni, chksum
           Write(6,*) 'Freezing extra orbitals in CASPT2 stops.'
           Call Abend
          Endif
*         Add diagonal elements that belong to selected atoms
          selch=zero
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
           Do iname=1,lnfro
            if(name(1,ib+np).eq.namfro(iname)) selch=selch+DPQ(ipp)
           Enddo
          Enddo
          If(abs(selch).lt.thrfr) labfro(ni)=1
*         write(6,*) selch
          imo=imo+nbi
         Enddo
*        Sort the inactive CMO's such that frozen orbitals are first.
         nfro1=nfro(isym)
         Do  ni=1,nin
          If(labfro(ni).eq.1) then
*          Exchange this orbital with the first inactive orbital
           ist1=nfro(isym)*nbi+imo0
           ist2=(nfro1+ni-1)*nbi+imo0
*          write(6,*)'nfro,nish',nfro(isym),nish(isym),ist1,ist2
           Do np=1,nbi
            Swap=CMO(ist1+np)
            CMO(ist1+np)=CMO(ist2+np)
            CMO(ist2+np)=Swap
           Enddo
           nfro(isym)=nfro(isym)+1
           nish(isym)=nish(isym)-1
          Endif
         Enddo
        Endif
        ipq0=ipq0+nbi*(nbi+1)/2
        imo0=imo0+nbi**2
        ib=ib+nbi
      Enddo
*     Now sort virtual orbitals
*     Orbitals with too low population on selected atoms will be deleted
      Do i=1,nb2
       DPQ(i)=zero
      Enddo
      ib=0
      imo0=0
      ipq0=0
      Do isym=1,nsym
       nbi=nbas(isym)
       ndi=ndel(isym)
       nsi=nssh(isym)
       nssh(isym)=0
       ndel(isym)=ndi+nsi
       imo=imo0+nbi*(nfro(isym)+nish(isym)+nash(isym))
       If(nsi.ne.0) then
        Do i=1,nsi
         labfro(i)=0
        Enddo
         Do ni=1,nsi
*         write(6,*) 'loop over sym and secondary orbitals',isym,ni
          ipq=ipq0
          ipq1=0
          Do np=1,nbi
           Do nq=1,np
            ipq=ipq+1
            ipq1=ipq1+1
            DPQ(ipq1)=
     &      CMO(imo+np)*CMO(imo+nq)*SMAT(ipq)
           Enddo
          Enddo
*         DPQ is the charge matrix for orbital ni in symmetry isym
*         Now add non-diagonal elements to the diagonal
          ipq1=0
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
           iqq=0
           Do nq=1,np
            iqq=iqq+nq
            ipq1=ipq1+1
            If(np.ne.nq) then
             DPQ(ipp)=DPQ(ipp)+DPQ(ipq1)
             DPQ(iqq)=DPQ(iqq)+DPQ(ipq1)
            Endif
           Enddo
          Enddo
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
*          write(6,*) 'diagonal element',ipp,DPQ(ipp)
          Enddo

*         The diagonal now contains the charges for each basis function
*         Add charges for basis functions centered on the selected atoms
*         First check that the sum is equal to one
          chksum=0.0d0
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
           chksum=chksum+DPQ(ipp)
          Enddo
          If(abs(chksum-1.d0).gt.1.d-08) then
           Write(6,*) 'Error on Checksum in Afreez.',
     &     'Value is not equal to 1:', isym, ni, chksum
           Write(6,*) 'Deleting extra orbitals in CASPT2 stops.'
           Call Abend
          Endif
*         Write(6,*) 'Checksum', isym, ni, chksum
*         Add diagonal elements that belong to selected atoms
          selch=zero
          ipp=0
          Do np=1,nbi
           ipp=ipp+np
           Do iname=1,lnfro
            if(name(1,ib+np).eq.namfro(iname)) selch=selch+DPQ(ipp)
           Enddo
          Enddo
          If(abs(selch).gt.thrde) labfro(ni)=1
*         write(6,*) selch
          imo=imo+nbi
         Enddo
*        Sort the CMO's such that secondary orbitals are first.
         Do  ni=1,nsi
          If(labfro(ni).eq.1) then
*          Exchange this orbital with the first deleted orbital
           ist1=(nfro(isym)+nish(isym)+nash(isym)+nssh(isym))*nbi+imo0
           ist2=(nfro(isym)+nish(isym)+nash(isym)+ni-1)*nbi+imo0
           Do np=1,nbi
            Swap=CMO(ist1+np)
            CMO(ist1+np)=CMO(ist2+np)
            CMO(ist2+np)=Swap
           Enddo
*          write(6,*)'Orbital number',ni,ist1,ist2
*           write(6,'(4E18.12)') (CMO(ist1+np),np=1,nbi)
*           write(6,'(4E18.12)') (CMO(ist2+np),np=1,nbi)
           ndel(isym)=ndel(isym)-1
           nssh(isym)=nssh(isym)+1
          Endif
         Enddo
        Endif
        ipq0=ipq0+nbi*(nbi+1)/2
        imo0=imo0+nbi**2
        ib=ib+nbi
      Enddo

*      imo=0
*      do isym=1,nsym
*       nbi=nbas(isym)
*       do ib=1,nbi
*        write(6,*) 'orbital', isym, ib
*        write(6,'(4E18.12)') (CMO(imo+i),i=1,nbi)
*       imo=imo+nbi
*       enddo
*      enddo
*     Write the resorted MO's back to JobIph
*
      IF (IFQCAN.NE.0) IFQCAN=0 ! MOs to be recanonicalized on exit

      CALL MMA_DEALLOCATE(SMAT)
      Return
*
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(NCMO)
      End
