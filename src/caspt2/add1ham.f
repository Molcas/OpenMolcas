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
      SUBROUTINE ADD1HAM(H1EFF)
* NOT TESTED (used for OFEmbed below)
!#define _OFEmbed_
#ifdef _OFEmbed_
      use RunFile_procedures, only: Get_dExcdRa
      use OFembed, only: Do_OFemb, FMAux, OFE_First
#endif
      use definitions, only: iwp, wp
#ifdef _DEBUGPRINT_
      use definitions, only: u6
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use OneDat, only: sNoNuc, sNoOri
      use caspt2_module, only: ERFSelf, NBTRI, nSym, PotNuc, RFpert,
     &                         nBas

      Implicit None

      real(kind=wp), intent(inout):: H1EFF(*)
* ----------------------------------------------------------------
* Purpose: Reads and adds one-electron naked Hamiltonian into H1EFF.
* Dress it with reaction field (if any).
* Also get POTNUC at the same time.
*
      character(len=8) :: Label
      Logical(kind=iwp) Found
      real(kind=wp), allocatable:: ONEHAM(:), Temp(:)
      integer(kind=iwp) ICOMP, IOPT, IRC, ISYLBL, iSym, nTemp
#ifdef _DEBUGPRINT_
      integer(kind=iwp) ISTLT
#endif
#ifdef _OFEmbed_
      real(kind=wp), allocatable:: Coul(:)
#endif

c Add naked one-el Hamiltonian in AO basis to H1EFF:
      CALL mma_allocate(ONEHAM,NBTRI,Label='OneHam')
      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICOMP=1
      ISYLBL=1
      Label='OneHam'
      CALL RDONE(IRC,IOPT,Label,ICOMP,ONEHAM,ISYLBL)
      CALL DAXPY_(NBTRI,1.0D0,ONEHAM,1,H1EFF,1)
      CALL mma_deallocate(ONEHAM)

c Read nuclear repulsion energy:
      IRC=-1
      IOPT=0
      ICOMP=0
      ISYLBL=1
      Call Get_dScalar('PotNuc',PotNuc)

c If this is a perturbative reaction field calculation then
c modify the one-electron Hamiltonian by the reaction field and
c the nuclear attraction by the cavity self-energy
      If ( RFpert ) then
         nTemp=0
         Do iSym=1,nSym
            nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call mma_allocate(Temp,nTemp,Label='Temp')
         Call Get_dScalar('RF Self Energy',ERFSelf)
         Call Get_dArray('Reaction field',Temp,nTemp)
         If (Found) Call NameRun('#Pop')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nTemp,1.0D0,Temp,1,H1EFF,1)
         Call mma_deallocate(Temp)
      End If

#ifdef _DEBUGPRINT_
         WRITE(u6,*)' 1-EL HAMILTONIAN (MAY INCLUDE REACTION FIELD)'
         ISTLT=0
         DO ISYM=1,NSYM
           IF ( NBAS(ISYM).GT.0 ) THEN
             WRITE(u6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
             CALL TRIPRT(' ',' ',H1EFF(1+ISTLT),NBAS(ISYM))
             ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
           END IF
         END DO
#endif

#ifdef _OFEmbed_
c If this is a perturbative Orbital-Free Embedding (OFE) calculation
c then modify the one-electron Hamiltonian by the OFE potential and
c the nuclear attraction by the Rep_EN
      If ( Do_OFEmb ) then
         nTemp=0
         Do iSym=1,nSym
            nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Call mma_allocate(Coul,nTemp,Label='Coul')
         Coul(:)=0.0D0
         If (OFE_First) Then
            Call mma_allocate(FMaux,nTemp,Label='FMaux')
            Call Coul_DMB(.true.,1,Rep_EN,FMaux,Coul,Coul,nTemp)
         EndIf
         Call DaXpY_(nTemp,1.0d0,FMaux,1,H1EFF,1)
         Call mma_deallocate(Coul)
         OFE_First=.false.
*
         Call NameRun('AUXRFIL') ! switch the RUNFILE name
         Call Get_dExcdRa(Vxc,nVxc)
         Call DaXpY_(nTemp,1.0d0,Vxc,1,H1EFF,1)
         If (nVxc.eq.2*nTemp) Then ! but fix for Nuc Attr added twice
            Call DaXpY_(nTemp,1.0d0,Vxc(1+nTemp:2*nTemp),1,H1EFF,1)
            Call Get_dArray('Nuc Potential',Vxc,nTemp)
            Call DaXpY_(nTemp,-1.0d0,Vxc,1,H1EFF,1)
         EndIf
         Call mma_deallocate(Vxc)
         Call NameRun('#Pop')    ! switch back to old RUNFILE
#ifdef _DEBUGPRINT_
             WRITE(u6,*)' 1-EL HAMILTONIAN INCLUDING OFE POTENTIAL'
             ISTLT=0
             DO ISYM=1,NSYM
               IF ( NBAS(ISYM).GT.0 ) THEN
                 WRITE(u6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
                 CALL TRIPRT(' ',' ',H1EFF(1+ISTLT),NBAS(ISYM))
                 ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
               END IF
             END DO
#endif
      End If
#endif

      END SUBROUTINE ADD1HAM
