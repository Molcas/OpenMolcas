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
      Implicit real*8 (a-h,o-z)
      Dimension H1EFF(*)
* ----------------------------------------------------------------
* Purpose: Reads and adds one-electron naked Hamiltonian into H1EFF.
* Dress it with reaction field (if any).
* Also get POTNUC at the same time.
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "ofembed.fh"
*
* NOT TESTED (used for OFEmbed below)
#if 0
      Character(Len=16) NamRfil
#endif
      Logical Found

c Add naked one-el Hamiltonian in AO basis to H1EFF:
      CALL GETMEM('ONEHAM','ALLO','REAL',LONEHAM,NBTRI)
      IRC=-1
      IOPT=6
      ICOMP=1
      ISYLBL=1
      CALL RDONE(IRC,IOPT,'OneHam  ',ICOMP,WORK(LONEHAM),ISYLBL)
      CALL DAXPY_(NBTRI,1.0D0,WORK(LONEHAM),1,H1EFF,1)
      CALL GETMEM('ONEHAM','FREE','REAL',LONEHAM,NBTRI)

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
         Call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
         Call Get_dScalar('RF Self Energy',ERFSelf)
         Call Get_dArray('Reaction field',Work(lTemp),nTemp)
         If (Found) Call NameRun('RUNFILE')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nTemp,1.0D0,Work(lTemp),1,H1EFF,1)
         Call GetMem('RFFLD','Free','Real',lTemp,nTemp)
      End If

#ifdef _DEBUGPRINT_
         WRITE(6,*)' 1-EL HAMILTONIAN (MAY INCLUDE REACTION FIELD)'
         ISTLT=0
         DO ISYM=1,NSYM
           IF ( NBAS(ISYM).GT.0 ) THEN
             WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
             CALL TRIPRT(' ',' ',H1EFF(1+ISTLT),NBAS(ISYM))
             ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
           END IF
         END DO
#endif

* NOT TESTED
#if 0
c If this is a perturbative Orbital-Free Embedding (OFE) calculation
c then modify the one-electron Hamiltonian by the OFE potential and
c the nuclear attraction by the Rep_EN
      If ( Done_OFEmb ) then
         nTemp=0
         Do iSym=1,nSym
            nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Call GetMem('DCoul','Allo','Real',ipCoul,nTemp)
         Call FZero(Work(ipCoul),nTemp)
         If (First_OFE) Then
            Call GetMem('FMaux','Allo','Real',ipFMaux,nTemp)
            Call Coul_DMB(.true.,1,Rep_EN,Work(ipFMaux),Work(ipCoul),
     &                             Work(ipCoul),nTemp)
         EndIf
         Call DaXpY_(nTemp,1.0d0,Work(ipFMaux),1,H1EFF,1)
         Call GetMem('DCoul','Free','Real',ipCoul,nTemp) ! used as Dum
         First_OFE=.false.
*
         Call Get_NameRun(NamRfil) ! save the old RUNFILE name
         Call NameRun('AUXRFIL')   ! switch the RUNFILE name
         Call Get_dExcdRa(ipVxc,nVxc)
         Call DaXpY_(nTemp,1.0d0,Work(ipVxc),1,H1EFF,1)
         If (nVxc.eq.2*nTemp) Then ! but fix for Nuc Attr added twice
            Call DaXpY_(nTemp,1.0d0,Work(ipVxc+nTemp),1,H1EFF,1)
            Call Get_dArray('Nuc Potential',Work(ipVxc),nTemp)
            Call DaXpY_(nTemp,-1.0d0,Work(ipVxc),1,H1EFF,1)
         EndIf
         Call Free_Work(ipVxc)
         Call NameRun(NamRfil)   ! switch back to old RUNFILE
#ifdef _DEBUGPRINT_
             WRITE(6,*)' 1-EL HAMILTONIAN INCLUDING OFE POTENTIAL'
             ISTLT=0
             DO ISYM=1,NSYM
               IF ( NBAS(ISYM).GT.0 ) THEN
                 WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
                 CALL TRIPRT(' ',' ',H1EFF(1+ISTLT),NBAS(ISYM))
                 ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
               END IF
             END DO
#endif
      End If
#endif

*
      RETURN
      END
