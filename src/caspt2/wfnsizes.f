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
      subroutine wfnsizes
************************************************************************
*
* Compute various orbital sizes
*
************************************************************************
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

      Integer NASHT2
      Integer NI, NR1, NR2, NR3, NS, N123
      Integer I, iSym

* Table sizes
      Integer IIABS, ITABS, IAABS
      Integer II, IA, IS, IO
      Integer ITOT, IINA, IACT, IEXT

      NFROT=0
      NISHT=0
      NASHT=0
      NRAS1T=0
      NRAS2T=0
      NRAS3T=0
      NOSHT=0
      NSSHT=0
      NDELT=0
      NORBT=0
      NBAST=0
      NOSQT=0
      NBSQT=0
      NIMX=0
      NAMX=0
      NSMX=0
      NOMX=0
      NBMX=0
      DO ISYM=1,NSYM
        NIES(ISYM)=NISHT
        NAES(ISYM)=NASHT
        NSES(ISYM)=NSSHT
        NOSH(ISYM)=NISH(ISYM)+NASH(ISYM)
        NSSH(ISYM)=NBAS(ISYM)-NFRO(ISYM)-NOSH(ISYM)-NDEL(ISYM)
        NORB(ISYM)=NOSH(ISYM)+NSSH(ISYM)
        NORBT=NORBT+NORB(ISYM)
        NOSQT=NOSQT+NORB(ISYM)**2
        NBSQT=NBSQT+NBAS(ISYM)**2
        NFROT=NFROT+NFRO(ISYM)
        NISHT=NISHT+NISH(ISYM)
        NASHT=NASHT+NASH(ISYM)
        NOSHT=NOSHT+NOSH(ISYM)
        NRAS1T=NRAS1T+NRAS1(ISYM)
        NRAS2T=NRAS2T+NRAS2(ISYM)
        NRAS3T=NRAS3T+NRAS3(ISYM)
        NSSHT=NSSHT+NSSH(ISYM)
        NDELT=NDELT+NDEL(ISYM)
        NBAST=NBAST+NBAS(ISYM)
        NIMX=MAX(NIMX,NISH(ISYM))
        NAMX=MAX(NAMX,NASH(ISYM))
        NSMX=MAX(NSMX,NSSH(ISYM))
        NOMX=MAX(NOMX,NORB(ISYM))
        NBMX=MAX(NBMX,NBAS(ISYM))
      END DO
C Set RHS Boxes to maximum size
      NINABX=NIMX
      NSECBX=NSMX
      NBTRI=(NBSQT+NBAST)/2
      NOTRI=(NOSQT+NORBT)/2
* Size of orbital transformation arrays:
      NTORB=0
      NTAT=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        N123=NR1**2+NR2**2+NR3**2
        NTAT=NTAT+N123
        NTORB=NTORB+NI**2+N123+NS**2
      END DO
C Sizes of DREF and PREF arrays:
      NDREF=1
      NPREF=1
C Sizes of GAMMA1, GAMMA2, and GAMMA3:
      NG1=1
      NG2=1
      NG3TOT=1
      IF(NASHT.GT.0) THEN
        NDREF=(NASHT**2+NASHT)/2
        NASHT2=NASHT**2
        NPREF=(NASHT2**2+NASHT2)/2
        NG1=NASHT**2
        NG2=NG1**2
        NG3TOT=((NG1+2)*(NG1+1)*NG1)/6
      END IF

C  Identify the wave function type
      ISCF=0
      IF((ISPIN.EQ.NACTEL+1).AND.(NACTEL.EQ.NASHT)) ISCF=2
      IF(NASHT.EQ.0) ISCF=1
      IF(NACTEL.EQ.2*NASHT) ISCF=1

************************************************************************
*
* Create orbital name vector
*
************************************************************************
      II=0
      IA=0
      IS=0
      ITOT=0
      IINA=0
      IACT=0
      IEXT=0
      DO ISYM=1,NSYM
        IO=0
        DO I=1,NFRO(ISYM)
          ITOT=ITOT+1
          IO=IO+1
          WRITE(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)')
     &      'Fr',ISYM,'.',IO
        END DO
        DO I=1,NISH(ISYM)
          ITOT=ITOT+1
          IINA=IINA+1
          IINAIS(IINA)=ITOT
          IO=IO+1
          WRITE(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)')
     &      'In',ISYM,'.',IO
          II=II+1
          IINAM(II)=ORBNAM(ITOT)
        END DO
        DO I=1,NASH(ISYM)
          ITOT=ITOT+1
          IACT=IACT+1
          IACTIS(IACT)=ITOT
          IO=IO+1
          WRITE(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)')
     &      'Ac',ISYM,'.',IO
          IA=IA+1
          IANAM(IA)=ORBNAM(ITOT)
        END DO
        DO I=1,NSSH(ISYM)
          ITOT=ITOT+1
          IEXT=IEXT+1
          IEXTIS(IEXT)=ITOT
          IO=IO+1
          WRITE(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)')
     &      'Se',ISYM,'.',IO
          IS=IS+1
          ISNAM(IS)=ORBNAM(ITOT)
        END DO
        DO I=1,NDEL(ISYM)
          ITOT=ITOT+1
          IO=IO+1
          WRITE(ORBNAM(ITOT),'(A2,I1,A1,I3.3,1X)')
     &      'De',ISYM,'.',IO
        END DO
      END DO

************************************************************************
*
* Precompute table sizes
*
************************************************************************
      IIABS=0
      ITABS=0
      IAABS=0
      DO ISYM=1,NSYM
        DO I=1,NORB(ISYM)
          IF(I.LE.NISH(ISYM)) THEN
            IIABS=IIABS+1
            IISYM(IIABS)=ISYM
          ELSE IF(I.LE.NISH(ISYM)+NASH(ISYM)) THEN
            ITABS=ITABS+1
            IASYM(ITABS)=ISYM
          ELSE
            IAABS=IAABS+1
            IESYM(IAABS)=ISYM
          END IF
        END DO
      END DO

*---  Check consistency of the orbitals
      If ( NISHT.gt.MXINA ) Then
        Call WarningMessage(2,'Too many inactive orbitals.')
        WRITE(6,'(a,2i8)')' NISHT >  MXINA:',NISHT,MXINA
        Call Quit_OnUserError
      End If
      If ( NASHT.gt.MXACT ) Then
        Call WarningMessage(2,'Too many active orbitals.')
        WRITE(6,'(a,2i8)')' NASHT > MXACT:',NASHT,MXACT
        Call Quit_OnUserError
      End If
      If ( NSSHT.gt.MXEXT ) Then
        Call WarningMessage(2,'Too many secondary orbitals.')
        WRITE(6,'(a,2i8)')' NSSHT > MXEXT:',NSSHT,MXEXT
        Call Quit_OnUserError
      End If
      If ( NBAST.gt.MXORB ) Then
        Call WarningMessage(2,'Too many basis functions.')
        WRITE(6,'(a,2i8)')' NBAST > MXORB:',NBAST,MXORB
        Call Quit_OnUserError
      End If

*
* GG-Nov04  The following informations must be passed to the Cholesky
* transformation section through RunFile. COMMON blocks cannot be
* used due to several conflicts.
      Call Put_iScalar('nSym',nSym)
      Call Put_iArray('nFroPT',nFro,nSym)
      Call Put_iArray('nIsh'  ,nIsh,nSym)
      Call Put_iArray('nAsh'  ,nAsh,nSym)
      Call Put_iArray('nDelPT',nDel,nSym)
      Call Put_iArray('nBas'  ,nBas,nSym)

      RETURN
      END
