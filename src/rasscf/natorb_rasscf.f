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
      Subroutine NatOrb_RASSCF(CMOO,SCR1,SCR2,SMAT,CMON,OCCN)
C
C     RASSCF program: version IBM-3090: Output section
C
C     PURPOSE: Calculation of natural orbitals from the
C              density matrix. These orbitals are written onto JOBIPH
C              in position IADR15(12), followed by the occupation
C              numbers in IADR15(13).
C              The calculation is performed for each of the NROOT
C              density matrices obtained in an average CASSCF calc.
C              Called from MAIN before OUTCTL
C
C          ****** IBM 3090 MOLCAS Release: 90 02 22 ******
C

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
#include "splitcas.fh"

      Parameter (ROUTINE='NATORB  ')

      DIMENSION CMOO(*),SCR1(*),SCR2(*),SMAT(*),CMON(*),OCCN(*)

      IPRLEV=IPRLOC(7)
      iDisk=IADR15(12)
      jDisk=IADR15(3)

*        Write(LF,*)
*        Write(LF,*) ' CMO in NATORB_RASSCF very beginning'
*        Write(LF,*) ' ---------------------'
*        Write(LF,*)
*        ioff=0
*        Do iSym = 1,nSym
*         iBas = nBas(iSym)
*         if(iBas.ne.0) then
*           write(6,*) 'Sym =', iSym
*           do i= 1,iBas
*             write(6,*) (CMOO(ioff+iBas*(i-1)+j),j=1,iBas)
*           end do
*           iOff = iOff + (iBas*iBas)
*         end if
*        End Do

      if(.not.DoSplitCAS) then
        Do kRoot = 1,lRoots
          If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
            Write(LF,*)
            Write(LF,'(6X,A,I3)')
     &          'Natural orbitals and occupation numbers for root',kRoot
          End If
          Call DDaFile(JOBIPH,2,SCR1,NACPAR,jDisk)
          Call DDaFile(JOBIPH,0,SCR1,NACPAR,jDisk)
          Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
          Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
          Call DBLOCK(SCR1)

          Call dCopy_(NTOT,(0.0d0),0,OCCN,1)
          Call dCopy_(NTOT2,CMOO,1,CMON,1)

          ID=0
          ISTMO1=0
          IB=0
          DO ISYM=1,NSYM
            NB=NBAS(ISYM)
            NFI=NFRO(ISYM)+NISH(ISYM)
            NAO=NASH(ISYM)
            NB2=NB**2
            NA1=ITRI(NAO+1)
            IO=IB+NFI
            ISTMO=ISTMO1+NB*NFI
C
C  set occupation number of frozen and inactive orbitals
C
            Call dCopy_(NFI,(2.0d0),0,OCCN(IB+1),1)
C
C  Diagonalize the density matrix and transform orbitals
C
            IF(NAO.GT.0) THEN
              Call dCopy_(NAO*NAO,(0.0d0),0,SCR2,1)
              Call dCopy_(NAO,(1.0d0),0,SCR2,NAO+1)
              CALL JACOB(SCR1(ID+1),SCR2,NAO,NAO)
              II=0
              DO I=1,NAO
                II=II+I
                OCCN(IO+I)=SCR1(II+ID)
              END DO
              IST=IO+1
              IEND=IO+NAO
              If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
                Write(LF,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')
     &                'sym',iSym,':',(OCCN(I),I=IST,IEND)
              End If
              CALL DGEMM_('N','N',
     &                    NB,NAO,NAO,
     &                    1.0d0,CMOO(ISTMO+1),NB,
     &                    SCR2,NAO,
     &                    0.0d0,CMON(ISTMO+1),NB)
            END IF
C
            ID=ID+NA1
            ISTMO1=ISTMO1+NB2
            IB=IB+NB
          END DO
C
C ORTHOGONALIZE NEW MO'S
C
          CALL SUPSCH(SMAT,CMOO,CMON)
          CALL ORTHO_RASSCF(SMAT,SCR1,CMON,SCR2)
C
C Place NOs wrt decreasing order of OCCN
c***  GLM commented it off for this reordering is not consistent
cGLM throughout Molcas
C
cGLM          iOff=0
cGLM          jOff=0
cGLM          Do iSym=1,nSym
cGLM            NB=nBas(iSym)
cGLM            NAO=nAsh(iSym)
cGLM            IA=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
cGLM            JA=1+jOff+nFro(iSym)+nIsh(iSym)
cGLM            Call Order_Arrays('decr',CMON(IA),NB,NAO,OCCN(JA),SCR1)
cGLM            iOff=iOff+NB**2
cGLM            jOff=jOff+NB
cGLM          End Do

       If(IPRLEV.ge.DEBUG) Then
           Write(LF,*)
           Write(LF,*) ' CMON in NATORB_RASSCF after ORDER_ARRAYS'
           Write(LF,*) ' ---------------------'
           Write(LF,*)
           ioff=0
           Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
              end do
              iOff = iOff + (iBas*iBas)
            end if
           End Do
       end if
C
C Write new molecular orbitals and occupation numbers to JOBIPH
C
C
          CALL DDAFILE(JOBIPH,1,CMON,NTOT2,iDisk)
          CALL DDAFILE(JOBIPH,1,OCCN,NTOT,iDisk)
        End Do

      else   ! if DoSplitCAS (GLMJ)...
        If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
          Write(LF,*)
          Write(LF,'(6X,A,I3)')
     &     'Natural orbitals and occupation numbers for root',lRootSplit
        End If
        Call DDaFile(JOBIPH,2,SCR1,NACPAR,jDisk)
        Call DDaFile(JOBIPH,0,SCR1,NACPAR,jDisk)
        Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
        Call DDaFile(JOBIPH,0,SCR1,NACPR2,jDisk)
        Call DBLOCK(SCR1)

        Call dCopy_(NTOT,(0.0d0),0,OCCN,1)
        Call dCopy_(NTOT2,CMOO,1,CMON,1)

        ID=0
        ISTMO1=0
        IB=0
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          NFI=NFRO(ISYM)+NISH(ISYM)
          NAO=NASH(ISYM)
          NB2=NB**2
          NA1=ITRI(NAO+1)
          IO=IB+NFI
          ISTMO=ISTMO1+NB*NFI
C
C  set occupation number of frozen and inactive orbitals
C
          Call dCopy_(NFI,(2.0d0),0,OCCN(IB+1),1)
C
C  Diagonalize the density matrix and transform orbitals
C
          IF(NAO.GT.0) THEN
            Call dCopy_(NAO*NAO,(0.0d0),0,SCR2,1)
            Call dCopy_(NAO,(1.0d0),0,SCR2,NAO+1)
            CALL JACOB(SCR1(ID+1),SCR2,NAO,NAO)
            II=0
            DO I=1,NAO
              II=II+I
              OCCN(IO+I)=SCR1(II+ID)
            END DO
            IST=IO+1
            IEND=IO+NAO
            If(KSDFT.eq.'SCF'.and.IPRLEV.ge.USUAL) Then
               Write(LF,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')
     &               'sym',iSym,':',(OCCN(I),I=IST,IEND)
            End If
            CALL DGEMM_('N','N',
     &                  NB,NAO,NAO,
     &                  1.0d0,CMOO(ISTMO+1),NB,
     &                  SCR2,NAO,
     &                  0.0d0,CMON(ISTMO+1),NB)
          END IF
C
         ID=ID+NA1
         ISTMO1=ISTMO1+NB2
         IB=IB+NB
       END DO
C
C ORTHOGONALIZE NEW MO'S
C
       CALL SUPSCH(SMAT,CMOO,CMON)
       CALL ORTHO_RASSCF(SMAT,SCR1,CMON,SCR2)
C
C Place NOs wrt decreasing order of OCCN
C
       iOff=0
       jOff=0
       Do iSym=1,nSym
          NB=nBas(iSym)
          NAO=nAsh(iSym)
          IA=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
          JA=1+jOff+nFro(iSym)+nIsh(iSym)
          Call Order_Arrays('decr',CMON(IA),NB,NAO,OCCN(JA),SCR1)
          iOff=iOff+NB**2
          jOff=jOff+NB
       End Do
C
C Write new molecular orbitals and occupation numbers to JOBIPH
C
       CALL DDAFILE(JOBIPH,1,CMON,NTOT2,iDisk)
       CALL DDAFILE(JOBIPH,1,OCCN,NTOT,iDisk)

      end if

      RETURN
      END
