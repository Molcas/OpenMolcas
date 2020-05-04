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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Readin_vb_m()
************************************************************************
*                                                                      *
*     Read the input                                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "warnings.fh"
#include "rasscf.fh"
#include "general.fh"
#include "jobiph_j.fh"
#include "gas.fh"
*
      Character*80       Blank
*---  set INVEC -> get MOs from JOBIPH file ----------------------------*
      INVEC=3
*---  Initialize -------------------------------------------------------*
      Do i=1,80
         Blank(i:i)=' '
      End Do
         Do j = 1, 18
            Title(j) = ' '
         End Do
*
*---  Read input from standard input ----------------------------------*
*---  process TITLE    command ----------------------------------------*
      ntit=0
      do ii=1,mxtit
        if(len_trim_cvb(title_j(ii)).gt.0)then
          ntit=ntit+1
          title(ntit)=title_j(ii)
        endif
      end do
*---  process NACT command --------------------------------------------*
      nactel=nactel_j
      nhole1=nhole1_j
      nelec3=nelec3_j
*---  process SPIN command --------------------------------------------*
      ispin=ispin_j
*---  process SYMM command --------------------------------------------*
      lsym=lsym_j
*---  process FROZ command --------------------------------------------*
      call imove_cvb(nfro_j,nfro,mxsym)
*---  process INAC command --------------------------------------------*
      call imove_cvb(nish_j,nish,mxsym)
*---  process RAS1 command --------------------------------------------*
      call imove_cvb(nrs1_j,nrs1,mxsym)
*---  process RAS2 command --------------------------------------------*
      call imove_cvb(nrs2_j,nrs2,mxsym)
*---  process RAS3 command --------------------------------------------*
      call imove_cvb(nrs3_j,nrs3,mxsym)
*---  process DELE command --------------------------------------------*
      call imove_cvb(ndel_j,ndel,mxsym)
*---
*TRS
*      If (nroots.gt.1.and.irlxroot.eq.0) iRlxRoot=iroot(nroots)
*      If (nroots.eq.1) iRlxRoot=0
*---  complete orbital specifications ---------------------------------*
      Do iSym=1,mxsym
         NASH(ISYM)=NRS1(ISYM)+NRS2(ISYM)+NRS3(ISYM)
         NORB(ISYM)=NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
         NSSH(ISYM)=NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
      End Do

* SVC: convert CAS/RAS to general GAS description here, then we only
* need to copy it for lucia later, which always uses GAS description.
      NGSSH(1,1:NSYM)=NRS1(1:NSYM)
      NGSSH(2,1:NSYM)=NRS2(1:NSYM)
      NGSSH(3,1:NSYM)=NRS3(1:NSYM)
      IGSOCCX(1,1) = MAX(2*SUM(NRS1(1:NSYM))-NHOLE1,0)
      IGSOCCX(1,2) = 2*SUM(NRS1(1:NSYM))
      IGSOCCX(2,1) = NACTEL - NELEC3
      IGSOCCX(2,2) = NACTEL
      IGSOCCX(3,1) = NACTEL
      IGSOCCX(3,2) = NACTEL
*
*---  Compute IZROT. IZROT is a matrix (lower triangular over the -----*
*     active space), which specifies which t,u rotations should be
*     avoided, since the orbitals belong to the same RAS space.
*     This is the only way the RAS concept is explicitly used in the
*     SX section of the program.
      ITU=0
      DO ISYM=1,mxsym
        NAO=NASH(ISYM)
        IF(NAO.GT.1) THEN
          NRAS1=NRS1(ISYM)
          NRS12=NRS2(ISYM)+NRAS1
          DO NT=2,NAO
            DO NU=1,NT-1
              ITU=ITU+1
              IZROT(ITU)=0
CSVC: check if NU<NT are included in the same gas space
              NGSSH_LO=0
              DO IGAS=1,NGAS
                NGSSH_HI=NGSSH_LO+NGSSH(IGAS,ISYM)
                IF (NU.GT.NGSSH_LO.AND.NT.LE.NGSSH_HI) THEN
                  IZROT(ITU)=1
                END IF
                NGSSH_LO=NGSSH_HI
              END DO
            END DO
          END DO
        END IF
      END DO
*---  complete the input processing -----------------------------------*
      NTOT=0
      NTOT1=0
      NTOT2=0
      NO2M=0
      NIN=0
      NAC=0
      NDELT=0
      NFROT=0
      NSEC=0
      NORBT=0
      NTOT3=0
      NTOTSP=0
      NTOT4=0
      NRS1T=0
      NRS2T=0
      NRS3T=0
      DO ISYM=1,NSYM
         NTOT=NTOT+NBAS(ISYM)
         NTOT1=NTOT1+(NBAS(ISYM)*(NBAS(ISYM)+1))/2
         NTOT2=NTOT2+NBAS(ISYM)**2
         NO2M=MAX(NO2M,NBAS(ISYM)**2)
         NRS1T=NRS1T+NRS1(ISYM)
         NRS2T=NRS2T+NRS2(ISYM)
         NRS3T=NRS3T+NRS3(ISYM)
         NFROT=NFROT+NFRO(ISYM)
         NIN=NIN+NISH(ISYM)
         NAC=NAC+NASH(ISYM)
         NDELT=NDELT+NDEL(ISYM)
         NSEC=NSEC+NSSH(ISYM)
         NORBT=NORBT+NORB(ISYM)
         NTOT3=NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
         NTOTSP=NTOTSP+(NASH(ISYM)*(NASH(ISYM)+1)/2)
         NTOT4=NTOT4+NORB(ISYM)**2
      END DO
      NACPAR=(NAC+NAC**2)/2
      NACPR2=(NACPAR+NACPAR**2)/2
*
      Call Put_iArray('nIsh',nIsh,nSym)
      Call Put_iArray('nAsh',nAsh,nSym)

*---  Exit, end -------------------------------------------------------*
      Return
      End
      Subroutine RdPAM_m(Line,iNumber,rNumber)
      Implicit Real*8 (A-H,O-Z)
#include "output_ras.fh"
      Parameter (ROUTINE='RDPAM')
*
*     Subroutine read character line and
*     read one integer and later one real
*
      Character*72 Line
      iStart=1
      Do While ( Line(iStart:iStart).eq.' ')
         iStart=iStart+1
      End Do
      iEnd=iStart
      Do While ( Line(iEnd:iEnd)    .ne.' ')
         iEnd=iEnd+1
      End Do
c
      Read(Line(iStart:iEnd-1),*,err=10,end=20) iNumber
c
      iStart=iEnd
      Do While ( Line(iStart:iStart).eq.' ')
         iStart=iStart+1
      End Do
c
      Read(Line(iStart:72),*,err=10,end=20) rNumber
c
      Return
c
 10   continue
      Write(LF,*) 'RdPAM: I/O error while reading input file'
      call quit(_RC_INPUT_ERROR_)
 20   continue
      Write(LF,*) 'RdPAM: end of file while reading input file'
      call quit(_RC_INPUT_ERROR_)
      End
