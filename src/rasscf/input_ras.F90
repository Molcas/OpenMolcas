!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module Input_RAS
use definitions, only: iwp
Private
#include "rasdim.fh"
! Logical unit number for reading input:
      Integer(kind=iwp), Public:: LuInput
! Used for input processing
      Integer(kind=iwp), Parameter, Public :: NKeys=150
!------------------------------------------------------
! Logical flags, to check whether a keyword has been used
! in the input:
      Logical(kind=iwp), Public:: KeyFlags(0:NKeys)
      Logical(kind=iwp), Public:: KeyAAAA,                              &
     &     KeyALTE, KeyATOM,          KeyCHAR, KeyCHOI,                 &
     &     KeyCHOL, KeyCIMX, KeyCION, KeyCIRE, KeyCIRO,                 &
     &     KeyCISE, KeyCLEA, KeyCORE, KeyDELE, KeyEND,                  &
     &     KeyFILE, KeyFROZ, KeyHOME, KeyINAC, KeyINPO,                 &
     &     KeyIPHN, KeyITER, KeyJOBI, KeyKSDF, KeyLEVS,                 &
     &     KeyLINE, KeyLOWD, KeyLOWM, KeyLUMO, KeyMAXO,                 &
     &     KeyNACT, KeyNEWI, KeyNONE, KeyNOQU, KeyOPTO,                 &
     &     KeyORBA, KeyORBL, KeyORBO, KeyORDE, KeyOUTO,                 &
     &     KeyOUTP, KeyPRIN, KeyPROR, KeyPRSP, KeyPRWF,                 &
     &     KeyQUNE, KeyRAS1, KeyRAS2, KeyRAS3, KeyGASS,                 &
     &     KeyRFPE, KeyCIRF, KeyRFRO, KeyRLXR, KeyRASS,                 &
     &     KeySDAV, KeySPIN, KeySUPS, KeySXDA, KeySYMM,                 &
     &     KeyTHRS, KeyTIGH, KeyTITL, KeyTYPE, KeyVB  ,                 &
     &     KeyEXPE, KeySPLI, KeyNUSP, KeyENSP, KeyPESP,                 &
     &     KeyFOSP, KeyMDRL, KeyOFEM, KeyFTHA, KeyDFMD,                 &
     &     KeyBKAP, KeyALPH, KeyFARO, KeyDMRG, Key3RDM,                 &
     &     KeyNECI, KeyTOTA, KeyTIME, KeyNMCY, KeyCALC,                 &
     &     KeyRDMS, KeyREAL, KeyDEFI, KeyDIAG, KeyEMBD,                 &
     &     KeyBLOK, KeySOCC, KeyRGIN, KeyPRSD, KeyFCID,                 &
     &     KeyNOCA, KeySAVE, KeyEXPA, KeyH5OR, KeyH5CI,                 &
     &     KeyHEXS, KeyHEUR, KeyDMPO, KeyNEVP, KeyHFOC,                 &
     &     KeyDAVT, KeyCHRE, KeyCHBL, KeyMXSW, KeyNOIS,                 &
     &     KeyDMRE, KeyMXCA, KeyDEXS, KeyHROO, KeyTDM,                  &
     &     KeyDFCF, KeyNKEE, KeyREOR, KeyTRIA, KeyPOPS,                 &
     &     KeySEMI, KeyMEMO, KeyIVO , KeyCRPR, KeyRDML,                 &
     &     KeyORTH, KeyCCCI, KeyROST, KeyXMSI, KeyCMSI,                 &
     &     KeyCMMA, KeyCMMI, KeyCMTH, KeyGUGA, KeyCMSS,                 &
     &     KeyCMSO, KeyPERI, KeySSCR, KeyMCM7, KeyWRMA,                 &
     &     KeyDICE, KeySTOC, KeyEPSI, KeySAMP, KeyDITE,                 &
     &     KeyDIRE, KeyDIOC, KeyPPT2, KeyNDPT, KeyRGRA,                 &
     &     KeySTAV

! we need the common block for the equivalence statement below to make sense.
      Common /InputFlags/ KeyAAAA,                                      &
     &     KeyALTE, KeyATOM,          KeyCHAR, KeyCHOI,                 &
     &     KeyCHOL, KeyCIMX, KeyCION, KeyCIRE, KeyCIRO,                 &
     &     KeyCISE, KeyCLEA, KeyCORE, KeyDELE, KeyEND,                  &
     &     KeyFILE, KeyFROZ, KeyHOME, KeyINAC, KeyINPO,                 &
     &     KeyIPHN, KeyITER, KeyJOBI, KeyKSDF, KeyLEVS,                 &
     &     KeyLINE, KeyLOWD, KeyLOWM, KeyLUMO, KeyMAXO,                 &
     &     KeyNACT, KeyNEWI, KeyNONE, KeyNOQU, KeyOPTO,                 &
     &     KeyORBA, KeyORBL, KeyORBO, KeyORDE, KeyOUTO,                 &
     &     KeyOUTP, KeyPRIN, KeyPROR, KeyPRSP, KeyPRWF,                 &
     &     KeyQUNE, KeyRAS1, KeyRAS2, KeyRAS3, KeyGASS,                 &
     &     KeyRFPE, KeyCIRF, KeyRFRO, KeyRLXR, KeyRASS,                 &
     &     KeySDAV, KeySPIN, KeySUPS, KeySXDA, KeySYMM,                 &
     &     KeyTHRS, KeyTIGH, KeyTITL, KeyTYPE, KeyVB  ,                 &
     &     KeyEXPE, KeySPLI, KeyNUSP, KeyENSP, KeyPESP,                 &
     &     KeyFOSP, KeyMDRL, KeyOFEM, KeyFTHA, KeyDFMD,                 &
     &     KeyBKAP, KeyALPH, KeyFARO, KeyDMRG, Key3RDM,                 &
     &     KeyNECI, KeyTOTA, KeyTIME, KeyNMCY, KeyCALC,                 &
     &     KeyRDMS, KeyREAL, KeyDEFI, KeyDIAG, KeyEMBD,                 &
     &     KeyBLOK, KeySOCC, KeyRGIN, KeyPRSD, KeyFCID,                 &
     &     KeyNOCA, KeySAVE, KeyEXPA, KeyH5OR, KeyH5CI,                 &
     &     KeyHEXS, KeyHEUR, KeyDMPO, KeyNEVP, KeyHFOC,                 &
     &     KeyDAVT, KeyCHRE, KeyCHBL, KeyMXSW, KeyNOIS,                 &
     &     KeyDMRE, KeyMXCA, KeyDEXS, KeyHROO, KeyTDM,                  &
     &     KeyDFCF, KeyNKEE, KeyREOR, KeyTRIA, KeyPOPS,                 &
     &     KeySEMI, KeyMEMO, KeyIVO , KeyCRPR, KeyRDML,                 &
     &     KeyORTH, KeyCCCI, KeyROST, KeyXMSI, KeyCMSI,                 &
     &     KeyCMMA, KeyCMMI, KeyCMTH, KeyGUGA, KeyCMSS,                 &
     &     KeyCMSO, KeyPERI, KeySSCR, KeyMCM7, KeyWRMA,                 &
     &     KeyDICE, KeySTOC, KeyEPSI, KeySAMP, KeyDITE,                 &
     &     KeyDIRE, KeyDIOC, KeyPPT2, KeyNDPT, KeyRGRA,                 &
     &     KeySTAV

      Equivalence(KeyAAAA,KeyFlags(0))
!------------------------------------------------------
! Actual keywords, note: Order matters!
      Character(LEN=4), Parameter, Public::                             &
     &    CMD(nKeys) = ['ALTE','ATOM',       'CHAR','CHOI',             &
     &                  'CHOL','CIMX','CION','CIRE','CIRO',             &
     &                  'CISE','CLEA','CORE','DELE','END ',             &
     &                  'FILE','FROZ','HOME','INAC','INPO',             &
     &                  'IPHN','ITER','JOBI','KSDF','LEVS',             &
     &                  'LINE','LOWD','LOWM','LUMO','MAXO',             &
     &                  'NACT','NEWI','NONE','NOQU','OPTO',             &
     &                  'ORBA','ORBL','ORBO','ORDE','OUTO',             &
     &                  'OUTP','PRIN','PROR','PRSP','PRWF',             &
     &                  'QUNE','RAS1','RAS2','RAS3','GASS',             &
     &                  'RFPE','CIRF','RFRO','RLXR','RASS',             &
     &                  'SDAV','SPIN','SUPS','SXDA','SYMM',             &
     &                  'THRS','TIGH','TITL','TYPE','VB  ',             &
     &                  'EXPE','SPLI','NUSP','ENSP','PESP',             &
     &                  'FOSP','MDRL','OFEM','FTHA','DFMD',             &
     &                  'BKAP','ALPH','FARO','DMRG','3RDM',             &
     &                  'NECI','TOTA','TIME','NMCY','CALC',             &
     &                  'RDMS','REAL','DEFI','DIAG','EMBD',             &
     &                  'BLOK','SOCC','RGIN','PRSD','FCID',             &
     &                  'NOCA','SAVE','EXPA','H5OR','H5CI',             &
     &                  'HEXS','HEUR','DMPO','NEVP','HFOC',             &
     &                  'DAVT','CHRE','CHBL','MXSW','NOIS',             &
     &                  'DMRE','MXCA','DEXS','HROO','TDM ',             &
     &                  'DFCF','NKEE','REOR','TRIA','POPS',             &
     &                  'SEMI','MEMO','IVO ','CRPR','RDML',             &
     &                  'ORTH','CCCI','ROST','XMSI','CMSI',             &
     &                  'CMMA','CMMI','CMTH','GUGA','CMSS',             &
     &                  'CMSO','PERI','SSCR','MCM7','WRMA',             &
     &                  'DICE','STOC','EPSI','SAMP','DITE',             &
     &                  'DIRE','DIOC','PPT2','NDPT','RGRA',             &
     &                  'STAV']
ENd Module Input_RAS
