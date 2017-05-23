#! /usr/bin/env python
# $ this file belongs to the Molcas repository $

ATOMIC_NUMBERS={
"H": 1,
"He": 2,
"Li": 3,
"Be": 4,
"B": 5,
"C": 6,
"N": 7,
"O": 8,
"F": 9,
"Ne": 10,
"Na": 11,
"Mg": 12,
"Al": 13,
"Si": 14,
"P": 15,
"S": 16,
"Cl": 17,
"Ar": 18,
"K": 19,
"Ca": 20,
"Sc": 21,
"Ti": 22,
"V": 23,
"Cr": 24,
"Mn": 25,
"Fe": 26,
"Co": 27,
"Ni": 28,
"Cu": 29,
"Zn": 30,
"Ga": 31,
"Ge": 32,
"As": 33,
"Se": 34,
"Br": 35,
"Kr": 36,
"Rb": 37,
"Sr": 38,
"Y": 39,
"Zr": 40,
"Nb": 41,
"Mo": 42,
"Tc": 43,
"Ru": 44,
"Rh": 45,
"Pd": 46,
"Ag": 47,
"Cd": 48,
"In": 49,
"Sn": 50,
"Sb": 51,
"Te": 52,
"I": 53,
"Xe": 54,
"Cs": 55,
"Ba": 56,
"La": 57,
"Hf": 72,
"Ta": 73,
"W": 74,
"Re": 75,
"Os": 76,
"Ir": 77,
"Pt": 78,
"Au": 79,
"Hg": 80,
"Tl": 81,
"Pb": 82,
"Bi": 83,
"Po": 84,
"At": 85,
"Rn": 86,
"Fr": 87,
"Ra": 88,
"Ac": 89,
"Ku": 104,
"Ha": 105
}

class Parser:

   def __init__(self,infile_name):
      self.infile_name=infile_name
      try:
         i=open(self.infile_name,'r')
      except IOError:
         print "No such file or directory:",self.infile_name
         import sys
         sys.exit()
      input=i.readlines()
      i.close()
      self.molecules=self._parse(input)

   def _parse(self,input):
      molecules={}
      count=0
      for line in input:
         if line[0] == '#':
            if count > 0:
               molecules[name]=data
            name=line[1:].split()[0]
            if name in molecules.keys():
               print name,"used more than once"
               import sys
               sys.exit(1)
            data=[]
            count+=1
         else:
            l=line.split()
            a=l[3]
            atom=a[0].upper()+a[1:]
            d=[atom,float(l[0]),float(l[1]),float(l[2])]
            data.append(d)
      if count > 0:
         molecules[name]=data
      return molecules

   def Molecules(self):
      return self.molecules.keys()

   def NumberOfMolecules(self):
      return len(self.Molecules())

   def Atoms(self,molecule):
      return [x[0] for x in self.molecules[molecule]]

   def toXYZ(self,molecule,comment="",conversion_factor=0.52917720859):
      c=conversion_factor
      s=""
      if molecule in self.Molecules():
         atoms=self.molecules[molecule]
         s+=str(len(atoms))+"\n"
         s+=molecule+" "+comment+"\n"
         for atom in atoms:
            s+="%-6s  %20.14f  %20.14f  %20.14f \n" %(atom[0],c*atom[1],c*atom[2],c*atom[3])
      return s

   def AlltoXYZ(self,comment=""):
      s=""
      for molecule in self.Molecules():
         s+=self.toXYZ(molecule,comment=comment)
      return s

   def toSewardInput(self,molecule,basis="$Basis.....",unit="/Angstrom",conversion_factor=0.52917720859):
      c=conversion_factor
      s=""
      if molecule in self.Molecules():
         s+=" &SEWARD &END\n"
         s+="\nTitle\n"
         s+=" "+molecule+"\n\n"
         labels={}
         atoms=self.molecules[molecule]
         for atom in atoms:
            if atom[0] in labels.keys():
               labels[atom[0]]+=1
            else:
               labels[atom[0]]=1
            mylabel=atom[0]+"_"+str(labels[atom[0]])
            s+="\nBasis Set \n"
            s+=atom[0]+"."+basis+"\n"
            s+="%-6s  %20.14f  %20.14f  %20.14f %-10s \n" %(mylabel,c*atom[1],c*atom[2],c*atom[3],unit)
            s+="End Of Basis\n"
         s+="\nEnd Of Input"
      return s

   def toBlockedSewardInput(self,molecule,basis="$Basis",unit="/Angstrom",conversion_factor=0.52917720859):
      c=conversion_factor
      s=""
      if molecule in self.Molecules():
         s+=" &SEWARD &END\n"
         s+="\nTitle\n"
         s+=" "+molecule+"\n\n"
         labels={}
         atoms=self.molecules[molecule]
         for atom in atoms:
            if atom[0] in labels.keys():
               labels[atom[0]]+=1
            else:
               labels[atom[0]]=1
         for label in labels.keys():
            s+="\nBasis Set\n"
            s+=label+"."+basis+"\n"
            i=0
            for atom in atoms:
               if atom[0] == label:
                  i+=1
                  if labels[label] == 1:
                     mylabel=atom[0]
                  else:
                     mylabel=atom[0]+"_"+str(i)
                  s+="%-6s  %20.14f  %20.14f  %20.14f %-10s \n" %(mylabel,c*atom[1],c*atom[2],c*atom[3],unit)
            s+="End Of Basis\n"
         s+="\nEnd Of Input"
      return s

def SCFInput(*args):
    s="\n &SCF &END\n"
    for arg in args:
        s+=str(arg)+"\n"
    s+="End Of Input"
    return s

def addKeywordsToInput(input,*args):
    s=""
    for line in input.split("\n"):
       if line == "End Of Input":
          for arg in args:
             s+=str(arg)+"\n"
          s+="\n"
       s+=line+"\n"
    return s


if __name__=="__main__":
   import sys
   inp=sys.argv[1:]
   if len(inp)==0:
      print "No input file specified"
      sys.exit()
   # parse input file
   data=Parser(inp[0])
   # get list of molecules
   molecules=data.Molecules()
   # sort the list alphabetically (in place)
   molecules.sort()
   # for each molecule:
   for molecule in molecules:
      # get a seward input file (a string), assuming bohr as unit in input file
      input=data.toBlockedSewardInput(molecule,basis="$Basis",unit="/bohr",conversion_factor=1.0e0)
      # add CD keywords to seward input
      input=addKeywordsToInput(input,"$Seward_Options")
      # Count number of electrons in neutral system
      Ne=0
      for atom in data.Atoms(molecule):
         Ne+=ATOMIC_NUMBERS[atom]
      # add an SCF input section (set charge if needed)
      if molecule[-1] == "-":
         # add 1 to number of electrons
         Ne+=1
         # UHF for odd number of electrons, else RHF
         if Ne%2==1:
            model="UHF"
            input+=SCFInput("UHF","Charge"," -1","$SCF_Options")
         else:
            model="RHF"
            input+=SCFInput("Charge"," -1","$SCF_Options")
      else:
         # UHF for odd number of electrons, else RHF
         if Ne%2==1:
            model="UHF"
            input+=SCFInput("UHF","$SCF_Options")
         else:
            model="RHF"
            input+=SCFInput("$SCF_Options")
      # write input to file molecule.inp
      file="input/"+molecule+".input"
      f=open(file,"w")
      print >>f, input
      print molecule,"\t",file,"\t",model
      f.close()
