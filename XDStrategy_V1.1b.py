#!/usr/bin/python

import re, sys
from copy import deepcopy
from collections import OrderedDict
standard = [
'01 SCALE',
'02 CC M',
'03 CC D Q O H',
'04 CC M D Q O H',
'05 CC U2',
'06 CC M D Q O H',
'07 CC M D Q O H U2',
'08 CC XYZ',
'09 CC M D Q O H XYZ',
'10 CC M D Q O H XYZ U2',
'11 CC KAPPA[]',
'12 CC M',
'13 CC M KAPPA[]',
'14 CC M D Q O H XYZ U2',
'15 CC KAPPA[]',
'16 CC M D Q O H XYZ U2 KAPPA[]',
'17 CC HXYZ',
'18 CC M D Q O H XYZ U2',
'19 CC D Q O H XYZ U2 KAPPA[]',
'20 CC M D Q O H XYZ U2 KAPPA[]',
'21 CC KAPPAP[]',
'22 CC M D Q O H XYZ U2 KAPPA[]',
'23 CC SIGOBS[2] M D Q O H XYZ U2 KAPPA[]',
'24 CC SIGOBS[1] M D Q O H XYZ U2 KAPPA[]',
'25 CC SIGOBS[0] M D Q O H XYZ U2 KAPPA[]'
]

#
#
#
atom_table_regex = re.compile('(?P<ATOM>[a-zA-Z\(\)0-9]+)\s+(?P<ATOM0>[a-zA-Z\(\)0-9]+)\s+(?P<AX1>[XYZxyz])\s+(?P<ATOM1>[a-zA-Z\(\)0-9]+)\s+(?P<ATOM2>[a-zA-Z\(\)0-9]+)\s+(?P<AX2>[XYZxyz])\s+(?P<RL>[RLrl])\s+(?P<TP>\d)\s+(?P<TBL>\d+)\s+(?P<KAP>\d+)\s+(?P<LMX>\d)\s+(?P<SITESYM>[0-9a-zA-Z_]*){0,1}(\s+(?P<CHEMCON>.*)\s*)*\n')
key_table_regex = re.compile('(?P<ATOM>[a-zA-Z\(\)0-9]+)\s+(?P<XYZ>[01]{3})\s+(?P<U2>[01]{6})\s+(?P<U3>[01]{10})\s+(?P<U4>[01]{15})\s+(?P<M>[01]{2})\s+(?P<D>[01]{3})\s+(?P<Q>[01]{5})\s+(?P<O>[01]{7})\s+(?P<H>[01]{9}).*\n')
scat_table_regex = re.compile('(?P<SCAT>[a-zA-Z\+\-]+)\s+(?P<CORE>\w+)\s+(?P<SPHV>\w+)\s+(?P<DEFV>\w+)\s+(?P<S1>[0-9\-]+)\s+(?P<S2>[0-9\-]+)\s+(?P<S3>[0-9\-]+)\s+(?P<S4>[0-9\-]+)\s+(?P<P2>[0-9\-]+)\s+(?P<P3>[0-9\-]+)\s+(?P<P4>[0-9\-]+)\s+(?P<D3>[0-9\-]+)\s+(?P<D4>[0-9\-]+)\s+(?P<F4>[0-9\-]+)\s+(?P<S5>[0-9\-]+)\s+(?P<P5>[0-9\-]+)\s+(?P<S6>[0-9\-]+)\s+(?P<P6>[0-9\-]+)\s+(?P<D5>[0-9\-]+)\s+(?P<S7>[0-9\-]+)\s+(?P<D6>[0-9\-]+)\s+(?P<F5>[0-9\-]+)\s+(?P<DELF>[0-9\.\-]+)\s+(?P<DELFf>[0-9\.\-]+)\s+(?P<NSCTL>[0-9\.\-]+)\s*')
kappa_regex = re.compile('KAPPA\s+[01]+\s*')
four_regex = re.compile('FOUR\s+fmod1\s+4\s+[23]\s+0\s+0\s+fmod2\s+-1\s+[23]\s+0\s+0')
model_regex = re.compile('SELECT\s+\**model\s+4\s+[23]\s+1\s+0\s+based_on\s+F\^2\s+test\s+verbose\s+1\s+')
skip_regex = re.compile('SKIP\s+[\*]{0,1}obs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sigobs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sinthl\s+[d\d.]+\s+[d\d.]+\s+')
vibcon_regex = re.compile('CON(\s+[\d.-]+\s+[U\d/]+)+\s+=\s+0\s+.*')
reset_regex = re.compile('RESET\s+BOND(\s+[a-zA-Z\(\)0-9]+){2}\s+[\d.]+\s+')
cycles_regex = re.compile('SELECT\s+cycle\s+[0-9\-]+\s+dampk\s+[0-9\.]+\s+cmin\s+[0-9\.]+\s+cmax\s+[0-9\.]+\s+eigcut\s+[d0-9\.\-]+\s+convcrit\s+[d0-9\.\-]+\s+')
#
#
#
def writefile(number, OutMaster):
  outfile = open('xd{}.mas'.format(number),'w')
  for i in OutMaster:
    outfile.write('{}\n'.format(i.strip()))
  outfile.close()
  
#
#
#
def getinstructions(line):
  if line.startswith('#'):
    return None
  elif line == '':
    return None
  else:
    instruction = line.split()[1:]
    number = line.split()[0]
    dictionary = {}
    for i in instruction:
      re.sub('[\[\];]',' ', i).split()
      dictionary[re.sub('[\[\];]',' ', i).split()[0]] = re.sub('[\[\];]',' ', i).split()[1:]
    return number, dictionary
#
#
#
def getatomsandlines(MasterFile, atomclassdict):
  kappalist = []
  fourline = 0
  modelline = 0
  skipline = 0
  cyclesline = 0
  scatdict = OrderedDict()
  for i in range(len(MasterFile)):
    if atom_table_regex.search(MasterFile[i]):
      atom = atom_table_regex.match(MasterFile[i]).groupdict()
      atomclassdict[atom['ATOM']] = Atom(scatdict.keys()[int(atom['TBL'])-1],atom['KAP'], atom['CHEMCON'])
      atomclassdict[atom['ATOM']].atomtable = MasterFile[i]
      atomclassdict[atom['ATOM']].atomline = i  
    elif key_table_regex.search(MasterFile[i]):
      key = key_table_regex.match(MasterFile[i]).groupdict()
      atomclassdict[key['ATOM']].dict = key
      atomclassdict[key['ATOM']].keyline = i
    elif scat_table_regex.search(MasterFile[i]):
      scat = scat_table_regex.match(MasterFile[i]).groupdict()
      #
      # Falls der scattering factor schon existiert (C, Cv fehler proudly presented by C. Schuermann)
      #
      if scat['SCAT'] in scatdict:
        scatdict[scat['SCAT']+'dublette'] = scat
      else: scatdict[scat['SCAT']] = scat
    elif kappa_regex.match(MasterFile[i]): kappalist.append(i)
    elif four_regex.match(MasterFile[i]): fourline = i
    elif skip_regex.match(MasterFile[i]): skipline = i
    elif model_regex.match(MasterFile[i]): modelline = i
    elif cycles_regex.match(MasterFile[i]): cyclesline = i
    elif MasterFile[i] == 'END ATOM': break
  return kappalist, fourline, skipline, modelline, cyclesline
#
#
#
def removeshit(MasterFile):
  removecons = []
  removeresets = []
  for i in range(len(MasterFile)):
    if re.match('^\!{0,1}CON(\s+[\d.-]+\s+[U\d/]+)+\s+=\s+0\s+.*', MasterFile[i]): removecons.append(i)
  for i in reversed(removecons):
    MasterFile.pop(i)
  for i in range(len(MasterFile)):
    if re.match('^\!{0,1}RESET\s+', MasterFile[i]): removeresets.append(i)
  for i in reversed(removeresets):
    MasterFile.pop(i)
#
#
#
def putinstructionsinaction(MasterFile, resets, vibcons, atomclassdict, line, number, instruction):
  kappalist, fourline, skipline, modelline, cyclesline = getatomsandlines(MasterFile, atomclassdict)
  OutMaster = deepcopy(MasterFile)
  OutMaster[skipline] = re.sub('\*', '', OutMaster[skipline])
  for i in kappalist:
    OutMaster[i] = re.sub('1', '0', OutMaster[i])
  if 'KAPPA' in instruction or 'KAPPAP' in instruction:
    if 'KAPPA' in instruction and 'KAPPAP' in instruction:
      both = [x for x in instruction['KAPPA'] if x in instruction['KAPPAP']]
      for i in both:
        instruction['KAPPA'].remove(i)
        instruction['KAPPAP'].remove(i)
        for i in instruction['KAPPA']:
          OutMaster[kappalist[int(i)-1]] = re.sub('[0,1]{2}','10',OutMaster[kappalist[int(i)-1]], count = 1)
        for i in instruction['KAPPAP']:
          OutMaster[kappalist[int(i)-1]] = re.sub('[0,1]{2}','01',OutMaster[kappalist[int(i)-1]], count = 1)
        for i in both:
          OutMaster[kappalist[int(i)-1]] = re.sub('[0,1]{2}','11',OutMaster[kappalist[int(i)-1]], count = 1)
    elif 'KAPPA' in instruction:
      for i in instruction['KAPPA']:
        OutMaster[kappalist[int(i)-1]] = re.sub('[0,1]{2}','10',OutMaster[kappalist[int(i)-1]], count = 1)
    elif 'KAPPAP' in instruction:
      for i in instruction['KAPPAP']:
        OutMaster[kappalist[int(i)-1]] = re.sub('[0,1]{2}','01',OutMaster[kappalist[int(i)-1]], count = 1)
  if 'HXYZ' not in instruction or 'SIGOBS' not in instruction or 'SINTHL' not in instruction:
    OutMaster[skipline] =  re.sub('sigobs', '*sigobs', OutMaster[skipline]) 
  if 'HXYZ' in instruction:
    OutMaster[skipline] =  re.sub('sinthl\s+[d\d.]+\s+[d\d.]+', '*sinthl 0.0 0.5', OutMaster[skipline])
  if 'SIGOBS' in instruction:
    OutMaster[skipline] =  re.sub('sigobs\s+[d\d.]+', '*sigobs '+ instruction['SIGOBS'][0], OutMaster[skipline])
  if 'SINTHL' in instruction:
    OutMaster[skipline] =  re.sub('sinthl\s+[d\d.]+\s+[d\d.]+', '*sinthl {} {}'.format(instruction['SINTHL'][0], instruction['SINTHL'][1]), OutMaster[skipline])
  if 'SCALE' in instruction:
    OutMaster[cyclesline] = re.sub('cycle\s+-*','cycle -',OutMaster[cyclesline])        
  for i in atomclassdict:
    for j in ['U2', 'XYZ', 'HXYZ']:
      try:
        if i not in instruction[j] and atomclassdict[i].atomtype == 'H' and j == 'HXYZ':
          atomclassdict[i].instruction.append(j)
        else:
          atomclassdict[i].instruction.append(j)
      except: pass 
    for j in ['M', 'D', 'Q', 'O', 'H', 'NOSYMM', 'CC']:
      try:
        if i not in instruction[j]:
          atomclassdict[i].instruction.append(j)
      except: pass
    for j in ['U3', 'U4']:
      try:
        if i in instruction[j]:
          if j == 'U3': 
            atomclassdict[i].U3 = True
            atomclassdict[i].instruction.append(j)
          else: 
            atomclassdict[i].U4 = True
            atomclassdict[i].instruction.append(j)
      except: pass
      if atomclassdict[i].U3: 
        OutMaster[modelline] = re.sub('model\s+4\s+[23]\s+1\s+0', 'model 4 3 1 0',OutMaster[modelline])
        OutMaster[fourline] = re.sub('fmod2\s+-1\s+[23]\s+0\s+0', 'fmod2 -1 3 0 0',re.sub('fmod1\s+4\s+[23]\s+0\s+0\s+', 'fmod1 4 3 0 0 ', OutMaster[fourline]))
      if atomclassdict[i].U4: 
        OutMaster[modelline] = re.sub('model\s+4\s+[23]\s+1\s+0', 'model 4 4 1 0',OutMaster[modelline])
        OutMaster[fourline] = re.sub('fmod2\s+-1\s+[23]\s+0\s+0', 'fmod2 -1 4 0 0',re.sub('fmod1\s+4\s+[23]\s+0\s+0\s+', 'fmod1 4 4 0 0 ', OutMaster[fourline]))
    OutMaster[atomclassdict[i].keyline] =  atomclassdict[i].print_key(atomclassdict)
    OutMaster[atomclassdict[i].atomline] = atomclassdict[i].printatomtable()
  atomclassdict[i].instruction=[]
  if 'XYZ' in instruction:
    for i in range(1,len(resets)+1):
      OutMaster.insert(fourline + i, resets[i-1])
  if 'U2' in instruction:
    for i in range(1,len(vibcons)+1):
      OutMaster.insert(fourline + i, vibcons[i-1])
  return OutMaster
#  except TypeError: pass
#
# jedes atom im masterfile ist eine instanz der klasse Atom.
# die klasse hat eine funktion print_key die aus der klassenvariable 'instruction' vom typ liste 
# die fuer das atom spezifische zeile der keytable.
#
class Atom():
  def __init__(self, atomtype, kappa, chemcon = None):
    self.dict = {}
    self.atomtable = ''
    self.keyline = 0
    self.atomline = 0
    self.chemcon = chemcon.strip()
    self.atomtype = atomtype
    self.U3 = False
    self.U4 = False
    self.instruction = []
    self.order = ['XYZ', 'HXYZ', 'U2', 'U3', 'U4', 'M', 'D', 'Q', 'O', 'H']
  def print_key(self, atomclassdict):
    keydict = {}
    keydict['ATOM'] = self.dict['ATOM']
    if self.chemcon:
      for i in self.order:
        if i == 'HXYZ': pass
        elif i in ['U3', 'U2', 'U4', 'XYZ']:
          keydict[i] = self.dict[i]
        else: 
          if 'NOSYMM' in self.instruction:
            keydict[i] = re.sub('0','1', atomclassdict[self.chemcon].dict[i])
          else: keydict[i] = self.dict[i]
    else:
      if 'NOSYMM' in self.instruction:
        for i in self.order:
          if i == 'HXYZ': pass
          elif i in ['M', 'D', 'Q', 'O', 'H']:
           keydict[i] = re.sub('0','1',self.dict[i])
          else: keydict[i] = self.dict[i]
      else:
        keydict = self.dict
    out = '{:<7}'.format(keydict['ATOM'])
    for i in self.order:
      if i in self.instruction:
        if i == 'HXYZ' and self.atomtype != 'H':
          self.instruction.remove('HXYZ')
        if i == 'HXYZ' and self.atomtype == 'H':
          out = '{} {}'.format(out, keydict['XYZ'])
          continue
        elif i == 'XYZ' and self.atomtype != 'H':
          out = '{} {}'.format(out, keydict['XYZ'])
          continue
        elif i == 'XYZ' and self.atomtype == 'H' and 'HXYZ' not in self.instruction: 
          out = '{} {}'.format(out, re.sub('1', '0', keydict[i]))
          continue
        elif i in ['M', 'D', 'Q', 'O', 'H'] and 'CC' not in self.instruction:
          out = '{} {}'.format(out, keydict[i])
        elif 'CC' in self.instruction and self.chemcon and i in ['M', 'D', 'Q', 'O', 'H']:
          out = '{} {}'.format(out, re.sub('1', '0', keydict[i]))
        else: 
          if i == 'HXYZ': pass
          else: out = '{} {}'.format(out, keydict[i])
      else:
        if i == 'HXYZ': pass
        elif i == 'XYZ' and 'HXYZ' in self.instruction and self.atomtype == 'H': pass
        else: 
          out = '{} {}'.format(out, re.sub('1','0', keydict[i]))
    return out
  def printatomtable(self):
    outstring = self.atomtable
    if 'CC' not in self.instruction and self.chemcon: outstring = ''.join(outstring.rsplit(self.chemcon,1))
    if self.U3: outstring = re.sub('\s+[124]', '   3',outstring, count = 1)
    elif self.U4: outstring = re.sub('\s+[124]', '   3',outstring, count = 1)
    return outstring

def main():
  filename = ''
  while not filename:
    filename = raw_input('Enter name of source XD master file:[xd.mas] ') or 'xd.mas'
    try:
      MasterFile = open(filename,'r').readlines()
    except:
      print 'Source XD master file not found.'
  filename = ''
  while not filename:
    filename = raw_input('Enter name of file that contains suitable XD constraints for hydrogen atoms:[xd.const] ') or 'xd.const'
    try:
      ConstFile = open(filename,'r').readlines()
    except:
      print 'XD constraints file not found.'
      print 'Please execute XDConstraints.py'
      sys.exit()
  filename = ''
  while not filename:
    filename = raw_input('Enter name of file that contains instructions for strategy generation:[xd.const] ') or 'xd.inst'
    try:
      InstFile = open(filename,'r')
    except:
      print 'XD constraints file not found.'
      print 'Writing instructions file with standard strategy'
      TempFile = open('xd.inst', 'w') 
      for i in standard:
        TempFile.write(i)
      TempFile.close()
      print 'File written! Please edit instructions and restart this program!'
      sys.exit
  removeshit(MasterFile)
  atomclassdict = OrderedDict()
  vibcons = []
  resets = []
  for i in ConstFile:
    if vibcon_regex.match(i):
      vibcons.append(i.strip())
    if reset_regex.match(i):
      resets.append(i.strip())
  for line in InstFile:
    try:
      number, instruction = getinstructions(line.strip())
      OutMaster = putinstructionsinaction(MasterFile, resets, vibcons, atomclassdict, line, number, instruction)
      writefile(number, OutMaster)
    except: pass

if __name__ == "__main__":
    main()
