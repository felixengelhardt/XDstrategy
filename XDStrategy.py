# README:
# this script automatically generates xdXX.mas files according to a given strategy.
# this is work in progress and far from finished!
#
# INPUT:
#  xd.mas	- correctly setup (see section below)
#  xd.const	- constraints file (generated with XDConstraint.py!)
#  xd.inst	- instruction file (see section below)
# OUTPUT:
#  xdXX.mas files
#
# >> xd.mas:
#  will be changed by the script!   F or F^2
#                  | |              |
#  > SELECT *model 4 2 1 0 based_on F^2 test verbose 1
#  > SELECT  cycle 25 dampk 1. cmin 0.6 cmax 1. eigcut 1.d-09 convcrit 1.d-06
#                  |
#  do not enter a minus sign! will be added by the script in the SCALE step! 
#  
#	 													the CHEMCONS need to be setup correctly!
#	 																								  |
#  > ATOM     ATOM0    AX1 ATOM1    ATOM2   AX2 R/L TP  TBL KAP LMX SITESYM  CHEMCON
#  > S(1)     P(1)      Z  S(1)     C(1)     Y   R   2   1   1   4  _6            
#  > P(1)     S(1)      Z  P(1)     C(14)    X   R   2   2   2   4  _3my          
#  > C(1)     P(1)      Z  C(1)     C(2)     X   R   2   3   3   4  _mm2          
#  > C(2)     C(1)      X  C(2)     C(7)     Y   R   2   3   4   4  _mz           
#  > C(3)     C(2)      X  C(3)     C(4)     Y   R   2   3   5   4  _mz           
#  > C(4)     C(3)      X  C(4)     C(5)     Y   R   2   3   6   4  _mz           
#  > C(5)     C(4)      X  C(5)     C(6)     Y   R   2   3   6   4  _mz      C(4)
#  > C(6)     C(5)      X  C(6)     C(7)     Y   R   2   3   5   4  _mz      C(3)
#  > C(7)     C(2)      X  C(7)     C(6)     Y   R   2   3   4   4  _mz      C(2)
#  
#	 															    assume maximum symmetry!
#	 														   	|	 |		|		 |			 |
#  > KEY     XYZ --U2-- ----U3---- ------U4------- M- -D- --Q-- ---O--- ----H----
#  > S(1)    111 111111 0000000000 000000000000000 10 001 10000 1000000 100000000      
#  > P(1)    111 111111 0000000000 000000000000000 10 001 10000 1000010 100001000      
#  > C(1)    111 111111 0000000000 000000000000000 10 001 10010 1001000 100100010      
#  > C(2)    111 111111 0000000000 000000000000000 10 110 10011 0110011 100110011      
#  > C(3)    111 111111 0000000000 000000000000000 10 110 10011 0110011 100110011      
#  > C(4)    111 111111 0000000000 000000000000000 10 110 10011 0110011 100110011      
#  > C(5)    111 111111 0000000000 000000000000000 00 000 00000 0000000 000000000 !C(4)	- skip multipole
#  > C(6)    111 111111 0000000000 000000000000000 00 000 00000 0000000 000000000 !C(3)	  parameters of
#  > C(7)    111 111111 0000000000 000000000000000 00 000 00000 0000000 000000000 !C(2)	  CHEMCONed atoms
#  
#  Add Kappa constraints as demonstrated below, the script will activate them when Kappas are refined
#  > !CON  num1 par1/iat1 num2 par2/iat2 ... = num0
#  > !CON  1 KS/10 -1 KS/11 = 0
#  > !CON  1 KS/11 -1 KS/12 = 0
#  
#  enter values for '*sinthl 0. 1.22' but SIGOBS will be edited by the script!
#  > SKIP  obs 0. 1.d10 *sigobs 0 1.d06 *sinthl 0. 1.22
# 
# >> xd.inst
#  valid expressions:
#  CC				- use chemical constraints
#  SIGOBS[0]	- use sigma cutoff at [INTEGER]
#  NOSYMM[]		- release the multipol symmetry constraints for ALL (non H) atoms (NOSYMM[C(1);S(1)] will release the symmetry for ALL but C(1) and S(1)!)
#  M D Q O H	- refine Mono-, Di-, Quadro-, Okta- and/or Hexadecapoles
#  XYZ			- refine coordinates of non H atoms
#  HXYZ			- refine coordinates of H atoms (adds sintl cut at 0.5 and adds RESET BOND to the next cycle)
#  U2				- refine vibrational parameters of ALL atoms
#  U3[]			- refine 3rd order gram-charlier parameters for given atoms (syntax: U3[C(1);C(2)] refines GC3 for C(1) and C(2) ONLY)
#  U4[]			- refine 4th order gram-charlier parameters for given atoms (syntax: U4[C(1);C(2)] refines GC4 for C(1) and C(2) ONLY)
#  KAPPA[1;2;3;4;5;6;7;8;9]		- refine Kappa (syntax: KAPPA[1;2;3;7] will refine Kappas 1, 2, 3 and 7)
#  KAPPAP[1;2;3;4;5;6;7;8;9]		- refine Kappa Prime (syntax: KAPPAP[1;2;3;7] will refine Kappa Prime 1, 2, 3 and 7)
#
# >> standard strategy
standard = [
'SCALE',
'CC M',
'CC D Q O H',
'CC M D Q O H',
'CC U2',
'CC M D Q O H',
'CC M D Q O H U2',
'CC XYZ',
'CC M D Q O H XYZ',
'CC M D Q O H XYZ U2',
'CC KAPPA[1;2;3;4;5;6;7;8;9]',
'CC M',
'CC M KAPPA[1;2;3;4;5;6;7;8;9]',
'CC M D Q O H XYZ U2',
'CC KAPPA[1;2;3;4;5;6;7;8;9]',
'CC M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC HXYZ',
'CC M D Q O H XYZ U2',
'CC D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC KAPPAP[1;2;3;4;5;6;7;8;9]',
'CC M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC SIGOBS[2] M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC SIGOBS[1] M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC SIGOBS[0] M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'CC SIGOBS[0] NOSYMM[] M',
'CC SIGOBS[0] NOSYMM[] M D Q O H XYZ U2',
'CC SIGOBS[0] NOSYMM[] D Q O H XYZ U2',
'CC SIGOBS[0] NOSYMM[] M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'SIGOBS[0] NOSYMM[] M',
'SIGOBS[0] NOSYMM[] D Q O H',
'SIGOBS[0] NOSYMM[] M D Q O H XYZ U2',
'SIGOBS[0] NOSYMM[] D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'SIGOBS[0] NOSYMM[] M D Q O H XYZ U2 KAPPA[1;2;3;4;5;6;7;8;9]',
'SIGOBS[0] NOSYMM[] U3[]',
'SIGOBS[0] NOSYMM[] XYZ',
'SIGOBS[0] NOSYMM[] M D Q O H XYZ U2 U3[]',
'SIGOBS[0] NOSYMM[] D Q O H XYZ U2 U3[] KAPPA[1;2;3;4;5;6;7;8;9]',
'SIGOBS[0] NOSYMM[] M D Q O H XYZ U2 U3[] KAPPA[1;2;3;4;5;6;7;8;9]',
'SIGOBS[0] NOSYMM[] U4[]',
'SIGOBS[0] NOSYMM[] U2',
'SIGOBS[0] NOSYMM[] M D Q O H XYZ U2 U3[] U4[]',
'SIGOBS[0] NOSYMM[] D Q O H XYZ U2 U3[] U4[] KAPPA[1;2;3;4;5;6;7;8;9]',
'SIGOBS[0] NOSYMM[] M D Q O H XYZ U2 U3[] U4[] KAPPA[1;2;3;4;5;6;7;8;9]'
'\n# valid expressions:',
'# CC          - use chemical constraints',
'# SIGOBS[0]   - use sigma cutoff at [INTEGER]',
'# NOSYMM[]    - release the multipol symmetry constraints for ALL (non H) atoms (NOSYMM[C(1);S(1)] will release the symmetry for ALL but C(1) and S(1)!)',
'# M D Q O H   - refine Mono-, Di-, Quadro-, Okta- and/or Hexadecapoles',
'# XYZ         - refine coordinates of non H atoms',
'# HXYZ        - refine coordinates for H atoms (adds sintl cut at 0.5 and RESET BOND (next cycle!))',
'# U2          - refine vibrational parameter',
'# U3[]        - refine 3rd order gram-charlier parameters for given atoms (syntax: U3[C(1);C(2)] refines GC3 for C(1) and C(2))',
'# U4[]        - refine 4th order gram-charlier parameters for given atoms (syntax: U4[C(1);C(2)] refines GC4 for C(1) and C(2))',
'# KAPPA[1;2]  - refine Kappa (syntax: KAPPA[1;2;3;7] will refine Kappas 1, 2, 3 and 7)',
'# KAPPAP[1;2] - refine Kappa Prime (syntax: KAPPAP[1;2;3;7] will refine Kappa Prime 1, 2, 3 and 7)',
]

import sys,os,re,fnmatch
from copy import deepcopy
from collections import OrderedDict

print ' __  __ ____    ____  _____  ____      _   _____  _____   ____ __   __'
print ' \ \/ /|  _ \  / ___||_   _||  _ \    / \ |_   _|| ____| / ___|\ \_/ /'
print '  \  / | | | | \___ \  | |  | |_) |  / _ \  | |  |  _|  | |  _  \   / '
print '  /  \ | |_| |  ___) | | |  |  _ <  / ___ \ | |  | |___ | |_| |  | |  '
print ' /_/\_\|____/  |____/  |_|  |_| \_\/_/   \_\|_|  |_____| \____|  |_|  '

def Keytable(line,Instructions, atom_class_list):
	CCInstructions = deepcopy(Instructions)
	for i in ['M', 'D', 'Q', 'O', 'H']:
		if i in CCInstructions:
			CCInstructions.remove(i)
		else: pass	
	if atom_class_list[key_table_regex.match(line).group('ATOM')].chemcon == '' or atom_class_list[key_table_regex.match(line).group('ATOM')].chemcon == None:
		out = atom_class_list[key_table_regex.match(line).group('ATOM')].print_key(Instructions)
	elif atom_class_list[key_table_regex.match(line).group('ATOM')].chemcon.strip() in atom_class_list.keys() and 'CC' not in Instructions:
		out = atom_class_list[key_table_regex.match(line).group('ATOM')].print_key(Instructions, atom_class_list[atom_class_list[key_table_regex.match(line).group('ATOM')].chemcon.strip()].dict)
	elif atom_class_list[key_table_regex.match(line).group('ATOM')].chemcon.strip() in atom_class_list.keys() and 'CC' in Instruct:
		out = atom_class_list[key_table_regex.match(line).group('ATOM')].print_key(CCInstructions, atom_class_list[atom_class_list[key_table_regex.match(line).group('ATOM')].chemcon.strip()].dict)
	else:
		pass
	return out + '\n'

def contains(list, filter):
	for x in list:
		if filter(x):
			return True
	return False

class Atom():
	def __init__(self,chemcon, atomtype, kappa):
		self.dict = {}
		self.chemcon = chemcon
		self.atomtype = atomtype
		self.special = False
		self.U3 = False
		self.U4 = False
		self.kappa = kappa
		self.order = ['XYZ', 'U2', 'U3', 'U4', 'M', 'D', 'Q', 'O', 'H']
	def print_key(self, Instructions, keydict=None):
		if keydict==None:
			keydict = self.dict
		line = ''
		line = self.dict['ATOM'].ljust(7)
		for i in self.order:
			if self.atomtype == 'H' and i == 'XYZ' and 'HXYZ' in Instructions:
				line = '{} {}'.format(line,keydict[i].strip().replace('0','1')) 
			elif 'NOSYMM' in Instructions and i in Instructions and not i == 'U3' and not i == 'U4':
				if self.special:
					if i in Instructions:
						if self.atomtype == 'H' and i == 'XYZ':
							line = '{} {}'.format(line,re.sub('1', '0', keydict[i]).strip())
						else:
							line = '{} {}'.format(line,keydict[i].strip())
				else:
					if i == 'M':
						line = '{} {}'.format(line,keydict[i].strip())
					elif self.atomtype == 'H' and i == 'XYZ':
						line = '{} {}'.format(line,re.sub('1', '0', keydict[i]).strip())
					else:
						line = '{} {}'.format(line,re.sub('0', '1', keydict[i]).strip())
			elif not 'NOSYMM' in Instructions and i in Instructions and not i == 'U3' and not i == 'U4':
				if i == 'M':
					line = '{} {}'.format(line,'10')
				elif self.atomtype == 'H' and i == 'XYZ':
					line = '{} {}'.format(line,re.sub('1', '0', keydict[i]).strip())
				else:
					line = '{} {}'.format(line,keydict[i].strip())
			elif i == 'U3' and i in Instructions and self.U3:
				line = '{} {}'.format(line,re.sub('0', '1', keydict[i]).strip())
			elif i == 'U4' and i in Instructions and self.U4:
				line = '{} {}'.format(line,re.sub('0', '1', keydict[i]).strip())
			elif not i in Instructions or i == 'U3' or i == 'U4':
				line = '{} {}'.format(line,re.sub('1', '0', keydict[i]).strip())
		return line

atom_table_regex = re.compile('(?P<ATOM>[a-zA-Z\(\)0-9]+)\s+(?P<ATOM0>[a-zA-Z\(\)0-9]+)\s+(?P<AX1>[XYZxyz])\s+(?P<ATOM1>[a-zA-Z\(\)0-9]+)\s+(?P<ATOM2>[a-zA-Z\(\)0-9]+)\s+(?P<AX2>[XYZxyz])\s+(?P<RL>[RLrl])\s+(?P<TP>\d)\s+(?P<TBL>\d+)\s+(?P<KAP>\d+)\s+(?P<LMX>\d)\s+(?P<SITESYM>[0-9a-zA-Z_]*){0,1}(\s+(?P<CHEMCON>.*)\s*)*\n')
key_table_regex = re.compile('(?P<ATOM>[a-zA-Z\(\)0-9]+)\s+(?P<XYZ>[01]{3})\s+(?P<U2>[01]{6})\s+(?P<U3>[01]{10})\s+(?P<U4>[01]{15})\s+(?P<M>[01]{2})\s+(?P<D>[01]{3})\s+(?P<Q>[01]{5})\s+(?P<O>[01]{7})\s+(?P<H>[01]{9}).*\n')
scat_table_regex = re.compile('(?P<SCAT>[a-zA-Z\+\-]+)\s+(?P<CORE>\w+)\s+(?P<SPHV>\w+)\s+(?P<DEFV>\w+)\s+(?P<S1>[0-9\-]+)\s+(?P<S2>[0-9\-]+)\s+(?P<S3>[0-9\-]+)\s+(?P<S4>[0-9\-]+)\s+(?P<P2>[0-9\-]+)\s+(?P<P3>[0-9\-]+)\s+(?P<P4>[0-9\-]+)\s+(?P<D3>[0-9\-]+)\s+(?P<D4>[0-9\-]+)\s+(?P<F4>[0-9\-]+)\s+(?P<S5>[0-9\-]+)\s+(?P<P5>[0-9\-]+)\s+(?P<S6>[0-9\-]+)\s+(?P<P6>[0-9\-]+)\s+(?P<D5>[0-9\-]+)\s+(?P<S7>[0-9\-]+)\s+(?P<D6>[0-9\-]+)\s+(?P<F5>[0-9\-]+)\s+(?P<DELF>[0-9\.\-]+)\s+(?P<DELFf>[0-9\.\-]+)\s+(?P<NSCTL>[0-9\.\-]+)\s*')
vib_con_regexp = re.compile('CON(\s+[\d.-]+\s+[U\d/]+){4,5}\s+=\s+0\s+.*')#{4,5} statt {5}!
reset_regexp = re.compile('RESET\s+BOND(\s+[a-zA-Z\(\)0-9]+){2}\s+[\d.]+\s+')

atom_class_list = OrderedDict()
scat_list = OrderedDict()
vibcons = []
resets = []
atom_list = []
key_list = []
KAPPAlist = []

####################################################
#                  file handler                    #
####################################################
def ReadFile(File_Name):
	try:
		Out = []
		with open(File_Name, 'r') as File_Name:
			File = File_Name.readlines()
			for line in File:
				if not re.match('^#',line) or not re.match('^\n|^\s',line):
					Out.append(line)
				else: pass
			return Out
	except IOError:
		if File_Name == 'xd.inst':
			xdinst = open('xd.inst','w')
			for i in standard:
				xdinst.write('{}\n'.format(i))
			print ' No instruction file found!\n xd.inst written! Please edit and run XDStratgen.py again!'
		elif File_Name == 'xd.const':
			try:
				print ' No constraints file found!\n > running XDConstraint.py!'
				os.system('x:\programme\skripte\XDConstraint.py')
			except:
				print ' No constraints file found!\n Please run XDConstraint.py!'
		else:
			print ' error: {} not found!'.format(File_Name)
		sys.exit()

Master_Name = raw_input(' enter filename [xd.mas]: ') or 'xd.mas'
if not re.search('\.mas', Master_Name):
	Master_Name = Master_Name + '.mas'
MasterFile = ReadFile(Master_Name)

Instructions_raw = ReadFile('xd.inst')

# insert constraints from xd.const?
insert_const = raw_input(' insert xd.const [Y/N] N? ').upper() or 'N'
if insert_const == 'Y':
	Constraints  = ReadFile('xd.const')
	insert_const = True
else:
	Constraints  = []
	insert_const = False


####################################################
#                                                  #
####################################################
kappas = []
for line in MasterFile:
	if atom_table_regex.search(line):
		atom = atom_table_regex.match(line).groupdict()
		atom_class_list[atom['ATOM']] = Atom(atom['CHEMCON'],scat_list.keys()[int(atom['TBL'])-1],atom['KAP'])
		# was kappas is never used?!?
		# if not re.match('^H\(.+\)',atom['ATOM']):# only append non H kappas
		kappas.append(int(atom['KAP']))
	elif key_table_regex.search(line):
		key = key_table_regex.match(line).groupdict()
		atom_class_list[key['ATOM']].dict = key
	elif scat_table_regex.search(line):
		scat = scat_table_regex.match(line).groupdict()
		scat_list[scat['SCAT']] = scat
	elif line == 'END ATOM':
		break
kappas = sorted(list(set(kappas)))

for con in Constraints:
	if vib_con_regexp.search(con):
		vibcons.append(con)
	elif reset_regexp.search(con): 
		resets.append(con)

# UGLY BUT NECESSARY
Instructions = []
for j in Instructions_raw:
	if re.match('^#|^\n|^\s', j):
		pass
	else:
		Instructions.append(j)

re_kappa  = re.compile('KAPPA\[.*\]')
re_kappap = re.compile('KAPPAP\[.*\]')
SELECT = re.compile('SELECT\s+cycle\s+[0-9\-]+\s+dampk\s+[0-9\.]+\s+cmin\s+[0-9\.]+\s+cmax\s+[0-9\.]+\s+eigcut\s+[d0-9\.\-]+\s+convcrit\s+[d0-9\.\-]+\s+')
for j in range(len(Instructions)):
	outfile = open('xd%02d'%(j+1)+'.mas','w')
	TempMaster = deepcopy(MasterFile)
	KAPPAIns = {}
	NOSYMM = ['H']
	U3 = []
	U4 = []
	SIGOBS = '3'
	Instruct = Instructions[j].split()
	keytable_instructions = []
	constraints_instructions = []

	for i in Instruct:
		if i == 'SCALE':
			for line in TempMaster:
				if re.match(SELECT, line):
					cycles = str(abs(int(line.split()[2])))
					TempMaster[TempMaster.index(line)] = line = (re.sub('cycle\s+-*\d+','cycle -'+cycles,line))
				if re.match('^!RESET\s+', line) and 'HXYZ' not in Instruct:
					for reset in reversed(resets):
						TempMaster.insert(TempMaster.index(line)+1,reset)
				####################################### das ist hier falsch! ##########################################
				#if re.match('!CON\s+num1\s+par1/iat1\s+num2\s+par2/iat2\s+...\s+=\s+num0',line) and 'U2' in Instruct:
				#if re.match('!CON\s+',line) and 'U2' in Instruct:
				#	for con in reversed(vibcons):
				#		TempMaster.insert(TempMaster.index(line)+1,con)
			break
		
		if i.startswith('SIGOBS'):
			SIGOBS = re.sub('[\[\];]',' ', i).split()[1]
		elif i.startswith('NOSYMM'):
			for j in re.sub('[\[\];]',' ', i).split()[1:]:
				NOSYMM.append(j)
			i = 'NOSYMM'
		elif i.startswith('U3'):
			U3 = re.sub('[\[\];]',' ', i).split()[1:]
			i = 'U3'
		elif i.startswith('U4'):
			U4 = re.sub('[\[\];]',' ', i).split()[1:]
			i = 'U4'
		elif i.startswith('KAPPA'):
			KAPPA = re.sub('[\[\];]',' ', i).split()
			KAPPAIns[KAPPA[0]] = KAPPA[1:]
		if i in ['CC', 'HXYZ', 'XYZ', 'U2', 'U3', 'U4','M', 'D', 'Q', 'O', 'H', 'NOSYMM']:
			keytable_instructions.append(i)
	
	for atom in atom_class_list:
		if atom_class_list[atom].atomtype in NOSYMM:
			atom_class_list[atom].special = True
		if atom_class_list[atom].dict['ATOM'] in U3:
			atom_class_list[atom].U3 = True
		if atom_class_list[atom].dict['ATOM'] in U4:
			atom_class_list[atom].U4 = True 		

	set_reset = False
	set_const = False
	try:
		if cycles: pass
	except:
		print ' SCALE step is missing in instruction file!'
		sys.exit()
	firstkappa = 0# -> finde nur das erste kappa
	for line in TempMaster:
		if atom_table_regex.search(line):
			#U3, U4 TP in Atom List
			if atom_class_list[atom_table_regex.match(line).groupdict()['ATOM']].U3:
				TempMaster[TempMaster.index(line)] = line = re.sub('\s+[124]', '   3',line, count = 1)
			if atom_class_list[atom_table_regex.match(line).groupdict()['ATOM']].U4:
				TempMaster[TempMaster.index(line)] = line = re.sub('\s+[123]', '   4',line, count = 1)
		if 'CC' not in keytable_instructions:
			if atom_table_regex.search(line):
				atom = atom_table_regex.match(line).groupdict()
				if atom['CHEMCON']:	
					TempMaster[TempMaster.index(line)] = line = ''.join(line.rsplit(atom['CHEMCON'],1))
		if re.match(SELECT, line):
			TempMaster[TempMaster.index(line)] = line = re.sub('cycle\s+\d+','cycle '+cycles,line)
		#if re.match('!RESET    bond C\(1\) H\(1\) 1.09 ...', line) and 'XYZ' in Instruct and 'HXYZ' not in Instruct:
		if re.match('^!RESET\s+', line) and 'XYZ' in Instruct and 'HXYZ' not in Instruct and not set_reset and insert_const:
				set_reset = True
				for reset in reversed(resets):
					TempMaster.insert(TempMaster.index(line)+1,reset)
		#if re.match('!CON      num1 par1/iat1 num2 par2/iat2 ... = num0',line) and 'U2' in Instruct:
		if re.search('^!CON\s+',line) and 'U2' in Instruct and not set_const and insert_const:
			set_const = True
			for con in reversed(vibcons):
				TempMaster.insert(TempMaster.index(line)+1,con)

		####################################################
		#        very sloppy KAPPA const workaround!       #
		####################################################
		if re.search('^!CON\s+1\s+KS\s*/\s*\d+',line) and [True for x in Instruct if re_kappa.match(x)]:#sloppy! Another loop in the loop!
			TempMaster[TempMaster.index(line)] = line = re.sub('^!CON','CON',line)
		if re.search('^!CON\s+1\s+K0\s*/\s*\d+',line) and [True for x in Instruct if re_kappap.match(x)]:#sloppy! Another loop in the loop!
			TempMaster[TempMaster.index(line)] = line = re.sub('^!CON','CON',line)

		if re.match('SKIP\s+[\*]{0,1}obs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sigobs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sinthl\s+[d\d.]+\s+[d\d.]+\s+', line) and SIGOBS:
			#TempMaster[TempMaster.index(line)] = line = re.sub('sigobs\s+[d\d.]+', '*sigobs '+SIGOBS, re.sub('\*', '',line))
			TempMaster[TempMaster.index(line)] = line = re.sub('sigobs\s+[d\d.]+', '*sigobs '+SIGOBS, re.sub('\*sigobs','sigobs',line))
		if re.match('SKIP\s+[\*]{0,1}obs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sigobs\s+[d\d.]+\s+[d\d.]+\s+[\*]{0,1}sinthl\s+[d\d.]+\s+[d\d.]+\s+', line) and 'HXYZ' in keytable_instructions:
			TempMaster[TempMaster.index(line)] = line = re.sub('sinthl\s+[d\d.]+\s+[d\d.]+', '*sinthl 0.0 0.5', re.sub('\*sinthl', 'sinthl',line))
		if re.match('SELECT\s+\**model\s+4\s+[23]\s+1\s+0\s+based_on\s+F\^2\s+test\s+verbose\s+1\s+', line) and contains(atom_class_list, lambda x: atom_class_list[x].U3 == True):
			TempMaster[TempMaster.index(line)] = line = re.sub('model\s+4\s+[23]\s+1\s+0', 'model 4 3 1 0',line)
		if re.match('FOUR\s+fmod1\s+4\s+[23]\s+0\s+0\s+fmod2\s+-1\s+[23]\s+0\s+0', line) and contains(atom_class_list, lambda x: atom_class_list[x].U3 == True):
			TempMaster[TempMaster.index(line)] = line = re.sub('fmod2\s+-1\s+[23]\s+0\s+0', 'fmod2 -1 3 0 0',re.sub('fmod1\s+4\s+[23]\s+0\s+0\s+', 'fmod1 4 3 0 0 ',line))
		if re.match('SELECT\s+\**model\s+4\s+[23]\s+1\s+0\s+based_on\s+F\^2\s+test\s+verbose\s+1\s+', line) and contains(atom_class_list, lambda x: atom_class_list[x].U4 == True):
			TempMaster[TempMaster.index(line)] = line = re.sub('model\s+4\s+[23]\s+1\s+0', 'model 4 4 1 0',line)
		if re.match('FOUR\s+fmod1\s+4\s+[23]\s+0\s+0\s+fmod2\s+-1\s+[23]\s+0\s+0\s+', line) and contains(atom_class_list, lambda x: atom_class_list[x].U4 == True):
			TempMaster[TempMaster.index(line)] = line = re.sub('fmod2\s+-1\s+[23]\s+0\s+0', 'fmod2 -1 4 0 0',re.sub('fmod1\s+4\s+[23]\s+0\s+0', 'fmod1 4 4 0 0 ',line))
		if key_table_regex.search(line):
			TempMaster[TempMaster.index(line)] = line = Keytable(line, keytable_instructions, atom_class_list)
		#outfile.write(line)
		if re.match('^KAPPA\s+[01]+\s*',line) and firstkappa == 0:
			firstkappa = TempMaster.index(line)

	####################################################
	#                    KAPPA                         #
	####################################################
	if not 'KAPPA' in KAPPAIns.keys() and not 'KAPPAP' in KAPPAIns.keys():# mach kappas aus wenn nicht in inst angeschalted
		for i in kappas:
			TempMaster[firstkappa+int(i)-1] = re.sub('[0,1]{2}','00',TempMaster[firstkappa+int(i)-1],1)
	if 'KAPPA' in KAPPAIns.keys() and not 'KAPPAP' in KAPPAIns.keys():
		for i in KAPPAIns['KAPPA']:
			TempMaster[firstkappa+int(i)-1] = re.sub('[0,1]{2}','10',TempMaster[firstkappa+int(i)-1],1)
	if 'KAPPAP' in KAPPAIns.keys() and not 'KAPPA' in KAPPAIns.keys():
		for i in KAPPAIns['KAPPAP']:
			TempMaster[firstkappa+int(i)-1] = re.sub('[0,1]{2}','01',TempMaster[firstkappa+int(i)-1],1)
	if 'KAPPA' in KAPPAIns.keys() and 'KAPPAP' in KAPPAIns.keys():
		for i in KAPPAIns['KAPPA']:
			if i in KAPPAIns['KAPPAP']:
				TempMaster[firstkappa+int(i)-1] = re.sub('[0,1]{2}','11',TempMaster[firstkappa+int(i)-1],1)
			else:
				TempMaster[firstkappa+int(i)-1] = re.sub('[0,1]{2}','10',TempMaster[firstkappa+int(i)-1],1)
		for i in KAPPAIns['KAPPAP']:
			if i not in KAPPAIns['KAPPA']:
				TempMaster[firstkappa+int(i)-1] = re.sub('[0,1]{2}','01',TempMaster[firstkappa+int(i)-1],1)
			else: pass

	for line in TempMaster:
		outfile.write(line)
	print ' -> {} written!\r'.format(outfile.name),
	outfile.close()

#print ' __  __ ____    ____  _____  ____      _   _____   ____  _____  _   _ '
#print ' \ \/ /|  _ \  / ___||_   _||  _ \    / \ |_   _| / ___|| ____|| \ | |'
#print '  \  / | | | | \___ \  | |  | |_) |  / _ \  | |  | |  _ |  _|  |  \| |'
#print '  /  \ | |_| |  ___) | | |  |  _ <  / ___ \ | |  | |_| || |___ | |\  |'
#print ' /_/\_\|____/  |____/  |_|  |_| \_\/_/   \_\|_|   \____||_____||_| \_|'