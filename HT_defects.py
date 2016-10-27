############################################
#
# High-throughput work flow example 
# System: In2O3
# Point-defects: Intrinsic vacancies, anti-sites and interstitials
# Functional: DFT GGA_PBE
#
############################################
# import from python and pylada
from pylada.vasp.relax import Relax
from pylada.crystal import read, supercell, neighbors
import numpy as np
import pylada_defects
from random import choice, random
from copy import deepcopy
######################################
# find nghs 
def find_ngh_indices(atomind,structure):
    """ function to return ngh indices from given atomind
   
    """
    
    nghs=neighbors(structure,1,structure[atomind],0.2)
    ngh_indices=[]

    for ii in range(len(structure)):
        if any([ all(structure[ii].pos==ngh[0].pos) for ngh in nghs]) and any([ structure[ii].type==ngh[0].type for ngh in nghs]):
           ngh_indices.append(ii)

    return ngh_indices

# condensed magmom
def gather_magmom(mgm):
    mgm=mgm.split()
    hlp=[ [ mgm[0] ] ]

    for ii in range(1,len(mgm)):
        if mgm[ii]==mgm[ii-1]:
            hlp[len(hlp)-1].append(mgm[ii])
        else:
            hlp.append([mgm[ii]])

    out_str=''

    for hh in hlp:
        out_str=out_str+'%s*%1.1f ' %(len(hh),float(hh[0]))

    return out_str

# generating k-point 
def gen_kpts(structure,density):
    import numpy as np
    import numpy.linalg as lalg

    nkpts=density/len(structure)
    direct_cell=np.transpose(structure.cell)
    rec_cell=2*np.pi*np.transpose(lalg.inv(direct_cell))

    b1=np.sqrt(np.dot(rec_cell[0],rec_cell[0]))
    b2=np.sqrt(np.dot(rec_cell[1],rec_cell[1]))
    b3=np.sqrt(np.dot(rec_cell[2],rec_cell[2]))

    step=(b1*b2*b3/nkpts)**(1./3)

    n1=int(round(b1/step))
    if np.mod(n1,2)!=0: n1=n1-1
    n2=int(round(b2/step))
    if np.mod(n2,2)!=0: n2=n2-1
    n3=int(round(b3/step))
    if np.mod(n3,2)!=0: n3=n3-1
    print n1,n2,n3

    if n1==0:n1=1
    if n2==0:n2=1
    if n3==0:n3=1

    return "\n0\nGamma\n%2i %2i %2i\n 0. 0. 0.\n" %(n1,n2,n3)

############### setting up the functional
vasp=Relax()

vasp.program = '/home/vstevano/bin/vasp'

pseudoDir = '/home/vstevano/software/pseudos'
vasp.add_specie = "In", pseudoDir + "/In"
vasp.add_specie = "O", pseudoDir + "/O"

vasp.prec       = "accurate"
vasp.encut      = 340.
vasp.ismear     = 0
vasp.sigma      = 0.05
vasp.ediff      = 1.0e-6
vasp.ediffg     = -0.01
vasp.convergence= 1.0e-6
vasp.minrelsteps= 4
vasp.nsw        = 75
vasp.lwave      = True
vasp.lorbit     = 10
vasp.lplane     = True
vasp.addgrid    = True
vasp.npar       = 8
vasp.isym       = 0
vasp.lcharg     = True
vasp.lwave      = True
vasp.lmaxmix    = 4
vasp.loptics    = False
vasp.lpead      = False
vasp.algo = "Normal"
vasp.relaxation = "ionic"
vasp.maxiter = 10
vasp.keep_steps = True
vasp.first_trial = { "kpoints": "\n0\nAuto\n10", "encut": 0.9 }

############### setting up the structures

# input bulk primitive cell
In2O3prim = read.poscar('POSCAR_In2O3')

# create bulk supercell
In2O3_sc = supercell(In2O3prim,np.diag([In2O3prim.cell[0][0]*2., In2O3prim.cell[1][1]*2., In2O3prim.cell[2][2]*2.]))

# create list of job folders
calcs = ['epsilon', 'SC', 'Ini', 'VIn', 'OIn']

structures = {}

for calc in calcs:
    if calc=='epsilon':
        structure=deepcopy(In2O3prim)
        structures[calc]=structure

    else:
        structure=deepcopy(In2O3_sc)

        # supercell
        if calc=='SC':
            structures[calc]=structure

        # interstitial
        if calc=='Ini':
            ints_list = []
            ints_list = pylada_defects.get_interstitials(structure)
            for j in range(len(ints_list)):
                structure2 = deepcopy(structure)
                structure2.add_atom(ints_list[j][0], ints_list[j][1], ints_list[j][2], 'In')
                structure2[-1].spin = 1.
                key = 'Ini'+'_'+str(j)
                structures[key]=structure2

        # vacancy
        if calc=='VIn':
            vacancy_indices = pylada_defects.get_atom_indices('In', structure)
            for k in range(len(vacancy_indices)):
                structure2 = deepcopy(structure)
                for ii in find_ngh_indices(vacancy_indices[k], structure2):
                    structure2[ii].pos = structure2[ii].pos + 0.1*np.array([2*random()-1, 2*random()-1, 2*random()-1])
                    structure2[ii].spin = 1.
                del structure2[vacancy_indices[k]]
                key = 'VIn'+'_'+str(k)
                structures[key]=structure2

        # subsitution
        if calc=='OIn':
            vacancy_indices = pylada_defects.get_atom_indices('In', structure)
            for k in range(len(vacancy_indices)):
                structure2 = deepcopy(structure)
                for ii in find_ngh_indices(vacancy_indices[k], structure2):
                    structure2[ii].pos = structure2[ii].pos + 0.1*np.array([2*random()-1, 2*random()-1, 2*random()-1])
                    structure2[ii].spin = 1.
                structure2[vacancy_indices[k]].type = 'O'
                key = 'OIn'+'_'+str(k)
                structures[key]=structure2


############### setting up the jobfolder
from IPython.core.interactiveshell import InteractiveShell
from pylada.jobfolder import JobFolder
from pylada import interactive
from copy import deepcopy

# Job dictionary.
jobfolder = JobFolder()

# loop over material-lattice pairs.
for name in structures:
    if name=='epsilon':
        structure=structures[name]
  
        # job folder for this lattice.
        job = jobfolder / name
        vasp_individual = deepcopy(vasp)
  
        vasp_individual.relaxation = "ionic"
        vasp_individual.add_keyword('lepsilon',True)
        vasp_individual.add_keyword('lrpa',False)
        vasp_individual.ibrion = 7

        vasp_individual.kpoints=gen_kpts(structure,4000)

        job.functional = vasp_individual
        job.params["structure"] = structure

    else:
        structure=structures[name]
        
        # job folder for this lattice.
        job = jobfolder / 'defects' / name 
        vasp_individual = deepcopy(vasp)
        vasp_individual.ispin=2
  
        magmom=''
        for ii in range(len(structure)):
            if hasattr(structure[ii],'spin'):
                magmom=magmom+'%s  ' %(structure[ii].spin)
            else:
                magmom=magmom+'%s  ' %(0.)

        vasp_individual.magmom=gather_magmom(magmom)
        vasp_individual.kpoints="\n0\nGamma\n2  2  2\n0. 0. 0.\n"
        vasp_individual.ediffg = -0.01

        job.functional = vasp_individual
        job.params["structure"] = structure


%load_ext pylada
%savefolders defects.pickle jobfolder
#%launch scattered --account=NRELMatDB --walltime=36:00:00 --ppn=24  --queue=batch
