#################################################
# 
# anuj.goyal@nrel.gov
# Date: October 26 2016
#
##################################################

# imports from pylada
from pylada.crystal import Structure, primitive
from pylada.crystal import read, write, neighbors, supercell
from pylada.crystal.defects import reindex_sites
from pylada.vasp import Extract, MassExtract

# imports from python
from pyspglib import spglib
from glob import iglob
from itertools import *
from scipy.spatial import distance
from collections import Counter
from copy import deepcopy
from quantities import eV
import numpy as np
import sys
import math
import tess
import os

##########################################################
def compute_voronoi(points):                                                        
    """ Function to return the python container having information of Voronoi cells for            
        for given points in 3D space                                                               

    Parameters                                                                                          
        pts = numpy array of points coordinates                                                    

    Returns                                                                                         
        container = python contrainer having information of Voronoi cells                          
    """

    P = np.array(points)
    
    # box limits along x, y, and z axis
    Lx = 50
    Ly = 50
    Lz = 50

    cntr = tess.Container(P, ((-50, -50, -50),(50, 50, 50)), periodic=(False, False, False))

    return cntr

##########################################################
def calculate_midpoint(p1, p2):
    """ Calculate the midpoint given coordinates
    
    Parameters                                                                                          
        p1, p2 = numpy array of point coordinates                                                  
                                                                                                   
    Returns                                                                                         
        midpoint = numpy array of midpoint coordinates                                             
    """

    return((p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0, (p1[2]+p2[2])/2.0)

##########################################################
def calculate_polygon_centroid(poly_pts):
    """ Function to calculate the centroid of non-self-intersecting polygon

    Parameters
        pts = numpy array of coordinates of vertices of polygon

    Returns
        centroid = numpy array of centroid coordinates
    """

    P = np.array(poly_pts)
    C = np.mean(P, axis=0)

    return C

##########################################################
def neighbor_list(list):
    """ Function to form unique neighboring pairs along the polygon perimeter

    Parameters 
        list = list of indicies of Voronoi vertices forming the perimeter                          

    Returns
        list = list of neighboring pairs as tuples
    """

    i = 0
    while i + 1 < len(list):
        yield (list[i], list[i+1])
        i += 1
    else:
        yield (list[i], list[0])

##########################################################
def get_vertices(site_num, cntr):
    """ Function that returns vertices of the Voronoi associated with given site

    Parameters
        site_num = number for the lattice site of interest
        cntr = python contrainer having information of Voronoi cells

    Returns
        vertices = numpy array of Voronoi vertices coordinates
    """

    list_voronoi_vertices = cntr[site_num].vertices()
    V = list_voronoi_vertices

    # convert the list to numpy array
    V = np.asarray(V)

    return V

##########################################################
def get_edgecenter(site_num, cntr):
    """ Function that returns vertices unique edge centers of the Voronoi associated with specific\
 lattice site

    Parameters
        site_num = number for the lattice site of interest
        cntr = python contrainer having information of Voronoi cells

    Returns
        Edge center = numpy array of Voronoi edge center coordinates
    """

    list_face_vertices_indices = cntr[site_num].face_vertices()

    V_vertices = get_vertices(site_num, cntr)

    all_midpoint = []

    for face in list_face_vertices_indices:
        for(x,y) in neighbor_list(face):
            midpoint = calculate_midpoint(V_vertices[x], V_vertices[y])
            all_midpoint.append(midpoint)

    #using set so to choose only unique edge centers
    S = set(all_midpoint)

    #converting set to list
    Ec = list(S)

    #converting list to numpy array
    Ec = np.asarray(Ec)

    return Ec

##########################################################
def get_facecentroid(site_num, cntr):
    """Function the returns vertices of face centers of the Voronoi associated with specific latti\
ce site

    Parameters
        site_num = number for the lattice site of interest
        cntr = python contrainer having information of Voronoi cells

    Returns
        Face center = numpy array of Voronoi face center coordinates                               
    """

    list_face_vertices_indices = cntr[site_num].face_vertices()

    V_vertices = get_vertices(site_num, cntr)

    list_face_centroid = []

    for face in list_face_vertices_indices:
        l = []
        for j in face:
            vv  = V_vertices[j]
            l.append(vv)
        l = np.asarray(l)
        pc = calculate_polygon_centroid(l)
        list_face_centroid.append(pc.tolist())

    Fc = list_face_centroid

    # converting list to numpy array                                                               
    Fc = np.asarray(Fc)

    return Fc

##########################################################
def get_all_interstitials(prim_structure, positions):
    """ function to return list of all interstitial sites using Voronoi.py

    Parameters
        prim = pylada primitive structure
        positions = positions (numpy array) to compute interstitials for

    Returns
        Inst_list = list of list, with inner list containing ['atom type', [x,y,z]]
    """

    ints_list = []

    for site_num in range(len(positions)):
        site = np.array([positions[site_num][0], positions[site_num][1], positions[site_num][2]])

        site_ngh = neighbors(prim_structure, 13, site, 0.1)

        points = [site]

        ### creating list with site and its neighbors
        for i in range(len(site_ngh)):
                a = site + site_ngh[i][1]
                points.append(a)

        ### converting list to numpy array
        points = np.asarray(points)

        ### using tess object cntr to compute voronoi
        cntr = compute_voronoi(points)

        ### Voronoi vertices
        ### the first position in points is the site, therefore '0'
        v = get_vertices(0, cntr)

        for i in range(len(v)):
            ints_list.append(['B', v[i].tolist()])

        ### Voronoi face centers
        f = get_facecentroid(0, cntr)

        for j in range(len(f)):
            ints_list.append(['C', f[j].tolist()])

        ### Voronoi edge centers
        e = get_edgecenter(0, cntr)

        for k in range(len(e)):
            ints_list.append(['N', e[k].tolist()])

    ### return list of list ['Atom type', [x,y,z]]
    return ints_list

##########################################################
def get_pos_in_prim_cell(prim, a):
    """ Function to to map positions onto the primitive cell

    Parameters
        prim = pylada primitive cell
        a = cartesian coordinates of position

    Returns
        a2 = cartesian coordinates, such that fractional coordination = [0,1)
    """
    
    p = deepcopy(prim)
    a1 = np.array(a)
    inv_cell = np.linalg.inv(p.cell)

    frac_a1 = np.dot(inv_cell, a1)
    for m in range(len(frac_a1)):
        if frac_a1[m] < 0.: frac_a1[m] = frac_a1[m] + 1.
        if frac_a1[m] >= 0.999: frac_a1[m] = frac_a1[m] - 1.
    a2 = np.dot(p.cell, frac_a1)

    return a2

##########################################################
def get_ints_in_prim_cell(prim, positions):
    """ Function to to map positions onto the primitive cell
    
    Parameters
        prim = pylada primitive cell
        positions = list of sites (as list ['element',[x,y,z]])

    Returns
        int_pos_list1 = unique list of sites (as numpy array) within primitive cell
    """
    
    prim1 = deepcopy(prim)
    ints1 = deepcopy(positions)
    inverse_cell1 = np.linalg.inv(prim1.cell)
    int_pos_list1 = []
    
    for i in range(len(ints1)):
        a1 = np.array(ints1[i][1])
        frac_a1 = np.dot(inverse_cell1, a1)
        for m1 in range(len(frac_a1)):
            if frac_a1[m1] < 0.: frac_a1[m1] = frac_a1[m1]+1.
            if frac_a1[m1] >= 0.999: frac_a1[m1] = frac_a1[m1]-1.
        a2 = np.dot(prim1.cell, frac_a1)
        int_pos_list1.append(a2)
    
    return int_pos_list1

##########################################################
def get_unique_wyckoff(prim):
    """ Function to find unique wyckoff sites in the primitive cell
    
    Parameters
        prim = pylada primitive cell

    Returns
        unique_list = unique list of sites (as numpy array) in primitive cell
    """
    p = deepcopy(prim)

    # compute space group for the given primitive cell using spglib
    sym = spglib.get_symmetry(p, 0.1)

    # compute inverce cell of the primitive cell
    inverse_cell = np.linalg.inv(p.cell)
    
    dummy_list = []
    wyckoff_list = []

    for i in range(len(p)):
        a = p[i].pos
        a3 = get_pos_in_prim_cell(p, a)
        frac_a = np.dot(inverse_cell, a3)
        symm_list = []
        for j in range(len(sym['rotations'])):
            # rotation matrix from sym
            R = sym['rotations'][j]
            # translation vector from sym
            Tt = sym['translations'][j]
            frac_symm_a = np.dot(R, frac_a)+Tt
            symm_a = np.dot(p.cell, frac_symm_a)
            symm_list.append(symm_a)
            symm_a2 = get_pos_in_prim_cell(p, symm_a)
            symm_list.append(symm_a2)
        # loop to find symmetrical equivalent positions
        for k in range(i+1, len(p)):
            b = p[k].pos
            b2 = get_pos_in_prim_cell(p, b)
            # check distance between positions in pos and symm_list
            if any(distance.euclidean(b2,c) < 0.1 for c in symm_list):
                p[k].pos = p[i].pos
        dummy_list = a.tolist()
        dummy_list.append(p[i].type)
        wyckoff_list.append(dummy_list)
    
    ### getting unique array of positions
    unique_list = [list(t) for t in set(map(tuple, wyckoff_list))]
    
    return unique_list

##########################################################
def get_unique_ints(prim, int_pos, ttol=0.5):
    """ Function to find unique interstitial sites in the primitive cell
    
    Parameters
        prim = pylada primitive cell
        int_pos = list of interstitial positions in primitive cell
        ttol = tolerance, default = 0.5

    Returns
        unique_int_list = list of unique intestital sites (cartesian (x,y,z) as list) in primitive cell
    """
    prim2 = deepcopy(prim)
    pos2 = deepcopy(int_pos)
    sym2 = spglib.get_symmetry(prim2, 0.1)

    inverse_cell2 = np.linalg.inv(prim2.cell)
    int_list2 = []

    for i in range(len(pos2)):
    #   numpy 1D array into column vector with shape listed as (3,)
        a4 = pos2[i]
        frac_a4 = np.dot(inverse_cell2, a4)
        symm_list2 = []
        for j in range(len(sym2['rotations'])):
            R = sym2['rotations'][j]
            Tt = sym2['translations'][j]
            frac_symm_a4 = np.dot(R, frac_a4) + Tt
            symm_a4 = np.dot(prim2.cell, frac_symm_a4)
            symm_list2.append(symm_a4)
            symm_a5 = get_pos_in_prim_cell(prim2, symm_a4)
            symm_list2.append(symm_a5)
        # loop to find symmetrical equivalent positions    
        for k in range(i+1, len(pos2)):
            b2 = pos2[k]
            # check distance between positions in pos and symm_list
            if any(distance.euclidean(b2,c2) < ttol for c2 in symm_list2):
                pos2[k] = pos2[i]
        int_list2.append(pos2[i])

    ### getting the unique list of interstitials
    unique_int_list = [list(t) for t in set(map(tuple, int_list2))]
    
    return unique_int_list

##########################################################
def get_interstitials(structure, ttol=0.5):
    """ Function to return unique interstitial sites in the given structure
    
    Parameters
        structure = poscar file (using read.poscar)
        ttol = tolerance

    Returns
        unique_int_list = list of unique intestital sites (cartesian (x,y,z) as list) in the given structure
    """

    s = deepcopy(structure)
    prim = primitive(s)
    spg = spglib.get_spacegroup(prim, 0.1)

    ### Step 1: get unique sites in the primitive of the given structure
    uniq_sites = get_unique_wyckoff(prim)

    ### Step 2: get all interstitial sites from Voronoi method
    ints2 = get_all_interstitials(prim, uniq_sites)

    ### get interstital sites within primitive cell
    ints_prim = get_ints_in_prim_cell(prim, ints2)

    ### Step 3: get unique interstitials after symmetry analysis
    ints3 = get_unique_ints(prim, ints_prim, ttol=ttol)
        
    return ints3

##########################################################
def get_atom_indices(atomtype, structure):
    """ function to return all indices corresponding to unique Wyckoff positions for given atomtype

    Parameters
          atomtype = 'A'
          structure = pylada structure object

    Returns
          defect_indices = list with all unique indicies for given atomtype
    """

    unq_wyck = []
    atom_indices = []
    unq_wyck = get_unique_wyckoff(structure)
    for i in range(len(unq_wyck)):
        if unq_wyck[i][3]==atomtype:
            a = np.array([unq_wyck[i][0], unq_wyck[i][1], unq_wyck[i][2]])
            for j in range(len(structure)):
                b = structure[j].pos
                if all (np.isclose(a, b)):
                    atom_indices.append(j)

    return atom_indices

##########################################################
def get_duplicates(massextr, te_diff=0.01):
    """ function to identify unique intersitials
    
    Parameters
        massextr = Pylada mass extraction object (directory with all intersitials directories)
        te_diff = difference in total energy among interstitials, default = 0.01 eV
    Returns
        dict = {key='foldername', total_energy, index}
    
    Note"
    index = all the foldername (intersitials) with same index are effectively considered equivalent
    Criteria for duplicates
    1. Total_energy_difference <= 0.01 eV
    2. Space_group
    3. comparing first neighbor list for all the sites in the structure
    """

    dict_E = massextr.total_energies
    dict_Str = massextr.structure

    dict2 = {}
    g = 0

    for k in dict_E.keys():
        g = g+1
        dict2[k]=[dict_E[k]]
        dict2[k].append(g)

    folder_list = [l for l in massextr]
    folder_combs = combinations(folder_list, 2)

    for i in folder_combs:
        diff_E = abs(dict2[i[0]][0] - dict2[i[1]][0])
        str1 = dict_Str[i[0]]
        str2 = dict_Str[i[1]]
        spg1 = spglib.get_spacegroup(str1)
        spg2 = spglib.get_spacegroup(str2)
        ngh_list_1 = sorted([len(neighbors(str1, 1, atom.pos, 0.2)) for atom in str1])
        ngh_list_2 = sorted([len(neighbors(str2, 1, atom.pos, 0.2)) for atom in str2])
        if diff_E <= te_diff:
            if spg1 == spg2:
                if ngh_list_1 == ngh_list_2:
                    dict2[i[1]][1] = dict2[i[0]][1]

    return(dict2)

##############################################################
def ffirst_shell(structure, pos, tolerance):
    """ Function to iterate through the first neighbor shell
    
    Parameters
        structure = pylada structure object
        pos = position (numpy array) of site, around which neighbor shell need to be computed
        tolerance = for choosing neighbors whose distance from 'pos' is within first neighbor distance, d(1+tolerance)

    Returns
        list of neighbors, same format as pylada.crystal.neigbors
    """
    
    struct = deepcopy(structure)

    for i, a in enumerate(struct):
        a.index=i
        neighs = [n  for n in neighbors(struct, 12, pos)]
        d = neighs[0][2]
    
    return [n for n in neighs if abs(n[2] - d) < tolerance*d]

##############################################################
def explore_defect(defect, host, tol):
    """ Function to find the position, atom-type etc. of the defect species

    Parameters
        defect = defect structure
        host = host structure
        tol = tolerance to compare defect and host structure

    Returns
        Dictionary containing four items
        1. numpy position vectors of defect site
        2. Atom type
        3. index of defect in defect structure
        4. site index with respect to bulk(host) structure

        Assumptions: 
        1. user has created only single point defect (may work for more than one point defect, but not thoroughly tested)
        2. supercell shape, volume and atom odering in CONTCAR is same in defect and host

        Note:
        1. Adopted from Haowei Peng version in pylada.defects
    """
    
    result = {'interstitial':[], 'substitution':[], 'vacancy': []}
    
    ###look for vacancies
    hstr = deepcopy(host)
    dstr = deepcopy(defect)
    reindex_sites(hstr, dstr, tol)
    
    for atom in hstr:
        if atom.site!=-1: continue
        result['vacancy'].append(deepcopy(atom))
        
    ###look for intersitials and substitutionals
    hstr = deepcopy(host)
    dstr = deepcopy(defect)
    reindex_sites(dstr, hstr, tol)
    
    dstr_copy = deepcopy(dstr)
    for i, atom in enumerate(dstr_copy):
        atom.index = i
        ## Is intersitial always last? No, but reindex_sites, always labels interstitial as site=-1
        if atom.site== -1:
            result['interstitial'].append(deepcopy(atom))
        elif atom.type != hstr[atom.site].type:
            result['substitution'].append(deepcopy(atom))
            
    return result

##############################################################
def avg_electropot(host):
    """ function to compute average electrostatic potential of bulk or defect (DFT calculation)

    Parameters
        pylada.vasp.Extract object

    Returns
        dictionary object with key = atom type, value = avg_electrostatic_potential
    """

    hstr = host.structure

    atmtype = hstr[0].type
    c1 = 0
    dummy_list = []
    
    for n in range(len(hstr)):
        if hstr[n].type == atmtype:
            c1 = c1+1
        else:
            dummy_list.append([atmtype, c1])
            atmtype = hstr[n].type
            c1 = 1
        if n == len(hstr)-1:
            dummy_list.append([atmtype, c1])

    dict_host_e = {}
    dummy = 0

    for i in range(len(dummy_list)):
        avg = np.mean(host.electropot[dummy:dummy+dummy_list[i][1]])
        key = dummy_list[i][0]
        dict_host_e[key]=avg
        dummy = dummy + dummy_list[i][1]

    return dict_host_e

##############################################################
def avg_potential_alignment(defect, host, e_tol=0.2):
    """ function to compute potential alignment correction using average host electrostatic potential
        Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)
    
    Parameters
        defect = pylada.vasp.Extract object
        host = pylada.vasp.Extract object
        e_tol = energy tolerance (same concept to S. Lany Fortran codes)

    Returns
        list = e_tol, #_atoms_excluded, potal_align_corr *eV
        format = [e_tol, nex, pot_align]

    Note:
    1. Some components of the function adopted from Haowei Pengs potential alignment function in pylada.defects
    """

    # get defect structure from pylada Extract object
    d_str = defect.structure

    # list of indicies in defect structure, that are accpetable in pot_align corr
    acceptable = [True for a in d_str]

    #dictionary with average host electrostatic potential with atom type
    dict_hoste = avg_electropot(host)

    #dictionary with average defect electrostatic potential with atom type
    dict_defecte = avg_electropot(defect)

    # compute the difference between defect electrostatic potential and avg. defect pot for each atom
    # Needed to determine unacceptable atoms based on e_tol
    diff_dh = [abs(e - dict_defecte[a.type]).rescale(eV) for e, a in zip(defect.electropot, d_str)]

    # find max value of difference in electrostatic potential of defect and host
    maxdiff = max(diff_dh).magnitude

    # make atoms unacceptable with diff_dh larger than maxdiff or the user energy tolerance(e_tol)
    for ii in range(len(acceptable)):
        if float(diff_dh[ii].magnitude)>= maxdiff or float(diff_dh[ii].magnitude)>=e_tol:
            acceptable[ii] = False

    # check for impurity
    impurity = [k for k in dict_defecte.keys() if k not in set(dict_hoste.keys())]

    for jj in range(len(acceptable)):
        if d_str[jj].type in impurity:
            acceptable[jj] = False

    # count for unacceptable number of atoms in defect structure
    nex = 0

    for jj in range(len(acceptable)):
        if acceptable[jj] == False:
            nex = nex + 1
    
    # Avoid excluding all atoms
    if not any(acceptable):
        # if user give e_tol < 0.0001
        print("ERROR; e_tol is too small excluding all atoms !! Increase e_tol")
    else:
        # compute the difference in electrostatic of defect and host only for acceptable atoms
        diff_dh2 = [(e - dict_hoste[a.type]).rescale(eV) for e, a, ok in zip(defect.electropot, d_str, acceptable) if ok]
        
        # compute alignment_correction
        align_corr = np.mean(diff_dh2)

    return ["{:0.3f}".format(float(e_tol)), nex, "{:0.4f}".format(float(align_corr)), 'eV']

##############################################################
def potential_alignment(defect, host, ngh_shell=False, str_tol=0.4, e_tol=0.2):
    """ function to compute potential alignment correction
        Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)
        
    Parameters
        defect = pylada.vasp.Extract object
        host = pylada.vasp.Extract object
        ngh_shell = bool, choice to let pylada decide to remove defect neighbors from pot_align corr
        str_tol = tolerance to compare host and defect structures
        e_tol = energy tolerance (same concept to S. Lany Fortran codes)

    Returns
        list = str_tol, e_tol, #_atoms_excluded, potal_align_corr *eV
        format = [str_tol, e_tol, nex, pot_align]
        
      Returns average difference of the electrostatic potential of the
      unperturbed atoms in the defect structure with respect to the host.
      *Perturbed* atoms are those flagged as defect by `explore_defect`, their
      first coordination shell if ''first_shell'' is true, and atoms for which
      the electrostatic potential differ to far from the average electrostatic
      potential for each lattice site.

    """
    # get defect structure from pylada Extract object
    d_str = defect.structure
    
    # get host structure from pylada Extract object
    h_str = host.structure
    
    # find defect type, its coordinates, atom-type
    defects = explore_defect(d_str, h_str, str_tol)

    # re-index defect structure w.r.t host structure,
    # function in-built in pylada_crystal_defects
    # CAUTION - may cause problem if defect structure is very highly distorted
    # Works well for split and reasonably (?) distorted structures
    reindex_sites(d_str, h_str, str_tol)
    
    # list of indicies in defect structure, that are accpetable in pot_align corr
    acceptable = [True for a in d_str]
    
    # make intersitial and substitutional unacceptable
    for keys in defects:
        if keys != 'vacancy':
            for atom in defects[keys]:
                acceptable[atom.index] = False

    # Check
    if not any(acceptable):
        # there is a problem: possibly defect and host structure cannot be compared using reindex
        print("ERROR; Cannot compare defect and host !! Switch to avg_potential_alignment")
        sys.exit()

    # make neighbors (upto=tol=0.2) to defects unacceptable 
    if ngh_shell:
        for keys in defects:
            for atom in defects[keys]:
                for n in ffirst_shell(d_str, atom.pos, 0.1):
                    acceptable[n[0].index] = False
                    
    # copy of acceptable before trusting user inputs
    raw_acceptable = deepcopy(acceptable)
    
    # dictionary with average defect elec. pot
    dict_defecte = avg_electropot(defect)

    # list to store abs. differnce in the electrostatic potential of defect and host
    diff_dh = [ (0.0*eV if not ok else abs(e - dict_defecte[a.type]).rescale(eV))
           for e, a, ok in zip(defect.electropot, d_str, acceptable)]
    
    # find max value of difference in electrostatic potential of defect and host
    maxdiff = max(diff_dh).magnitude

    # make atoms with diff_dh larger than user energy tolerance(e_tol) unacceptable
    for ii in range(len(acceptable)):
        if acceptable[ii] == False: pass
        elif float(diff_dh[ii].magnitude)>= maxdiff or float(diff_dh[ii].magnitude)>=e_tol:
            acceptable[ii] = False
    
    # Avoid excluding all atoms
    if not any(acceptable):
        # if user give e_tol < 0.0001
        print("WARNING; e_tol is too small excluding all atoms !! Switching to defaults")
        # return to default, which accepts all atomic sites expect defect sites (and ngh sites)
        acceptable = deepcopy(raw_acceptable)

    # Check
    if not any(acceptable):
        # there is a problem: possibly defect and host structure cannot be compared
        print("ERROR; Cannot compare defect and host !! Switch to avg_potential_alignment")
        sys.exit()

    # count for unacceptable number of atoms in defect structure
    nex = 0
    
    for jj in range(len(acceptable)):
        if acceptable[jj] == False:
            nex = nex + 1
        
    # compute the difference in electrostatic of defect and host only for acceptable atoms
    diff_dh2 = [(e - host.electropot[a.site]).rescale(eV)
           for e, a, ok in zip(defect.electropot, d_str, acceptable) if ok]
    
    # compute alignment_correction
    align_corr = np.mean(diff_dh2)
    
    return ["{:0.3f}".format(float(str_tol)), "{:0.3f}".format(float(e_tol)), nex, "{:0.4f}".format(float(align_corr)), 'eV']   

##############################################################
def get_potential_alignment(defect, host, etol=0.15, etol_range=0.075, steps=10):
    """ Function to compute potential alignment correction for range of e_tol values
        Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)

    Parameters
        defect = pylada.vasp.Extract object
        host = pylada.vasp.Extract object
        etol = energy tolerance
        etol_range = range within which to vary tolerance
        steps = total values in the range

    Returns
        list = str_tol, e_tol, #_atoms_excluded, potal_align_corr *eV
        format = [str_tol, e_tol, nex, pot_align]
    """
    xx = np.linspace(etol-etol_range, etol+etol_range, steps)

    for ii in range(len(xx)):
        potal = potential_alignment(defect, host, e_tol=xx[ii])
        print(potal)

##############################################################
def get_avg_potential_alignment(defect, host, etol=0.15, etol_range=0.075, steps=10):
    """ Function to compute potential alignment correcton using avg. host electrostatic potential, for range of e_tol values
        Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)

    Parameters
        defect = pylada.vasp.Extract object
        host = pylada.vasp.Extract object
        etol = energy tolerance
        tol_range = range within which to vary tolerance
        steps = total values in the range

    Returns
        list = str_tol, e_tol, #_atoms_excluded, potal_align_corr
        format = [str_tol, e_tol, nex, pot_align]
    """
    xx = np.linspace(etol-etol_range, etol+etol_range, steps)

    for ii in range(len(xx)):
        potal = avg_potential_alignment(defect, host, e_tol=xx[ii])
        print(potal)

##############################################################
def get_band_filling(defect, host, potal):
    """ Function to compute band filling correction
    Reference: T.S. Moss, Proc. Phys. Soc. London, Sect. B 67, 75 (1954); E. Burstein, Phys. Rev. 93, 632 (1954)

    Parameters
        host = pylada.vasp.Extract object
        defect = pylada.vasp.Extract object
        potal = potential alignment correction, magnitude only

    Returns
        prints: band-filling correction to the total energy (eV), and band-filling 
    
    Note:
    1.Accounts for the Moss-Burnstein band-filling effect in case of shallow donors and acceptors
    2.Modified the old(or original) Haowei Pengs' version of function(band_filling) in pylada.defects module    
    """

    cbm = host.cbm + potal*eV

    #compute band filling correction to energy (for eigenvalues > cbm)                               
    if defect.eigenvalues.ndim == 3:
        dummy = np.multiply(defect.eigenvalues-cbm, defect.multiplicity[np.newaxis,:,np.newaxis])
        dummy = np.multiply(dummy, defect.occupations)
    elif defect.eigenvalues.ndim == 2:
        dummy = np.multiply(defect.eigenvalues-cbm, defect.multiplicity[:,np.newaxis])
        dummy = np.multiply(dummy, defect.occupations)

    result_n = -np.sum(dummy[defect.eigenvalues > cbm])/ np.sum(defect.multiplicity)

    #compute occupation corresponding to band filling (> cbm)                                        
    if defect.eigenvalues.ndim == 3:
        dummy2_up = np.multiply(defect.occupations[0], defect.multiplicity[:,np.newaxis])
        dummy2_down = np.multiply(defect.occupations[1], defect.multiplicity[:,np.newaxis])
    elif defect.eigenvalues.ndim == 2:
        dummy2 = np.multiply(defect.occupations, defect.multiplicity[:,np.newaxis])

    #compute band filling (> cbm)                                                                    
    if defect.eigenvalues.ndim == 3:
        occ_n_up = np.sum(dummy2_up[defect.eigenvalues[0] > cbm]) / np.sum(defect.multiplicity)
        occ_n_down = np.sum(dummy2_down[defect.eigenvalues[1] > cbm]) / np.sum(defect.multiplicity)
    elif defect.eigenvalues.ndim == 2:
        occ_n = np.sum(dummy2[defect.eigenvalues > cbm]) / np.sum(defect.multiplicity)

    vbm = host.vbm + potal*eV

    #compute band filling correction to energy (for eigenvalues < vbm)                               
    if defect.eigenvalues.ndim == 3:
        dummy = np.multiply(vbm - defect.eigenvalues, defect.multiplicity[np.newaxis,:,np.newaxis])
        dummy = np.multiply(dummy, 1e0 - defect.occupations)
    elif defect.eigenvalues.ndim == 2:
        dummy = np.multiply(vbm - defect.eigenvalues, defect.multiplicity[:,np.newaxis])
        dummy = np.multiply(dummy, 2e0 - defect.occupations)

    result_p = -np.sum(dummy[defect.eigenvalues < vbm])/ np.sum(defect.multiplicity)

    #compute occupation corresponding to band filling (< vbm)                                        
    if defect.eigenvalues.ndim == 3:
        dummy2_up = np.multiply(1e0-defect.occupations[0], defect.multiplicity[:,np.newaxis])
        dummy2_down = np.multiply(1e0-defect.occupations[1], defect.multiplicity[:,np.newaxis])
    elif defect.eigenvalues.ndim == 2:
        dummy2 = np.multiply(2e0-defect.occupations, defect.multiplicity[:,np.newaxis])

    #compute band filling (< vbm)                                                                    
    if defect.eigenvalues.ndim == 3:
        occ_p_up = np.sum(dummy2_up[defect.eigenvalues[0] < vbm]) / np.sum(defect.multiplicity)
        occ_p_down = np.sum(dummy2_down[defect.eigenvalues[1] < vbm]) / np.sum(defect.multiplicity)
    elif defect.eigenvalues.ndim == 2:
        occ_p = np.sum(dummy2[defect.eigenvalues < vbm]) / np.sum(defect.multiplicity)

    if defect.eigenvalues.ndim == 3:
        print("CBM: bf_corr(eV)", "{:0.4f}".format(float(result_n.rescale(eV).magnitude)), 'bf (up, down)', "{:0.4f}".format(float(occ_n_up)), "{:0.4f}".format(float(occ_n_down)))
        print("VBM: bf_corr(eV)", "{:0.4f}".format(float(result_p.rescale(eV).magnitude)), 'bf (up, down)', "{:0.4f}".format(float(occ_p_up)), "{:0.4f}".format(float(occ_p_down)))
    elif defect.eigenvalues.ndim == 2:
        print('CBM: bf_corr(eV)', "{:0.4f}".format(float(result_n.rescale(eV).magnitude)), 'bf', "{:0.4f}".format(float(occ_n)))
        print('VBM: bf_corr(eV)', "{:0.4f}".format(float(result_p.rescale(eV).magnitude)), 'bf', "{:0.4f}".format(float(occ_p)))

##############################################################
def get_madelungenergy(defect, charge=None, epsilon=1e0, cutoff=100.0):
    """ Function returns leading first order correction term, i.e.,
        screened Madelung-like lattice energy of point charge
    Reference: M. Leslie and M. J. Gillan, J. Phys. C: Solid State Phys. 18 (1985) 973

    Parameters
        defect = pylada.vasp.Extract object
        charge = charge of point defect. Default 1e0 elementary charge
        epsilon = dimensionless relative permittivity
        cutoff = Ewald cutoff parameter

    Returns
        Madelung (electrostatic) energy in eV                                                                                                                

    Note:
        1. Units in this function are either handled by the module Quantities, or\
        defaults to Angstrom and elementary charges
        2. Function is adopted from Haowei Peng's version in pylada.defects modules
    """

    from quantities import elementary_charge, eV
    from pylada.crystal import Structure
    from pylada.physics import Ry
    from pylada.ewald import ewald

    if charge is None: charge = 1
    elif charge == 0: return 0e0 *eV
    if hasattr(charge, "units"): charge = float(charge.rescale(elementary_charge))

    ewald_cutoff = cutoff * Ry
    
    structure = defect.structure

    struc = Structure()
    struc.cell = structure.cell
    struc.scale = structure.scale
    struc.add_atom(0., 0., 0., "P", charge=charge)
    
    result = ewald(struc, ewald_cutoff).energy / epsilon
    return -1*result.rescale(eV)

##############################################################
def get_3rdO_corr(defect, charge=None, n=100, epsilon=1e0):
    """ Function returns scaled 3rd order image charge correction
    Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)
               S. Lany and A. Zunger, Model. Simul. Mater. Sci. Eng. 17, 0842002 (2009)
               [Eq. 6, 7 and 11]
                                                                                                                         
    Parameters
        defect = pylada.vasp.Extract object
        charge = charge of point defect. Default 1e0 elementary charge
        epsilon = dimensionless relative permittivity

    Returns
        scaled (1-1/epsilon) third image correction in eV

    Note:
        1. Function is adopted from Haowei Peng's version in pylada.defects modules
    """

    from quantities import elementary_charge, eV, pi, angstrom
    from pylada.physics import a0, Ry

     ## Pgraf's port of the old c-version. (Haowei comments could be faster)
    from pylada.crystal import third_order_cc

    if charge is None: charge = 1e0
    elif charge == 0: return 0e0 *eV
    if hasattr(charge, "units"): charge = float(charge.rescale(elementary_charge))
    if hasattr(epsilon, "units"): epsilon = float(epsilon.simplified)

    structure = defect.structure

    cell = (structure.cell*structure.scale).rescale(a0)

    scaled_3rdO = third_order_cc(cell, n) * (4e0*pi/3e0)* Ry.rescale(eV) * charge * charge \
        *(1e0- 1e0/epsilon) / epsilon

    return scaled_3rdO

##############################################################
def thirdO(defect, charge=None, n=100):
    """ Function returns 3rd order image charge correction, same as LZ fortran script
    Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)
               S. Lany and A. Zunger, Model. Simul. Mater. Sci. Eng. 17, 0842002 (2009)
               [Eq. 6, 7]

    Parameters
        defect = pylada.vasp.Extract object
        charge = charge of point defect. Default 1e0 elementary charge
        n = precision in integral of Eq. 7 (LZ 2009), larger the better

    Returns
        third image correction in eV
    """

    from quantities import elementary_charge, eV, pi, angstrom
    from pylada.physics import a0, Ry

    ## Pgraf's port of the old c-version. Could be faster
    from pylada.crystal import third_order_cc

    if charge is None: charge = 1e0
    elif charge == 0: return 0e0 *eV
    if hasattr(charge, "units"): charge = float(charge.rescale(elementary_charge))

    structure = defect.structure

    cell = (structure.cell*structure.scale).rescale(a0)

    thirdO = third_order_cc(cell, n) * (4e0*pi/3e0) * Ry.rescale(eV) * charge * charge

    return thirdO

##############################################################
def charge_corr(defect, charge=None, epsilon=1e0, cutoff=100., n=100, **kwargs):
    """ Function returns complete image charge correction (Madelung + scaled 3rd Order)
        Reference: S. Lany and A. Zunger, Model. Simul. Mater. Sci. Eng. 17, 0842002 (2009)[Eq. 11]

    Parameters
        defect = pylada.vasp.Extract object
        charge = charge of point defect. Default 1e0 elementary charge
        epsilon = dimensionless relative permittivity
        cutoff = Ewald cutoff parameter
        n = precision in integral of Eq. 7 (LZ 2009), larger the better

    Returns
        Madelung + scaled 3rd order image charge correction in eV

    Note: Adopted from Haowei Peng's version(charge_correction) in pylada.defects                                                           
    """

    if charge is None: charge = 1e0
    elif charge == 0: return 0e0 *eV

    E1 = get_madelungenergy(defect, charge=charge, epsilon=epsilon, cutoff=cutoff)
    E3 = get_3rdO_corr(defect, charge=charge, n=n, epsilon=epsilon)

    result_cc = E1 - E3

    return result_cc

##############################################################
def get_imagecharge(defect, charge=None, epsilon=1e0, cutoff=100., n=100, verbose=False, **kwargs):
    """ Function returns complete image charge correction (Madelung + scaled 3rd Order)
        Reference: S. Lany and A. Zunger, Model. Simul. Mater. Sci. Eng. 17, 0842002 (2009)[Eq. 11]

    Parameters
        defect = pylada.vasp.Extract object
        charge = charge of point defect. Default 1e0 elementary charge
        epsilon = dimensionless relative permittivity
        cutoff = Ewald cutoff parameter
        n = precision in integral of Eq. 7 (LZ 2009), larger the better
        verbose = True or False

    Returns
        Non-verbose = Madelung + scaled 3rd order image charge correction in eV
        Verbose = Madelung_energy, 3rd Order, shape-factor csh, scaling f, final_image_correction in eV
    """

    if charge is None: charge = 1e0
    elif charge == 0: return 0e0 *eV

    E1 = get_madelungenergy(defect, charge=1e0, epsilon=1e0, cutoff=cutoff)
    E3 = -1.*thirdO(defect, charge=1e0, n=n)

    if epsilon == 1e0:
        # epsilon==1e0, meaning vacuum                                                                                 
        print("epsilon=1e0, defect in vacuum. really!!")
        return E1/epsilon
    else:
        scaled_E3 = E3*(1e0 - 1e0/epsilon)
        csh = E3/E1
        f = csh*(1e0 - 1e0/epsilon)
        
        E_ic = (E1 + scaled_E3) * charge * charge /epsilon
        
        if verbose == False:
            return E_ic
        else:
            return ["{:0.3f}".format(float(E1)), "{:0.3f}".format(float(E3)), "{:0.3f}".format(float(csh)) \
                        ,"{:0.3f}".format(float(f)), "{:0.3f}".format(float(E_ic))]

##############################################################
def write_interstitials(structure, ttol=0.5):
    """ function to write POSCAR with interstitials

    Parameters:
        host = pylada structure object, bulk supercell recommended
        ttol = tolerance for finding interstitials

    Return:
        POSCAR file with interstitials, atom type "B"
    """

    struct = deepcopy(structure)
    ints = get_interstitials(struct, ttol)
    test_str = deepcopy(struct)

    for nn in range(len(ints)):
        test_str.add_atom(np.array(ints[nn])[0], np.array(ints[nn])[1], np.array(ints[nn])[2], 'B')

    print('Writing file POSCAR_ints with "B" type atoms as interstitials')

    with open("POSCAR_ints", "w") as file: write.poscar(test_str, file, vasp5=True)

##############################################################
if  __name__=='__main__':
    """ Small code to test the output of above defined functions
    """

#structure = read.poscar('./Example_POSCARS/POSCAR_Si')
#ints = get_interstitials(structure)

#defect = Extract('./Input_files/Si/vac/0')
#host = Extract('./Input_files/Si/host/')

#get_potential_alignment(defect, host)
#get_band_filling(defect, host, potal=0.01)
#get_imagecharge(defect, charge=1., epsilon=13.30)


