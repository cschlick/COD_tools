from rdkit import Chem
from rdkit.Chem import AllChem
from iotbx import cif
from pathlib import Path
import copy
import os
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, cdist
import io
import itertools
import copy
import tqdm
import numpy as np
from io import StringIO
from contextlib import redirect_stderr
from rdkit import Chem
import copy
from collections import defaultdict

from guess_bond_order import guess_bond_order_mda
from rdkit_utils import validate_rdkit_mol

class CODFile:
  """
  A class to parse, filter, and process molecules from a COD database file.
  
  It moves through a sequence of filters. If everything passes, the end result are
  sanitized rdkit molecules.
  
  """
  @staticmethod
  def subtract_paths(path1,path2):
      p1 = list(path1.parts)
      p2 = list(path2.parts)
      for p in p1:
        p2.remove(p)
      return Path(*tuple(p2))
    
  # failure messages in the order they will be checked
  default_failure_messages = {
  -1:"Unknown failure",
  0:"Error parsing smiles from html file",
  1:"File does not exist",
  2:"IOTBX cif file parsing error",
  3:"Is powder diffraction",
  4:"Has no diffraction reflections",
  5:"Has no resolution information",
  6:"Resolution falls below threshold",
  7:"Has no R-factor information",
  8:"R-factor above threshold",
  9:"Failure building CCTBX xray-structure",
  10:"Contains disallowed elements",
  10.5:"Unable to extract coordinates",
  11:"Number of atoms falls below threshold ",
  11.5:"Number of atoms falls above threshold ",
  12:"Atoms pairs clash within threshold ",
  13:"No explicit bond information",
  14:"Mismatched atom/bond labels",
  15:"Occupancy values exist below threshold ",
  16:"Failure linking hydrogen atoms",
  16.5:"Hydrogen has >1 bonded partner",
  17:"Bond accross symmetry",
  18:"Failure building molecule as a networkx graph",
  19:"Failure building molecule as an rdkit molecule object",
  19.5:"Failure in preliminary rdkit sanitization",
  20:"Warning when building molecule as rdkit",
  23:"Failed assigning bond orders",
  24:"Failed in final rdkit sanitization"}
  
  
  def __init__(self,filepath,
               cod_id=None,
               debug=False,
               allowed_elements=[],
               dmax_thresh=1.0,
               rfactor_thresh=0.1,
               occ_thresh=1.0,
               min_dist_thresh=0.5,
               min_num_thresh=3,
               max_num_thresh=200):
    
    filepath = Path(filepath)
    self.allowed_elements = allowed_elements
    if cod_id == None:
      cod_id = filepath.name
    self.cod_id = cod_id
    self.filepath = filepath
    self.warnings = [] # a list of error/warning strings
    pt  = Chem.GetPeriodicTable()
    self.allowed_elements = set([e for e in self.allowed_elements])
    
    self.failure_messages = copy.copy(self.default_failure_messages)
    
    self.dmax_thresh = dmax_thresh
    self.failure_messages[6] = "Resolution falls below threshold "+"(%s)" %str(self.dmax_thresh)
    
    self.occ_thresh = occ_thresh
    self.failure_messages[15] = "Occupancy values exist below threshold "+"(%s)" %str(self.occ_thresh)
    
    self.rfactor_thresh = rfactor_thresh
    self.failure_messages[8] = "R-factor above threshold "+"(%s)" %str(self.rfactor_thresh)
    
    self.min_dist_thresh = min_dist_thresh
    self.failure_messages[12] = self.failure_messages[12]+"(%s)" %str(self.min_dist_thresh)
    
    self.min_num_thresh = min_num_thresh
    self.failure_messages[11] = self.failure_messages[11]+"(%s)" %str(self.min_num_thresh)
    
    self.max_num_thresh = max_num_thresh
    self.failure_messages[11.5] = self.failure_messages[11.5]+"(%s)" %str(self.max_num_thresh)
    failed = False
    fail_message = None
    rdkit_error = None
    

    ########### Check if file exists
    
    if not self.filepath.exists():
      failed = True
      fail_message = self.failure_messages[1]
      
      
    ########### Parse file
    if not failed:
      try:
          cif_model = cif.reader(self.filepath.as_posix()).model()
          code = cif_model.keys()[0]
          self.cif_model = cif_model[code]
          cif_key_set= set(self.cif_model.keys())

      except:
        self.cif_model = None
        failed = True
        fail_message = self.failure_messages[2]
        
    
    
    
    ########### Check Powder diffraction
    if not failed:
      
      if any(["_pd_" in key for key in cif_key_set]):
        failed = True
        fail_message = self.failure_messages[3]
    
    
    ########### Check Powder diffraction
    
    if not failed:
      if not any(["_diffrn_reflns_" in key for key in cif_key_set]):
        failed=True
        fail_message = self.failure_messages[4]
    
    ########### Check for Resolution info
    
    if not failed:
      resolution_keys = {"_diffrn_radiation_wavelength","_diffrn_reflns_theta_max"}
      if len(cif_key_set.intersection(resolution_keys)) !=2:
        failed = True
        fail_message = self.failure_messages[5]
    
    ########### Check Resolution info
    if not failed:
      
      try:
        wavelength_str = self.cif_model["_diffrn_radiation_wavelength"]
        theta_max_str = self.cif_model["_diffrn_reflns_theta_max"]
        wavelength = float(np.array(wavelength_str))
        theta_max = float(np.array(theta_max_str))
        sin_theta = np.sin(theta_max)
        if np.isclose(sin_theta,0) or not np.isfinite(sin_theta) or not np.isfinite(wavelength):
          failed = True
          fail_message =self.failure_messages[5]
        if not failed:
          dmax = wavelength/(2*np.sin(theta_max))
          if dmax > self.dmax_thresh:
            failed = True
            fail_message = self.failure_messages[6]
      except:
        failed = True
        fail_message = self.failure_messages[5]

        
    ########### Check R factor info
    if not failed:
      Rfactor_cif_keys = set(["_refine_ls_R_factor_all","_refine_ls_R_factor_gt"])

      if cif_key_set.intersection(Rfactor_cif_keys)==0:
        failed = True
        fail_message = self.failure_messages[7]
    
    
    ########### Check R factor values
    if not failed:   
      R_values = []
      if "_refine_ls_R_factor_gt" in self.cif_model.keys():
        R_values.append(self.cif_model["_refine_ls_R_factor_gt"])
      if "_refine_ls_R_factor_all" in self.cif_model.keys():
        R_values.append(self.cif_model["_refine_ls_R_factor_all"])
      try:
        R_values = np.array(R_values,dtype=float)

        if np.any(R_values>self.rfactor_thresh):
          failed=True
          fail_message = self.failure_messages[8]

      except:
        failed=True
        fail_message = self.failure_messages[8]

    
    ########### Try building xray structure
    if not failed:
      try:
        self.xs = cif.builders.crystal_structure_builder(self.cif_model).structure
      except:
        self.xs = None
        failed = True
        fail_message = self.failure_messages[9]

    
    ########### Reject elements
    if not failed:
      try:
        pt  = Chem.GetPeriodicTable()
        self.atomic_names = ["H" if sc.element_symbol()=="D" else sc.element_symbol() for sc in self.xs.scatterers()]
        self.atomic_numbers = [pt.GetAtomicNumber(element) for element in self.atomic_names]
        self.element_set = set(self.atomic_names)
        if len(self.allowed_elements)>0:
          if not self.element_set.issubset(self.allowed_elements):
            failed = True
            fail_message = self.failure_messages[10]

      except:
        failed = True
        fail_message = self.failure_messages[10]
    
    
    ########### Get coords
    if not failed:
      try:
        self.xyz = self.xs.sites_cart().as_numpy_array()
      except:
        failed = True
        fail_message = self.failure_messages[10.5]

        
    ########### Check minimum size
    if not failed:
      if len(self.atomic_numbers)<=self.min_num_thresh:
        failed = True
        fail_message = self.failure_messages[11]
    
    ########### Check maximum size
    if not failed:
      if len(self.atomic_numbers)>=self.max_num_thresh:
        failed = True
        fail_message = self.failure_messages[11.5]
    
    
    ########### Check minimum atom pair distance
    if not failed:
      dists = pdist(self.xyz)
      if len(dists) <1 or not np.all(np.isfinite(dists)):
        failed = True
        fail_message = self.failure_messages[12]
      if not failed:
        if np.min(dists)<=self.min_dist_thresh:
          failed = True
          fail_message = self.failure_messages[12]

    ########### Check for explicit bonds
    if not failed:
      if "_geom_bond_atom_site_label_1" not in self.cif_model.keys() or "_geom_bond_atom_site_label_2" not in self.cif_model.keys():
        failed = True
        fail_message = self.failure_messages[13]

    ########### Check mismatched atom/bond labels
    if not failed:
      cif_model = self.cif_model
      site_labels = list(cif_model["_atom_site_label"])
      labels = {label for pair in zip(cif_model["_geom_bond_atom_site_label_1"],cif_model["_geom_bond_atom_site_label_2"]) for label in pair}
      if not labels.issubset(set(site_labels)):
        failed = True
        fail_message = self.failure_messages[14]

    ########### Check low occupancy
    if not failed:
      occupancy = np.array([sc.occupancy for sc in self.xs.scatterers()])
      if np.any(occupancy<self.occ_thresh):
        failed = True
        fail_message = self.failure_messages[15]
       
    

    ########### check for bonds accross symmetry
    if not failed:
      key = '_geom_bond_site_symmetry_2'
      try:
        bond_sym = list(self.cif_model[key])
        if len(list(set(bond_sym)))>1:
          failed = True
          fail_message = self.failure_messages[17]
      except:
        pass
    
        
    ########### Build nx graph
    if not failed:
      try:
        cif_model = self.cif_model
        site_labels = list(cif_model["_atom_site_label"])
        label_dict_label_keys = {label:i for i,label in enumerate(site_labels)}
        label_dict_idx_keys = {i:label for i,label in enumerate(site_labels)}
        xs = self.xs
        G = nx.Graph()
        for i,(sc,xyz) in enumerate(zip(xs.scatterers(),xs.sites_cart())):
          G.add_node(i)

        for label1,label2 in zip(cif_model["_geom_bond_atom_site_label_1"],cif_model["_geom_bond_atom_site_label_2"]):
          idx1,idx2 = label_dict_label_keys[label1], label_dict_label_keys[label2]
          G.add_edge(idx1,idx2)

        self.atom_graph = G
      except:
        failed = True
        fail_message = self.failure_messages[18]  
           
    ########### check if H bonds are explicit
    if not failed:
      try:
        non_explicit_Hbonds = []
        excessive_valence_hbonds = []
        for atom in self.atom_graph:
          atomic_number = self.atomic_numbers[atom]
          if atomic_number==1:
            if len(self.atom_graph.edges(atom)) == 0:
              non_explicit_Hbonds.append(atom)
            elif len(self.atom_graph.edges(atom)) >1:
              excessive_valence_hbonds.append(atom)
        
        if len(excessive_valence_hbonds)>0:
          failed = True
          fail_message = self.failure_messages[16.5]
        
        if not failed:
          if len(non_explicit_Hbonds)>0:
            dist = cdist(self.xyz,self.xyz)
            for atomi in non_explicit_Hbonds:
              partner = np.argwhere(dist[atomi]==np.sort(dist[atomi])[1])[0][0]
              self.atom_graph.add_edge(atomi,partner)
      except:
        failed = True
        fail_message = self.failure_messages[16]     
      

      ########### Build RDKIT Molecules (1 for each connected component subgraph)      
        
      if not failed:
        try:
          pt  = Chem.GetPeriodicTable()
          G = self.atom_graph
          xs = self.xs    
          atom_ids = list(self.cif_model["_atom_site_label"])
          mols = []
          mol_infos = []
          pdb_strings = []
          for subgraph in list(nx.connected_components(G)):
            mol = Chem.Mol()
            mol_info = defaultdict(list)
            rwmol = Chem.RWMol(mol)
            conformer = Chem.Conformer(len(subgraph))

            new_idx_mapping = {} # maps index in file to index in molecule (subgraph)
            bond_pairs = set()

            for i_new,i_old in enumerate(subgraph):
              new_idx_mapping[i_old]=i_new
              sc = xs.scatterers()[i_old]
              xyz = xs.sites_cart()[i_old]
              sc = ("H" if sc.element_symbol()=="D" else sc.element_symbol())
              atom = Chem.Atom(pt.GetAtomicNumber(sc))
              atom.SetProp("atom_id",atom_ids[i_old])
              atom.SetProp("file_atom_index",str(i_old))
              mol_info["atom_id"].append(atom_ids[i_old])
              mol_info["file_atom_index"].append(str(i_old))
              atom = rwmol.AddAtom(atom)
              conformer.SetAtomPosition(i_new,xyz)

            for i_old in subgraph:
              for edge in G.edges(i_old):
                start,end = edge
                if start!=end:
                  start, end = new_idx_mapping[start], new_idx_mapping[end]
                  bond_pairs.add(tuple(sorted([start,end])))

            for bond_pair in bond_pairs:
              idx1,idx2 = bond_pair
              rwmol.AddBond(idx1,idx2,Chem.BondType.SINGLE)

            rwmol.AddConformer(conformer)  
            mol = rwmol.GetMol()

            mols.append(mol)
            mol_infos.append(mol_info)
          self.mols = mols
          self.mol_infos = mol_infos
        except:
          failed = True
          fail_message = self.failure_messages[19]
          if debug:
            raise
            
            
    ########### Check rdkit issues
        


#     def test_bond_order(mol):
#       """
#       Does some checks to see if the mol will fail later.
#       Explicit valence issues for example.
#       Returns:
#         (failed,error) where failed is bool and error is the stderr string
#       """
#       with redirect_stderr(StringIO()) as err:
#         try:
#           frags = Chem.GetMolFrags(mol, asMols=True)
#           return True,err.getvalue()
#         except:
#           return False,err.getvalue()
        
#     if not failed:
#       for mol in self.mols:
#         f,e = test_bond_order(mol)
#         if not f:
#           failed = True
#           fail_message = self.failure_messages[19.5]
#           rdkit_error = e
#           break
    
    
    ########### Assign bond order     
        
    if not failed:
      with redirect_stderr(StringIO()) as err:
        try:
          fixed_mols = []
          for mol in self.mols:
            newmol = guess_bond_order_mda(mol,debug=False)
            fixed_mols.append(newmol)
          self.fixed_mols = fixed_mols
          if self.fixed_mols.count(None)>0:
            failed = True
            fail_message = self.failure_messages[23]
        except:
          failed = True
          fail_message = self.failure_messages[23]
      

    ########### Final Sanitization 

    if not failed:
      #try:
      for mol in self.fixed_mols:
        f,e = validate_rdkit_mol(mol,debug=debug,return_err=True)
        rdkit_error = e
        if not f:
          failed = True
          fail_message = self.failure_messages[24]
          break
      # except:
      #   failed = True
      #   fail_message = self.failure_messages[24]
      #   if debug:
      #     raise
      
    ######## END
    self.failed = failed
    self.fail_message = fail_message
    self.rdkit_error = rdkit_error
    if failed:
      self.mols = []
      self.fixed_mols = []
    