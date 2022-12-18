import copy
import numpy as np
import MDAnalysis as mda
import tempfile
from io import StringIO
from contextlib import redirect_stderr
from rdkit import Chem
import copy

def guess_bond_order_mda(rdkit_mol,debug=False,reterr=False):
  """
  Use MDAnalysis to guess bond orders
  Uses tempfile intermediate.
  """
  with redirect_stderr(StringIO()) as err:
    try:
      m = copy.deepcopy(rdkit_mol)
      params = Chem.rdmolops.AdjustQueryParameters()
      params.makeBondsGeneric = True
      modmol = Chem.rdmolops.AdjustQueryProperties(m, params)
      modmol.UpdatePropertyCache()
      pdb_string= Chem.MolToPDBBlock(modmol)
      with tempfile.NamedTemporaryFile("w+t",suffix=".pdb") as fh:
        fh.write(pdb_string)
        fh.seek(0)
        universe= mda.Universe(fh.name,format="PDB")
        fh.close()
      #elements = mda.topology.guessers.guess_types(universe.atoms.names)
      elements = [atom.GetSymbol().upper() for atom in m.GetAtoms()]
      universe.add_TopologyAttr('elements', elements)
      element_set = set([a.GetSymbol() for a in m.GetAtoms()])
      #force = ("H" not in element_set)
      force = True
      newmol = universe.atoms.convert_to("RDKIT",force=force)
      for i in range(m.GetNumAtoms()):
        atom_old = rdkit_mol.GetAtomWithIdx(i)
        atom_new = newmol.GetAtomWithIdx(i)
        if atom_old.HasProp("atom_id"):
          atom_new.SetProp("atom_id",atom_old.GetProp("atom_id"))
        if atom_old.HasProp("file_atom_index"):
          atom_new.SetProp("file_atom_index",atom_old.GetProp("file_atom_index"))
    except:
      if debug:
        raise
      newmol = None
    if not reterr:
      return newmol
    else:
      return newmol,err.getvalue()