from io import StringIO
from contextlib import redirect_stderr
from rdkit import Chem

def validate_rdkit_mol(rdkit_mol,roundtrip=False,return_err=False,debug=False):
  """
  Sanitize, check result, and optionally do
  round trip .mol conversion
  
  Will return False if:
    1. Chem.SanitizeMol() returns an exception
    2. Chem.SanitizeMol() succeeds but returns non-standard value
    3. Any text is sent to stderr
  """
  #import rdkit
  #rdkit.rdBase.LogToPythonStderr()
  with redirect_stderr(StringIO()) as err:
    try:
      # sanitize
      sanitize_ret1 = Chem.SanitizeMol(rdkit_mol)
      assert sanitize_ret1 == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
      if roundtrip:
        # round trip conversion
        mblock = Chem.MolToMolBlock(rdkit_mol)
        m = Chem.MolFromMolBlock(mblock)
        sanitize_ret2= Chem.SanitizeMol(m)
        assert sanitize_ret2 == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
      
      # Case1: everything ok
      if len(err.getvalue())==0: 
        if return_err:
          return True, None
        return True
      
      # Case2: No errors, but warnings
      else:
        if debug:
          raise
        if return_err:
          return False, err.getvalue()
        return False
    
    # Case 3: Errors
    except:
      if debug:
        raise
      if return_err:
        return False, ""
      return False
  