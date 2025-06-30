import numpy as np
import os
import re
import tempfile
import importlib.util
import sys

class Character_table(object):
    def __init__(self, point_group):
        point_group = point_group.lower()   # make it case insensitive
        tables = {
            "cs": {
                "classes": ["E", "σh"],
                "irreps": {
                    "A'":  [1, 1],
                    "A''": [1,-1],
                }
            },
            "c2": {
                "classes": ["E", "C2"],
                "irreps": {
                    "A'": [1, 1],
                    "B":  [1,-1],
                }
            },
            "c2v": {
                "classes": ["E", "C2", "σv(xz)", "σv(yz)"],
                "irreps": {
                    "A1": [1, 1, 1, 1],
                    "A2": [1, 1,-1,-1],
                    "B1": [1,-1, 1,-1],
                    "B2": [1,-1,-1, 1],
                }
            },
            "c2h": {
                "classes": ["E", "C2", "i", "σh"],
                "irreps": {
                    "Ag": [1, 1, 1, 1],
                    "Bg": [1,-1, 1,-1],
                    "Au": [1, 1,-1,-1],
                    "Bu": [1,-1,-1, 1],
                }
            },
            "d2h": {
                "classes": ["E", "C2(z)", "C2(y)","C2(x)","i", "σ(xy)", "σ(xz)", "σ(yz)"],
                "irreps": {
                    "Ag":  [1, 1, 1, 1, 1, 1, 1, 1],
                    "B1g": [1, 1,-1,-1, 1, 1,-1,-1],
                    "B2g": [1,-1, 1,-1, 1,-1, 1,-1],
                    "B3g": [1,-1,-1, 1, 1,-1,-1, 1],
                    "Au":  [1, 1, 1, 1,-1,-1,-1,-1],
                    "B1u": [1, 1,-1,-1,-1,-1, 1, 1],
                    "B2u": [1,-1, 1,-1,-1, 1,-1, 1],
                    "B3u": [1,-1,-1, 1,-1, 1, 1,-1]
                }
            }
        }
        
        # check if character table is present
        if point_group not in tables:
            raise ValueError (f"Point group '{point_group}' not found.")
        
        self.point_group = point_group
        self.classes = tables[point_group]["classes"]
        self.irreps = tables[point_group]["irreps"]
    
        #print(f"Symmetry Operations: {self.classes}")
        #print(f"irreps: {self.irreps}")

    @staticmethod
    def normalize_input(input_str):
        """
        Normalize common inputs: sigma → σ, prime → ', underscores removed.
        """
        input_str = input_str.lower().strip()
        input_str = input_str.replace("sigma", "σ")
        input_str = input_str.replace("prime", "'")
        input_str = input_str.replace("_", "")
        return input_str
    
    def get_matrix(self, symmetry_operation):
        """
        Return 3x3 transformation matrix for a named symmetry operation
        """
        symmetry_operation = self.normalize_input(symmetry_operation)

        # check if the operation exists
        normalized_op = [self.normalize_input(cls) for cls in self.classes]
        if symmetry_operation not in normalized_op:
            raise ValueError(f"Operation '{symmetry_operation}' is not in the character table of point group '{self.point_group}'.")

        # Use the original operation name as written in self.classes
        original_op = self.classes[normalized_op.index(symmetry_operation)]

        # Identity
        if original_op == "E":
            return np.identity(3)
        
        # Inversion
        if original_op == "i":
            return -np.identity(3)
        
        # Reflection
        if original_op in ["σh", "σ(xy)"]:
            return np.diag([1, 1, -1])
        if original_op in ["σv(xz)", "σ(xz)"]:
            return np.diag([1, -1, 1])
        if original_op in ["σv(yz)", "σ(yz)"]:
            return np.diag([-1, 1, 1])
        
        # Proper and improper rotation
        match_rot = re.fullmatch(r'(c|s)(\d+)(?:\((x|y|z)\))?', original_op, re.IGNORECASE)
        if match_rot:
            kind = match_rot.group(1)
            n = int(match_rot.group(2))
            axis = match_rot.group(3) or "z"

            # sanity check
            if n < 1:
                raise ValueError(f"n must be ≥ 1")
            if kind == 'S' and axis != 'z':
                raise ValueError(f"Sn improper rotations must be about the z-axis; got axis '{axis}'")

            theta = 2 * np.pi / n
            c, s = np.cos(theta), np.sin(theta)

            if axis == 'x':
                rot = np.array([
                    [1, 0,  0],
                    [0, c, -s],
                    [0, s,  c]
                ])
            elif axis == 'y':
                rot = np.array([
                    [ c, 0, s],
                    [ 0, 1, 0],
                    [-s, 0, c]
                ])
            elif axis == 'z':
                rot = np.array([
                    [c, -s, 0],
                    [s,  c, 0],
                    [0,  0, 1]
                ])
            else:
                raise ValueError(f"Invalid axis: '{axis}'")

            return rot if kind == 'C' else np.diag([1, 1, -1]) @ rot
        
        raise ValueError(f"Symmetry operation '{symmetry_operation}' not recognized.")
# print(Character_table("d2h").get_matrix("c2(z)"))

class Zmat(object):
    def __init__(self, base_dir=None):
        if base_dir is None:   # if base_dir is not provided, use where this script is located
            base_dir = os.path.dirname(os.path.abspath(__file__))
        self.base_dir = base_dir

    def read_zmat(self, index: str):
        """
        With the given molecule label, such as "4.01", find and extract information from the corresponding zmat file
        Returns: zmat_list, cartesian_list
        """
        try:
            test_set, molecule = index.split(".")
        except ValueError:
            raise ValueError("Index must be in the format 'X.YY'")
        
        # locate the test set directory
        test_set_dir = None
        for entry in os.listdir(self.base_dir):
            if entry.startswith(test_set + "_") and os.path.isdir(os.path.join(self.base_dir, entry)):
                test_set_dir = entry
                break
        if test_set_dir is None:
            raise FileNotFoundError(f"No directory starting with '{test_set}_' found in {self.base_dir}")
        test_set_path = os.path.join(self.base_dir, test_set_dir)

        # locate the molecule directory
        molecule_dir = None
        for entry in os.listdir(test_set_path):
            if re.match(rf'^{re.escape(molecule)}\D', entry) and os.path.isdir(os.path.join(test_set_path, entry)):
                molecule_dir = entry
                break
        if molecule_dir is None:
            raise FileNotFoundError(f"No molecule directory starting with '{molecule}' found in {test_set_path}")
        molecule_path = os.path.join(test_set_path, molecule_dir)

        # locate "zmat" file inside CCSD_T_TZ
        self.zmat_path = os.path.join(molecule_path, "CCSD_T_TZ", "zmat")
        if not os.path.isfile(self.zmat_path):
            raise FileNotFoundError(f"'zmat' file not found at expected location: {self.zmat_path}")
        
        # extract ZMAT into a list in list
        with open(self.zmat_path, "r") as zmat_txt:
            txt = zmat_txt.readlines()
        
        self.zmat_list = []
        in_zmat = False
        for line in txt:
            if "ZMAT begin" in line:
                in_zmat = True
                continue
            if "ZMAT end" in line:
                break
            if in_zmat:
                stripped = line.strip()
                if stripped:  # if stripped line contains string
                    self.zmat_list.append(stripped.split())

        # extract Cartesian coordinates (neglecting labels)
        self.cartesian_list = []
        in_cart = False
        for line in txt:
            if "cart begin" in line:
                in_cart = True
                continue
            if "cart end" in line:
                break
            if in_cart:
                stripped = line.strip()
                if stripped:
                    self.cartesian_list.append(list(map(float,line.strip().split()[1:4])))

        return self.zmat_list, self.cartesian_list
#print(Zmat().read_zmat("4.01"))

class Manual_projection(object):
    def __init__(self, base_dir=None):
        if base_dir is None:   # if base_dir is not provided, use where this script is located
            base_dir = os.path.dirname(os.path.abspath(__file__))
        self.base_dir = base_dir
    
    def get_proj(self, index: str):
        """
        With the given molecule label, such as "4.01",
        """
        try:
            test_set, molecule = index.split(".")
        except ValueError:
            raise ValueError("Index must be in the format 'X.YY'")
        
        # locate the test set directory
        test_set_dir = None
        for entry in os.listdir(self.base_dir):
            if entry.startswith(test_set + "_") and os.path.isdir(os.path.join(self.base_dir, entry)):
                test_set_dir = entry
                break
        if test_set_dir is None:
            raise FileNotFoundError(f"No directory starting with '{test_set}_' found in {self.base_dir}")
        test_set_path = os.path.join(self.base_dir, test_set_dir)

        # locate the molecule directory
        molecule_dir = None
        for entry in os.listdir(test_set_path):
            if re.match(rf'^{re.escape(molecule)}\D', entry) and os.path.isdir(os.path.join(test_set_path, entry)):
                molecule_dir = entry
                break
        if molecule_dir is None:
            raise FileNotFoundError(f"No molecule directory starting with '{molecule}' found in {test_set_path}")
        molecule_path = os.path.join(test_set_path, molecule_dir)

        # locate manual_projection.py
        self.manual_projection_path = os.path.join(molecule_path, "manual_projection.py")
        if not os.path.isfile(self.manual_projection_path):
            raise FileNotFoundError(f"'manual_projection.py' file not found at expected location: {self.manual_projection_path}")

        # read manual_projection.py and modify it
        with open(self.manual_projection_path, "r") as f:
            lines = f.readlines()

        modified_lines = []
        in_run_function = False

        for line in lines:
            stripped = line.lstrip()

            # detect start of run()
            if stripped.startswith("def run"):
                in_run_function = True
                modified_lines.append(line)
                continue
            # detect end
            elif in_run_function and stripped.startswith("self.Proj = Proj"):
                in_run_function = False
                modified_lines.append(line)
                continue
            # if in run(), get rid of normalize() call
            elif in_run_function:
                line = line.replace("normalize(", "")
                line = line.replace(".T)", "")
                modified_lines.append(line)
            else:
                modified_lines.append(line)

        # write to a temporary file
        tmp_dir = tempfile.gettempdir()
        tmp_path = os.path.join(tmp_dir, "manual_projection_modified.py")
        with open(tmp_path, "w") as f:
            f.writelines(modified_lines)

        # Dynamically import manual_projection.py
        spec = importlib.util.spec_from_file_location("manual_projection", tmp_path)
        manual_proj_module = importlib.util.module_from_spec(spec)
        sys.modules["manual_projection"] = manual_proj_module   # register
        spec.loader.exec_module(manual_proj_module)

        # Instantiate Projection
        proj_instance = manual_proj_module.Projection([])
        proj_instance.run()

        return proj_instance.Proj

#np.set_printoptions(linewidth=np.inf)   # don't wrap
#print(Manual_projection().get_proj("4.01"))

class Irrep(object):
    def __init__ (self, molecule, point_group):
        """
        example: Irrep("4.21", "c2v")
        """
        self.zmat = Zmat().read_zmat(molecule)[0]
        self.cart = Zmat().read_zmat(molecule)[1]
        self.proj = Manual_projection().get_proj(molecule)
        self.character_table = Character_table(point_group)
    
    def assign_irrep(self, tol=1e-5):
        
        # create a list of list to be filled with characters for each NIC later
        nic_chars = [ [] for _ in range(self.proj.shape[0]) ]

        N = len(self.cart)   # # of atoms
        # loop over each symmetry operation
        for op in self.character_table.classes:

            matrix = self.character_table.get_matrix(op)
            transformed = self.cart @ matrix.T   # transformed cartesian

            # generate atom permutation matrix
            P_atom = np.zeros((N, N), dtype = int)   # permutation matrix
            # for each transformed atom
            for i, trans_atom in enumerate(transformed):
                # compute distances to all original atoms
                distances = np.linalg.norm(self.cart - trans_atom, axis=1)
                # find closest original atom
                j = np.argmin(distances)            
                # check that the match is valid
                if distances[j] > tol:
                    raise ValueError(f"atom {i+1} could not be matched to original atom (distance {distances[j]} > tol={tol}) under operation {op}")
                # mark in permutation matrix
                P_atom[i, j] = 1
            
            # generate zmat permutation matrix
            perm_atom = np.argmax(P_atom, axis=1)   # atom permutation array

            """ DEBUG: check how atoms moved """
            # if op == "E":
            #     print(f"Permutation array under operation {op}:")
            #     print(matrix)
            #     print(perm_atom + 1)

            P_zmat = np.zeros((len(self.zmat), len(self.zmat)), dtype=int)
            for col, int_coord in enumerate(self.zmat):
                trans_int_coord = []
                for atom in int_coord:
                    if atom.isdigit():
                        trans_atom = str(perm_atom[int(atom) - 1] + 1)
                        trans_int_coord.append(trans_atom)
                    else:
                        trans_int_coord.append(atom)

                """ DEBUG: compare internal coords before and after transformation """
                # if op == "E":
                #     print(f"Compare internal coords under operation: {op}")
                #     print(f"before: {int_coord}")
                #     print(f"after:  {trans_int_coord}")

                # deal with atoms ordering and phase factor of int_coord
                matched = False
                for row, old_int_coord in enumerate(self.zmat):
                    # stretch
                    if len(old_int_coord) == 2 and len(trans_int_coord) == 2:
                        if set(old_int_coord) == set(trans_int_coord):
                            P_zmat[row, col] = 1
                            matched = True
                            break
                    # bend
                    elif len(old_int_coord) == 3 and len(trans_int_coord) == 3:
                        if old_int_coord == trans_int_coord or old_int_coord == trans_int_coord[::-1]:
                            P_zmat[row, col] = 1
                            matched = True
                            break
                    # torsion
                    elif len(old_int_coord) == 5 and len(trans_int_coord) == 5 and old_int_coord[-1] == 'T' and trans_int_coord[-1] == 'T':
                        if set(old_int_coord[1:3]) == set(trans_int_coord[1:3]) and (old_int_coord[0] == trans_int_coord[0] or old_int_coord[0] == trans_int_coord[3]):
                            P_zmat[row, col] = 1
                            if old_int_coord[1:3] == trans_int_coord[1:3][::-1]:
                                P_zmat[row, col] *= -1
                            if old_int_coord[0] == trans_int_coord[3] or old_int_coord[3] == trans_int_coord[0]:
                                P_zmat[row, col] *= -1
                            if "σ" in op or "S" in op or op == "i":
                                P_zmat[row, col] *= -1
                                matched = True
                                break
                            else:
                                matched = True
                                break
                    # out-of-plane
                    elif len(old_int_coord) == 5 and len(trans_int_coord) == 5 and old_int_coord[-1] == 'O' and trans_int_coord[-1] == 'O':
                        if old_int_coord == trans_int_coord:
                            P_zmat[row, col] = 1
                            matched = True
                        elif old_int_coord[0:2] == trans_int_coord[0:2] and old_int_coord[2:4] == trans_int_coord[3:1:-1]:
                            P_zmat[row, col] = -1
                            matched = True
                        # if operation involves reflection flip the sign
                        if matched:
                            if "σ" in op or "S" in op or op == "i":
                                P_zmat[row, col] *= -1
                                break
                            else:
                                break
                    # Lx
                    elif len(old_int_coord) == 5 and len(trans_int_coord) == 5 and old_int_coord[-1] == 'Lx' and trans_int_coord[-1] == 'Lx':
                        if old_int_coord == trans_int_coord:
                            P_zmat[row, col] = 1
                            matched = True
                            break
                    # Ly
                    elif len(old_int_coord) == 5 and len(trans_int_coord) == 5 and old_int_coord[-1] == 'Ly' and trans_int_coord[-1] == 'Ly':
                        if old_int_coord == trans_int_coord:
                            P_zmat[row, col] = 1
                            matched = True
                            if matched:
                                if "σ" in op or "S" in op or op == "i":
                                    P_zmat[row, col] *= -1
                                    break
                                else:
                                    break
                if not matched:
                    raise ValueError(f"Could not match transformed ZMAT entry {trans_int_coord} under operation {op}")
            """ DEBUG: print P_zmat """
            # print(f"Zmat permutation matrix for operation: {op}")
            # print(P_zmat)

            # transform projection matrix via P_zmat
            trans_proj = (P_zmat @ self.proj.T).T

            """ DEBUG: compare projection matrix before and after the transformation """
            # proj_idx = 29
            # print(f"Symmetry Operation: {op}")
            # print(self.proj[proj_idx])
            # print(trans_proj[proj_idx])

            # compare projection matrices before and after the transformation
            for i, (orig_nic, trans_nic) in enumerate(zip(self.proj, trans_proj)):
                with np.errstate(divide='ignore', invalid='ignore'):
                    quot_list = np.true_divide(trans_nic, orig_nic)   # divide element by element
                    quot_list = quot_list[~np.isnan(quot_list) & ~np.isinf(quot_list)]   # filter out nan or 0/0
                
                # sanity check
                if len(quot_list) > 0 and np.allclose(quot_list, quot_list[0], atol=1e-5):
                    char = quot_list[0]
                    #print(quot_list)
                else:
                    raise ValueError(
                        f"NIC #{i} does not transform consistently under operation '{op}'. \n"
                        f"Original:    {orig_nic}\n"
                        f"Transformed: {trans_nic}\n"
                    )
                nic_chars[i].append(char)

        # compare to character table and assign irreps
        character_table_syms = list(self.character_table.irreps.keys())
        character_table_chars = list(self.character_table.irreps.values())
        self.sym_sort = [[] for _ in character_table_syms]

        for nic_idx, nic_char in enumerate(nic_chars):
            for character_table_idx, character_table_char in enumerate(character_table_chars):
                if nic_char == character_table_char:
                    self.sym_sort[character_table_idx].append(nic_idx)
                    break
            else:
                raise ValueError(f"NIC {nic_idx} with characters {nic_char} does not match any irrep in the character table.")
        
        # print code that is ready to be copied
        print("=== sym_sort array generated ===")
        print("        self.sym_sort = np.array([")
        for row in self.sym_sort:
            print(f"            {row},")
        print("        ], dtype=object)\n")

        return
        

np.set_printoptions(linewidth=np.inf)   # don't wrap
np.set_printoptions(threshold=sys.maxsize)
print(Irrep("4.12","c2v").assign_irrep())
# print(Irrep("1.103","c2v").assign_irrep())   # testing on G2