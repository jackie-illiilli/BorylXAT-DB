"""Module specifically for reading Gaussian output files"""

from copy import deepcopy
import numpy as np
import os, glob
import re
import shutil
from rdkit import Chem

from . import Tool, FormatConverter

class bond_addition_function():

    def __init__(self):
        self.addition = {"atoms":{}, "value":{}, "type":{}}
        self.addition_idx = 0

    def compare(self, value1, value2, type):
        if type == "more":
            return value1 > value2
        elif type == "less":
            return value1 < value2
        elif type == "equal":
            return abs(value1 - value2) < 1E-6

    def add_function(self, atom_title_ids, value, type_):
        self.addition["atoms"][self.addition_idx] = atom_title_ids
        self.addition["type"][self.addition_idx] = type_
        self.addition["value"][self.addition_idx] = value
        assert type_ in ["more", "less", "equal"]

        self.addition_idx += 1
    
    def apply_addition(self, atom_idx_list, position):
        for addition_idx in range(self.addition_idx):
            atom_ids = [atom_idx_list[each] for each in self.addition["atoms"][addition_idx]]
            conf = Chem.rdchem.Conformer(len(atom_ids))
            for i, atom_id in enumerate(atom_ids):
                conf.SetAtomPosition(i, position[atom_id][:3])
            if len(atom_ids) == 2:
                bond_length = Chem.rdMolTransforms.GetBondLength(conf, 0, 1)
                if not self.compare(bond_length, self.addition["value"][addition_idx], self.addition["type"][addition_idx]):
                    return False
            elif len(atom_ids) == 3:
                bond_angle = Chem.rdMolTransforms.GetAngleDeg(conf, 0, 1, 2)
                if not self.compare(bond_angle, self.addition["value"][addition_idx], self.addition["type"][addition_idx]):
                    return False
            elif len(atom_ids) == 4:
                dihedral_angle = Chem.rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3)
                if not self.compare(dihedral_angle, self.addition["value"][addition_idx], self.addition["type"][addition_idx]):
                    return False
        return True

            


class Logfile():
    
    def __init__(self, file_dir, mol_file_dir=None,read_title=True, freq_warning=False, bond_attach_std='mol', bond_addition_function=None, bond_ignore_list=None, ignore_print=False):
        """Read Gaussian output file

        Args:
            file_dir (_type_): .log output file path
            mol_file_dir (_type_, optional): Corresponding .mol file path, some functions depend on the mol file. Defaults to None.
            read_title (bool, optional): Whether to read title information according to the specified template. Defaults to True.
            freq_warning (bool, optional): Ignore (imaginary frequency and spin multiplicity) warnings. Defaults to False.
            ignore_print (bool, optional): Whether to suppress warning print output from this Logfile. Defaults to False.

        Returns:
            _type_: _description_
        """        
        self.ignore_print = ignore_print
        self.file_dir = file_dir
        self.mol_file_dir = mol_file_dir
        with open(file_dir, "rt") as rf:
            filelines = rf.readlines()
            filelines = [line for line in filelines if line != ""]
        self.filelines = filelines
        self.normal_end = self.is_normal_end()
        if not self.normal_end:
            self.find_error_reason()
        self.S_2 = self.read_S_2()
        self.S_2_after_annihilation = self.read_S_2(after_annihilation=True)
        if read_title:
            self.title = self.read_title()
        self.charge, self.multiplicity = self.read_charge_multiplicity()
        # Special for Charge != 0
        if self.charge != 0:
            self.normal_end = False
            self.charge = 0
            self.error_reason = "link 9999"

        if self.multiplicity == -1:
            self._print("It's not a logfile")
            return None

        self.method = self.read_method().lower()
        self.file_type = 'SPE'
        if "irc" in self.method:
            self.file_type = 'IRC'
        elif "modredundant" in self.method:
            self.file_type = "OM"
        elif "readfc" in self.method or "calcfc" in self.method:
            self.file_type = "TS"
        elif "opt" in self.method:
            self.file_type = "OPT"

        self.symbol_list, self.first_atom_position = self.read_first_position()
        if self.symbol_list == None:
            self._print("It's a wrong file with unknown wrong")
            return None
        
        if self.file_type != "SPE":
            self.running_positions = self.read_running_position()
            if self.running_positions is not None and len(self.running_positions) != 0:
                self.running_rmsd = self.react_RMSD()
            if self.normal_end and "freq" in self.method:
                self.unreal_freq, self.unreal_freq_matrix, self.first_unreal_freq = self.read_unreal_freq(freq_warning)
            else:
                self.unreal_freq = -1
                self.unreal_freq_matrix = []
        else: 
            self.running_positions, self.unreal_freq, self.unreal_freq_matrix = [], 0, []
        
        if self.file_type == "OM":
            self.freeze, self.difreeze = self.read_freeze()
        else:
            self.freeze, self.difreeze = [], []


        self.all_engs, self.opt_engs = self.read_log_eng()
        if self.normal_end:
            self.running_time = self.read_log_time()
        else:
            self.engs = []
            self.running_time = 0
        if self.mol_file_dir != None and self.file_type != "SPE" and len(self.running_positions) != 0 and read_title:
            self.bond_attach = self.check_bond_attach(standard_file=bond_attach_std, bond_addition_function=bond_addition_function, bond_ignore_list=bond_ignore_list)
        else:
            self.bond_attach = True

        if self.file_type == "IRC" and read_title:
            self.irc_result = self.irc_check()

        if self.file_type == "TS" and read_title:
            if self.unreal_freq != 1:
                self._print("%s is not a right TS for unreal freq num of %d" % (self.file_dir, self.unreal_freq) )
                self.is_right_ts = False
            else:
                self.is_right_ts = self.check_om(False)
                pass

    def _print(self, *args, **kwargs):
        if not self.ignore_print:
            print(*args, **kwargs)

    def is_normal_end(self):
        """Detect "Normal termination of Gaussian"

        Args:
            file_dir (str): filedir

        Returns:
            bool: whether it end normal
        """    
        lastline = self.filelines[-1]
        if lastline.find(" Normal termination of Gaussian") == -1:
            self._print("%s didn't run successful" % self.file_dir)
            return False
        else: return True

    def read_title(self):
        """Read title information according to the $$$$Title####charge???? template

        Returns:
            list: [title, int of charge]
        """        
        allfile = "".join(self.filelines)
        title = ""
        charge = 0
        # title = $$$$Title####charge????
        if allfile.find("$$$$") == -1 or allfile.find("####") == -1 or allfile.find("????") == -1:
            title = ""
            self._print("Can't find title")
        else:
            title = allfile[allfile.find("$$$$") + 4: allfile.find("####")]
            charge = allfile[allfile.find("####") + 4: allfile.find("????")]
            title = title.split()
            if len(title) > 1:
                title = [int(each) for each in title]
        return title

    def read_charge_multiplicity(self):
        """Read charge and spin multiplicity

        Returns:
            list: [charge, multiplicity]
        """        
        line_id, line = Tool.find_first_line(self.filelines, 'Charge = ', 'in')
        if line_id is None:
            self._print("Can't find charge and multiplicity")
            return 0, -1
        charge = int(line.split()[2])
        multiplicity = int(line.split()[-1])
        return charge, multiplicity

    def read_S_2(self, after_annihilation=False):
        """Read the final S**2 value from Gaussian output.

        Args:
            after_annihilation (bool, optional): When Gaussian prints an
                annihilation-corrected value, return the final "after" value.
                Defaults to False.

        Returns:
            float: Final S**2 value, or -1 if it cannot be found.
        """
        s2_value = -1
        s2_after_annihilation = -1
        float_pattern = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
        s2_pattern = re.compile(r"<S\*\*2>\s*=\s*(%s)" % float_pattern)
        annihilation_pattern = re.compile(
            r"S\*\*2\s+before\s+annihilation\s+(%s),?\s+after\s+(%s)" %
            (float_pattern, float_pattern)
        )

        for line in self.filelines:
            s2_match = s2_pattern.search(line)
            if s2_match:
                s2_value = float(s2_match.group(1))
                s2_after_annihilation = -1

            annihilation_match = annihilation_pattern.search(line)
            if annihilation_match:
                s2_value = float(annihilation_match.group(1))
                s2_after_annihilation = float(annihilation_match.group(2))

        if after_annihilation and s2_after_annihilation != -1:
            return s2_after_annihilation
        return s2_value

    def read_first_position(self):
        """Read the element names and input coordinates of the molecule from the output file; can run even if Gaussian reports certain errors

        Returns:
            list: [atom list, position(n*3)]
        """        
        start_id = Tool.find_first_line(self.filelines, "Symbolic Z-matrix:", 'in')[0]
        if start_id != None:
            start_index = start_id + 2
            end_index = Tool.find_first_line(self.filelines[start_index:], ' \n', 'all')[0] + start_index
            if end_index is None:
                self._print("%s didn't have structure" % self.file_dir)
                return None, None
            # start_index = [i for i, line in enumerate(self.filelines) if line.find("Symbolic Z-matrix:") >= 0][0] + 2
            # end_index = [i + start_index for i, line in enumerate(self.filelines[start_index:]) if line == ' \n'][0]
            orientation = [line.split() for line in self.filelines[start_index: end_index]]
            symbol_list, position = [], []
            for each in orientation:
                if len(each)!= 4: return None, None
                symbol_list.append(each[0])
                position.append([float(each[1]), float(each[2]), float(each[3])]) 
            assert len(symbol_list) == len(position) 
            return symbol_list, np.array(position) 
        else:
            try:
                symbol, positions = self.read_running_position(read_first=1)
                return symbol, positions
            except:
                self._print("%s didn't have structure" % self.file_dir)
                return None, None

    def read_running_position(self, read_first=False):
        """Read coordinates for each step; the last coordinate is the stationary point coordinate

        Args:
            read_first (bool, optional): Whether only the first coordinate is needed. Defaults to False.

        Returns:
            np.array: position of atoms(n*3)
        """        
        orientation_sign = "Input orientation:"
        # if self.file_type == 'IRC' or read_first:
        #     orientation_sign = "Input orientation:"
        # else:
        #     orientation_sign = "Standard orientation:"
        start_indexs = [i for i, line in enumerate(self.filelines) if line.find(orientation_sign) >= 0]
        if len(start_indexs) == 0:
            orientation_sign = "Standard orientation:"
            start_indexs = [i for i, line in enumerate(self.filelines) if line.find(orientation_sign) >= 0]
        if len(start_indexs) == 0:
            self._print("%s Even not Input Structure, wrong file maybe" % self.file_dir)
            return None
        all_positions = []
        for each_start_index in start_indexs:
            start_index = each_start_index + 5
            end_index = Tool.find_first_line(self.filelines[start_index:], ' ---', 'in')[0] + start_index
            if end_index is None:
                break
            orientation = [line.strip("\n").split() for line in self.filelines[start_index: end_index]]
            position = []
            symbol = []
            for each in orientation:
                if len(each) != 6: return all_positions
                position.append([float(each[3]), float(each[4]), float(each[5])])
                symbol.append(int(each[1]))
            all_positions.append(position)
            
        if read_first:
            return symbol, all_positions[0]
        else:
            if len(all_positions) == 1:
                return np.array(all_positions)
            return np.array(all_positions)[1:]

    def read_method(self):
        method_id, method = Tool.find_first_line(self.filelines, ' #p', 'start')
        if method is None:
            method_id, method = Tool.find_first_line(self.filelines, ' #', 'start')
            if method is None:
                self._print("%s with not method" % self.file_dir)
                return None
        method_final_line_id, _ = Tool.find_first_line(self.filelines[method_id:], " -------", 'start')
        method = "".join([each.strip("\n").lstrip(" ") for each in self.filelines[method_id: method_id + method_final_line_id]])
        method = method.split('#p ')[-1]
        return method

    def read_unreal_freq(self, freq_warning=True):
        """Detect imaginary frequencies, output the number of imaginary frequencies, the vibrational coordinates of the largest imaginary frequency, and its frequency value

        Args:
            file_dir (_type_): _description_

        Returns:
            int, np.array, float: num of unreal freqences
        """    
        start_index = [i for i, line in enumerate(self.filelines) if '(negative Signs)' in line]
        if len(start_index) != 0: start_index = start_index[-1]
        else: start_index = 0
        smallest_freq_index, smallest_freq_line = Tool.find_first_line(self.filelines[start_index:], ' Frequencies --', "start")
        smallest_freq_index += start_index
        if smallest_freq_index == None:
            self._print("%s didn't calc freq" % self.file_dir)
            return -1, []
        smallest_freq_list = smallest_freq_line.strip("\n").split()[2:]
        num_unreal_freq = sum([1 for each in smallest_freq_list if float(each) < 0])
        # freq_fileline = [[lid, Tool.remove_space(line.strip("\n"))[2]] for lid, line in enumerate(self.filelines) if line.startswith(" Frequencies --")][0]
        matrix = []
        start_id = smallest_freq_index + 5
        while(1):
            line = self.filelines[start_id].strip("\n").split()
            if len(line) < 5:
                break
            line = [float(each) for each in line[2:5]]
            matrix.append(line)
            start_id += 1
        if freq_warning:
            self._print("%s have unreal freq" % self.file_dir)
        return num_unreal_freq, matrix, smallest_freq_list[0]

    def read_log_eng(self): 
        """Read energy information from the optimization file: electronic energy, zero-point energy, thermal energy correction, enthalpy correction, and Gibbs free energy correction

        Args:
            gjffile (str): *.log

        Returns:
            ee, zpc, cor_Energy, cor_Enthalpies, cor_Gibbs : [list with electronic energy, zero-point energy, thermal energy correction, enthalpy correction, and Gibbs free energy correction]
        """
        all_engs = []
        opt_engs = []
        ee_line = [[idx, each] for idx, each in enumerate(self.filelines) if " SCF Done: " in each] 
        if len(ee_line) == 0:
            self._print("%s, can't find any engs" % self.file_dir)
            return all_engs, opt_engs
        for line_id, each_ee_line in ee_line:
            # ee_line = [line for i, line in enumerate(
            #     self.filelines[start_index:]) if line.startswith(" SCF Done: ")][-1]
            if each_ee_line == None: 
                ee = -1
            else:
                ee = float(each_ee_line.strip("\n").split()[4])
            opt_engs.append(ee)
        try:
            start_idx = ee_line[-1][0]
            zpc_line = Tool.find_first_line(self.filelines[start_idx:], " Zero-point correction=", "in")[-1]
            zpc = zpc_line.strip("\n").split(" ")[-2]
            cor_ee = Tool.find_first_line(self.filelines[start_idx:], " Thermal correction to Energy=", "in")[-1].strip("\n").split(" ")[-1]
            cor_Enthalpies = Tool.find_first_line(self.filelines[start_idx:], " Thermal correction to Enthalpy=", "in")[-1].strip("\n").split(" ")[-1]
            cor_Gibbs = Tool.find_first_line(self.filelines[start_idx:], " Thermal correction to Gibbs Free Energy=", "in")[-1].strip("\n").split(" ")[-1]
            all_engs = [ee, zpc, cor_ee, cor_Enthalpies, cor_Gibbs]
            all_engs = [float(each) for each in all_engs]
        except:
            return [ee], opt_engs
        return all_engs, opt_engs

    def read_log_time(self):
        """Read the elapsed time in minutes for all calculation processes

        Returns:
            float: 
        """    
        alltime = 0
        for each in self.filelines:
            if each.startswith(" Elapsed time:"):
                times = list(each.strip("\n").split())
                alltime += float(times[2]) * 24 * 60 + float(times[4]) * 60 + float(times[6]) + float(times[8]) / 60
        return alltime

    def read_freeze(self):
        """Read the frozen bond lengths and dihedral angles in the input file, and their corresponding atom numbers

        Returns:
            list: [[[a, b], [c, d]], [[a,b,c,d], [e,f,g,h]]]
        """        
        freeze_start_line_id, _ = Tool.find_first_line(self.filelines, "The following ModRedundant", "in")
        freeze, difreeze = [], []
        for line in self.filelines[freeze_start_line_id + 1:]:
            line = line.strip('\n').split()
            if len(line) < 4:
                break
            if line[0] == "B":
                freeze.append([int(line[1]), int(line[2])])
            elif line[0] == "D":
                difreeze.append([int(line[1]), int(line[2]), int(line[3]), int(line[4])])
        return freeze, difreeze

    def find_error_reason(self):
        """Find the reason for the Gaussian error
        """        
        errorline_id = [i for i, line in enumerate(self.filelines) if " Error termination" in line]
        # errorline_id = Tool.find_first_line(self.filelines,"start")[0]
        if len(errorline_id) == 0 and not self.normal_end:
            self._print(self.file_dir, "Should have not finished")
            error_reason_line = "unfinished"
        else:
            errorline_id = errorline_id[-1]
            error_reason_line = self.filelines[errorline_id - 1]
            error_reason_line = error_reason_line.strip('.\n')
        self.error_reason = error_reason_line

    def modify_method(self, old_method, addition_part='opt', new_word="maxcycles=100"):
        method = old_method.split()
        opt_part = [[idx, each] for idx, each in enumerate(method) if addition_part in each]
        if len(opt_part) == 0:
            method = old_method
        else:
            opt_id = opt_part[0][0]
            opt_part = opt_part[0][1]
            if new_word.split("=")[0] in opt_part:
                method = old_method
            else:
                if "=" not in opt_part:
                    opt_part += f"=({new_word})"
                elif "(" in opt_part:
                    opt_part = opt_part.replace(")", f",{new_word})")
                else:
                    opt_part = opt_part.replace("=", "=(")
                    opt_part += f",{new_word})"
                method[opt_id] = opt_part
                method = " ".join(method)
        return method

    def solve_error_logfile(self, new_log_dir, Inv_dir='inv', yqc_dir='yqc', savechk=None, readchk=None, maxcycles=0, method=None):
        """Select the lowest energy structure and continue running with the same method; suitable for Link 9999, etc.

        Args:
            new_log_dir (str): Address for storing error output files and new input files
            move_file (bool, optional): Whether to move the error file. Defaults to True.
            savechk (_type_, optional): Whether to write savechk information in the input file. Defaults to None.
            readchk (_type_, optional): Whether to write readchk information in the input file. Defaults to None.

        """        
        reason = self.error_reason
        if not os.path.isdir(new_log_dir):
            os.mkdir(new_log_dir)
        new_gjf_name = new_log_dir + "/" + os.path.split(self.file_dir)[-1].split(".")[0] + '.gjf'
        if self.file_type in ["IRC"]:
            return 0
        if "FormBX" in reason or "Linear angle" in reason or "Tors failed for dihedral" in reason:
            new_position = self.solve_l103_problem()
            title = " ".join(str(each) for each in self.title)
            FormatConverter.block_to_gjf(self.symbol_list, new_position, new_gjf_name, self.charge, self.multiplicity, title, self.method, freeze=self.freeze, difreeze=self.difreeze, savechk=savechk, readchk=readchk)
            new_log_name = new_log_dir + "/" + os.path.split(self.file_dir)[-1]
            shutil.move(self.file_dir, new_log_name)
        elif "Convergence failure" in reason:
            title = " ".join(str(each) for each in self.title)
            new_log_name = yqc_dir + "/" + os.path.split(self.file_dir)[-1]
            new_gjf_name = yqc_dir + "/" + os.path.split(self.file_dir)[-1].split(".")[0] + '.gjf'
            if not os.path.isdir(yqc_dir):
                os.mkdir(yqc_dir) 
            method = self.method + " scf=yqc"
            new_log_name = yqc_dir + "/" + os.path.split(self.file_dir)[-1]
            shutil.move(self.file_dir, new_log_name)
            FormatConverter.block_to_gjf(self.symbol_list, self.running_positions[-1], new_gjf_name, self.charge, self.multiplicity, title, method, freeze=self.freeze, difreeze=self.difreeze, savechk=savechk, readchk=readchk)
            
        elif "Inv3" in reason:
            title = " ".join(str(each) for each in self.title)
            new_gjf_name = Inv_dir + "/" + os.path.split(self.file_dir)[-1].split(".")[0] + '.gjf'
            new_log_name = Inv_dir + "/" + os.path.split(self.file_dir)[-1] 
            if not os.path.isdir(Inv_dir):
                os.mkdir(Inv_dir) 
            shutil.move(self.file_dir, new_log_name)

        else:
            self._print("%s. Error Reason is link 9999 or unfinished." % self.file_dir)
            title = " ".join(str(each) for each in self.title)
            if len(self.running_positions) != 0:
                opt_engs = deepcopy(self.opt_engs) 
                idx = len(self.running_positions) - 1
                for _ in range(len(opt_engs)):
                    # try:
                    #     min_index =np.argmin(opt_engs)
                    #     if min_index < len(self.running_positions) - 1:
                    #         opt_engs[min_index] = 0
                    #         continue
                    #     if self.mol_file_dir != None:
                    #         bond_result = self.check_bond_attach('log', conf_id=min_index)
                    #         if bond_result:
                    #             target_position = self.running_positions[min_index]
                    #             break
                    #     else:
                    #         target_position = self.running_positions[min_index]
                    #         break
                    #     opt_engs[min_index] = 0
                    # except:
                    #     opt_engs[min_index] = 0
                    try:            
                        target_position = self.running_positions[idx]
                        break
                    except:
                        idx -= 1
                        target_position = self.first_atom_position
            else:
                target_position = self.first_atom_position
            if method == None:
                method = self.method
                if maxcycles != 0:
                    method = self.modify_method(method, new_word=f'maxcycles={maxcycles}')
            new_log_name = new_log_dir + "/" + os.path.split(self.file_dir)[-1]
            shutil.move(self.file_dir, new_log_name)
            if self.file_type == "TS":
                FormatConverter.block_to_gjf(self.symbol_list, target_position, new_gjf_name, self.charge, self.multiplicity, title, method, freeze=self.freeze, difreeze=self.difreeze, savechk=savechk, readchk=readchk)
            else:
                FormatConverter.block_to_gjf(self.symbol_list, target_position, new_gjf_name, self.charge, self.multiplicity, title, method, freeze=self.freeze, difreeze=self.difreeze, savechk=savechk, readchk=readchk)

    def l103_error_idx(self):
        """Find the atom numbers of the bond angles and dihedral angles reported in the L103 error

        Returns:
            list: [[a,b,c], [a,b,c,d]]
        """        
        angle_idx = np.zeros(3)
        dihedral_idx = np.zeros((0, 4))
        for line in self.filelines[-15:-5]:
            line = line.strip('\n').strip(' ')
            if 'Bend failed for angle' in line:
                ww = line.split()
                angle_idx = np.array([ww[4], ww[6], ww[8]], dtype=int)

            elif 'Tors failed for dihedral' in line:
                ww = line.split()
                tmp_list = np.array([[ww[4], ww[6], ww[8], ww[10]]], dtype=int)
                dihedral_idx = np.append(dihedral_idx, tmp_list, axis=0)
            
            elif 'Linear angle in Tors.' in line:
                dihedral_idx = np.zeros((1, 4))
        if not dihedral_idx.shape[0]:
            dihedral_idx = np.zeros((1, 4))
        dihedral_idx = dihedral_idx.astype('int',copy=False)
        angle_idx= angle_idx.astype('int',copy=False)
        if not angle_idx.all() and dihedral_idx.all():
            tmp_idx = dihedral_idx[0]
            if dihedral_idx.shape[0] == 1:
                angle_idx = tmp_idx[1:]
            else:
                if (dihedral_idx[:, :3] == tmp_idx[:3]).all():
                    angle_idx = tmp_idx[:3]
                else:
                    angle_idx = tmp_idx[1:]
        return angle_idx, dihedral_idx
                    
    def l103_adjust(self):
        """Appropriately adjust bond angles and dihedral angles
        """        
        def get_Rotation_M(axial_v, theta):
            v = np.array(axial_v[:3])
            # normalization
            u, v, w = v/np.linalg.norm(v)
            a = theta
            R_M = np.array([[u**2+(1-u**2)*np.cos(a),       u*v*(1-np.cos(a))-w*np.sin(a),  u*w*(1-np.cos(a))+v*np.sin(a),  0],
                            [u*v*(1-np.cos(a))+w*np.sin(a), v**2+(1-v**2) *
                            np.cos(a),        v*w*(1-np.cos(a))-u*np.sin(a),  0],
                            [u*w*(1-np.cos(a))-v*np.sin(a), v*w*(1-np.cos(a)) +
                            u*np.sin(a),  w**2+(1-w**2)*np.cos(a),        0],
                            [0,                             0,                              0,                              1]])
            return R_M

        if (self.dihedral_idx[:, :3] == self.angle_idx).all():
            atom_to_be_adjusted_idx = self.angle_idx[0]-1
            o_idx = self.angle_idx[1]-1
            v_idx = self.angle_idx[2]-1
        elif not self.dihedral_idx.all() or (self.dihedral_idx[:, 1:] == self.angle_idx).all():
            atom_to_be_adjusted_idx = self.angle_idx[-1]-1
            o_idx = self.angle_idx[-2]-1
            v_idx = self.angle_idx[-3]-1
        elif ((self.dihedral_idx[:, :3] == self.angle_idx)+(self.dihedral_idx[:, 1:] == self.angle_idx)).all():
            atom_to_be_adjusted_idx = self.angle_idx[-1]-1
            o_idx = self.angle_idx[-2]-1
            v_idx = self.angle_idx[-3]-1
        else:
            atom_to_be_adjusted_idx = self.angle_idx[-1]-1
            o_idx = self.angle_idx[-2]-1
            v_idx = self.angle_idx[-3]-1
        Coord = deepcopy(self.running_positions[-1])
        symbol_list = self.symbol_list
        angle_v = Coord[v_idx] - Coord[o_idx]
        changing_v = Coord[atom_to_be_adjusted_idx] - Coord[o_idx]
        axial_v = np.cross(angle_v, changing_v)
        axial_v = axial_v/np.linalg.norm(axial_v)
        if symbol_list[atom_to_be_adjusted_idx] == 'H':
            theta = np.pi/6.
        else:
            theta = np.pi/36.
    
        R_M = get_Rotation_M(axial_v, theta)
        o_coord = Coord[o_idx]
        tmp_v = np.append(changing_v, [1])
        tmp_v = np.dot(tmp_v, R_M)
        tmp_v = np.around(np.delete(tmp_v, 3), decimals=8)
        Coord[atom_to_be_adjusted_idx] = tmp_v + o_coord
        return Coord

    def solve_l103_problem(self):        
        self.angle_idx, self.dihedral_idx = self.l103_error_idx()
        if self.angle_idx.all():
            new_position = self.l103_adjust()
        else:
            self._print('!!! Warning file: %s; Unknown Error：%s'%(self.file_dir, self.error_reason))
            new_position = self.running_positions[-2]
        return new_position

    def check_bond_attach(self, standard_file='mol', print_num = False, conf_id=-1, bond_addition_function=None, bond_ignore_list=None):
        """With the help of the Mol file, check if atom connections are between 0.75 and 1.3 times the original bond length.

        Args:
            standard_file (str, optional): Mol file. Defaults to 'mol'.
            print_num (bool, optional): . Defaults to False.
            conf_id (int, optional): Conformation ID of Mol. Defaults to -1.

        Returns:
            _type_: _description_
        """        
        try:
            # assert standard_file == "mol"
            if self.running_positions is None:
                return False
            mol = Chem.MolFromMolFile(self.mol_file_dir, removeHs=False, sanitize=False)
            for atom_id, atom in enumerate(mol.GetAtoms()):
                if atom.GetSymbol() != self.symbol_list[atom_id] and atom.GetAtomicNum() != self.symbol_list[atom_id]:
                    self._print("wrong with", atom.GetIdx(), atom.GetAtomicNum(), atom_id, self.symbol_list[atom_id])
                    return False
            if standard_file == "mol":
                position = mol.GetConformer(0).GetPositions()
            else:
                position = self.first_atom_position
            new_position = self.running_positions[conf_id]
            except_idxs = []
            if self.file_type == 'TS' or self.file_type == "OM":
                title = [each - 1 for each in self.title]
                if bond_ignore_list == None:
                    except_idxs = []
                else:
                    except_idxs = [[title[idx] for idx in each] for each in bond_ignore_list]
            for bond in mol.GetBonds():
                ignore=False
                start_atom_id = bond.GetBeginAtomIdx()
                end_atom_id = bond.GetEndAtomIdx()
                distance_a = Tool.get_atoms_distance(position[start_atom_id], position[end_atom_id])
                distance_b = Tool.get_atoms_distance(new_position[start_atom_id], new_position[end_atom_id])
                num = distance_a / distance_b
                for except_idx in except_idxs:
                    if start_atom_id in except_idx and end_atom_id in except_idx:
                        ignore=True
                if print_num:
                    self._print("%d %d %.5f" % (start_atom_id, end_atom_id, num))
                if (num <= 0.75 or num >= 1.3) and not ignore:
                    self._print(os.path.split(self.file_dir)[-1], start_atom_id, end_atom_id, "with a wrong distance", num)
                    return False
            if bond_addition_function is not None:
                if not bond_addition_function.apply_addition(title, new_position):
                    return False
            return True
        except:
            return False
    
    def irc_check(self):
        title = [each - 1 for each in self.title]
        atom1, atom2, atom3 = title[:3]
        new_position = self.running_positions[-1]
        except_idxs = [[atom1, atom2], [atom2, atom3]]
        for id, except_idx in enumerate(except_idxs):
            start_atom_id = except_idx[0]
            end_atom_id = except_idx[1]
            distance = Tool.get_atoms_distance(new_position[start_atom_id], new_position[end_atom_id])
            if distance <= 2.2:
                self._print(os.path.split(self.file_dir)[-1], start_atom_id, end_atom_id, "atom may ircorrect", distance)
                return id
        return -1

    def unreal_freq_improve(self, new_log_dir, savechk=None, readchk=None, method=None):
        """Solve imaginary frequency problem: carry over 1.1 times the imaginary frequency vibration into the next optimization step

        Args:
            logfile (_type_): _description_
            newdir (str, optional): _description_. Defaults to 'err_imp'.
        """    
        position = self.running_positions[-1]
        title = " ".join(str(each) for each in self.title)
        new_position = np.array(self.unreal_freq_matrix) * 1.1 + np.array(position)
        if not os.path.isdir(new_log_dir):
            os.mkdir(new_log_dir)
        new_log_name = new_log_dir + "/" + os.path.split(self.file_dir)[-1]
        new_gjf_name = new_log_dir + "/" + os.path.split(self.file_dir)[-1].split(".")[0] + ".gjf"
        shutil.move(self.file_dir, new_log_name)
        if method == None:
            method = self.method
        FormatConverter.block_to_gjf(self.symbol_list, new_position, new_gjf_name, self.charge, self.multiplicity, title, method, freeze=self.freeze, difreeze=self.difreeze , savechk=savechk, readchk=readchk)

    def react_RMSD(self):
        """Determine the RMSD change during the optimization process

        Args:
            log_dir (_type_): _description_

        Returns:
            float: RMSD
        """    
        position_start = self.running_positions[0]
        position_end = self.running_positions[-1]
        sum_delta2 = 0
        for each_start, each_end in zip(position_start, position_end):
            delta = each_start - each_end
            sum_delta2 += delta @ delta
        sum_delta2 /= len(position_start)
        return np.sqrt(sum_delta2)

    def check_om(self, return_value=False, set_num=0.1, strict=True):
        """Read the transition state structure and determine if the imaginary frequency corresponds to the reaction site

        Args:
            file_dir (str): _description_
            assert_title (_type_, optional): whether define a title. Defaults to None.

        Returns:
            Bool: 
        """    
        title = [each - 1 for each in self.title]
        position = deepcopy(self.running_positions[-1])
        new_position = deepcopy(self.unreal_freq_matrix)
        assert self.unreal_freq != 0
        # return new_position, position
        new_position += position
        bond1_dist_change = Tool.get_atoms_distance(new_position[title[0]], new_position[title[1]]) - Tool.get_atoms_distance(position[title[0]], position[title[1]])
        bond2_dist_change = Tool.get_atoms_distance(new_position[title[1]], new_position[title[2]]) - Tool.get_atoms_distance(position[title[1]], position[title[2]])
        if return_value:
            return (bond1_dist_change, bond2_dist_change)
        if strict:
            TF = (bond1_dist_change > set_num and bond2_dist_change < set_num) or (bond1_dist_change < set_num and bond2_dist_change > set_num)
        else:
            TF = np.abs(bond1_dist_change) > set_num or np.abs(bond2_dist_change) > set_num
        if TF:
            return 1
        else:
            self._print("%s may find wrong TS" % self.file_dir)
        return 0

    def read_charge_spin_density(self):
        lines = self.filelines
        start_ids = [id for id, line in enumerate(lines) if line.startswith(" Mulliken charges:")]
        if len(start_ids) == 0:
            start_ids = [id for id, line in enumerate(lines) if line.startswith(" Mulliken charges and spin densities:")]
        start_id = start_ids[-1] + 2
        end_id = [id for id, line in enumerate(lines[start_id:]) if line.startswith(" Sum of Mulliken charges =")][-1] + start_id
        charges = []
        spin_density = []
        for line in lines[start_id:end_id]:
            line_split = line.strip("\n").split()
            if len(line_split) == 3:
                charges.append(line_split[-1])
            else:
                charges.append(line_split[-2])
                spin_density.append(line_split[-1])
        charges = [float(each) for each in charges]
        spin_density = [float(each) for each in spin_density]
        return charges, spin_density
    
    def orbit_str_split(self, str_):
        result_lists = []
        str_lists = str_.split()
        for each in str_lists:
            if each.count("-") >= 2:
                tmp_lists = [i for i in each.split("-") if i != ""]
                result_lists += ["-" + a for a in tmp_lists]
            else:
                result_lists.append(each)
        return result_lists

    def read_orbit_eng(self, HOMO_index = [-2, -1], LUMO_index=[0,1]):
        """Read the energies of HOMO and LUMO orbitals

        Args:
            HOMO_index (list, optional): Index corresponding to the HOMO orbital. Defaults to [-2, -1].
            LUMO_index (list, optional): Index corresponding to the LUMO orbital. Defaults to [0,1].

        Returns:
            list: 
        """        
        occ_eng = []
        virt_eng = []
        lines = self.filelines
        # start_id = [id for id, line in enumerate(lines) if line.startswith(" The electronic state is")][-1] + 1
        # end_id = [id for id, line in enumerate(lines[start_id:]) if line.startswith("          Condensed to atoms (all electrons):")][-1] + start_id
        occ_orbits = [line for line in lines if line.startswith(" Alpha  occ. eigenvalues")]
        virt_orbits = [line for line in lines if line.startswith(" Alpha virt. eigenvalues")]
        for occ_orbit in occ_orbits:
            occ_eng += [float(each) * 627.5 for each in self.orbit_str_split(occ_orbit.strip("\n").split("--")[-1])]
        for virt_orbit in virt_orbits:
            virt_eng += [float(each) * 627.5 for each in self.orbit_str_split(virt_orbit.strip("\n").split("--")[-1])]
        occ_eng = [occ_eng[each] for each in HOMO_index]
        virt_eng = [virt_eng[each] for each in LUMO_index]

        occ_eng_beta = []
        virt_eng_beta = []
        occ_orbits = [line for line in lines if line.startswith("  Beta  occ. eigenvalues")]
        virt_orbits = [line for line in lines if line.startswith("  Beta virt. eigenvalues")]
        if len(occ_orbits) == 0:
            return occ_eng + virt_eng
        for occ_orbit in occ_orbits:
            occ_eng_beta += [float(each) * 627.5 for each in self.orbit_str_split(occ_orbit.strip("\n").split("--")[-1])]
        for virt_orbit in virt_orbits:
            virt_eng_beta += [float(each) * 627.5 for each in self.orbit_str_split(virt_orbit.strip("\n").split("--")[-1])]
        occ_eng_beta = [occ_eng_beta[each] for each in HOMO_index]
        virt_eng_beta = [virt_eng_beta[each] for each in LUMO_index]
        return occ_eng + virt_eng + occ_eng_beta + virt_eng_beta
    
    def get_dipole(self):
        """Read the dipole moment of the molecule

        Returns:
            _type_: _description_
        """        
        lines = self.filelines
        dipole_line_id = [line_id for line_id, line in enumerate(lines) if line.startswith(" Dipole moment ")]
        if len(dipole_line_id) == 0:
            self._print(self.filelines, "did't have a dipole moment")
            return -1

        dipole_line = lines[dipole_line_id[-1] + 1]
        dipole_moment = float(dipole_line.strip("\n").split()[-1])
        return dipole_moment

if __name__ == "__main__":
    pwd = os.getcwd()
    log_files = glob.glob(pwd + "/*.log")
    # improve_dir = pwd + "/improve"
    for log_file in log_files:
        opt_log = Logfile(log_file)
        print("log_file", opt_log.read_orbit_eng([-1], [0]))
        # if not opt_log.normal_end:
        #     opt_log.solve_error_logfile(improve_dir)
        # elif opt_log.unreal_freq:
        #     opt_log.unreal_freq_improve(improve_dir)
