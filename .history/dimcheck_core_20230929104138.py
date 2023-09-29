'''A module for dimension check.'''


import copy
import json
import re
import hashlib
import pickle
import os
import sympy as sp
from sympy.physics.units.definitions import kg, m, s, A, mol, cd, K
from sympy.physics.units.quantities import Quantity
from sympy import sympify
from queue import Queue
from typing import Any
from .dimcheck_class import DIMCHECK_MODULE_PATH


class MultiSymDefError(Exception):
    '''Exception raised when users choose the same symbol with different definitions.
    '''
    pass

class ExprBktError(Exception):
    '''Exception raised when the left and right braket in the expression does not match.
    '''
    pass

class MissingBktError(Exception):
    '''Exception raised when users input a quantity without a braket.
    '''
    pass

class MissingQuantError(Exception):
    '''Exception raised when users input a quantity without a value.
    '''
    pass

class InvalidSymFormatError(Exception):
    '''Exception raised when users input a quantity or alias with incorrect format.
    '''
    pass

class DerivingQuantError(Exception):
    '''Exception raised when the quantities cannot be totally derived.
    '''
    pass

class UnknowQuantFoundError(Exception):
    '''Exception raised when a quantity not defined is found.
    '''
    pass







from scipy.optimize import newton
from queue import Queue
from numpy import ndarray



def _split_dir_file(s):
    dir_end = s.rfind(os.sep)
    if dir_end != -1:
        return s[:dir_end], s[dir_end + 1 :]
    else:
        return None, s
    

def _divide_file_path(path):
    dir_path,filename = _split_dir_file(path)
    filename_parts=filename.split(".")
    if len(filename_parts)!=2:
        print("filename ({}) is not in XXX.xxx form.".format(filename))
    filename=filename_parts[0]
    postfix=filename_parts[-1]
    return {"filename":filename,"dir_path":dir_path,"postfix":postfix}



class Dimcheck:
    '''Create a Dimcheck object based on the config_file.
    '''
    def __init__(self, config_file="./si.json", env="si", is_save=False):
        '''Initialize the object with configuration file

        Parameters
        ----------
        config_file : str, optional
            Filepath of definitions of quantities, by default "./si.json"
        env : str, optional
            Environment of the unit system. So far only "si" is allowed, by default "si"
        is_save : bool, optional
            Whether to save the object in a binary form or not. A hash method is used to check the configuration file , by default False
        '''
        self.config_file = config_file
        self.env = env
        self._is_save = is_save
        self._pickle_file = ""

        serialized_dimchecker = self._compare_serialized_hash_match(config_file=config_file)
        if serialized_dimchecker:
            self._dimchecker = serialized_dimchecker
            self._is_save = False
        else:
            self._dimchecker = Dimchecker(config_file=config_file, env=env)

        if self._is_save:
            pickle_file = self._get_serialized_path(config_file)
            data = {}
            with open(pickle_file,"wb") as f:
                data["hash"] = self._cal_config_file_hash(config_file)
                data["env"] = env
                data["config_file"] = config_file
                data["dimcheck_obj"] = self._dimchecker
                pickle.dump(data,f)

    def dimension(self,s: str) -> Any:
        '''Get the dimension of the input

        Parameters
        ----------
        s : str
            Any quantity or the combination of the quantities

        Returns
        -------
        Any
            Depending on the input, it could be number, float, Quantity, etc.
        '''
        return self._dimchecker.dimension(s)
    
    def dim(self,s: str) -> Any:
        '''Alias of Dimcheck.dimension.
        
        Get the dimension of the input.

        Parameters
        ----------
        s : str
            Any quantity or the combination of the quantities

        Returns
        -------
        Any
            Depending on the input, it could be number, float, Quantity, etc.
        '''
        return self.dimension(s)

    def is_dc(self,lhs: str,rhs: str) -> bool:
        '''To see whether the quantities on the two sides are equal or not.

        If return False, a straightforward relation based on the configuration file will be checked. 

        Parameters
        ----------
        lhs : str
            Quantities on the left hand side
        rhs : str
            Quantities on the right hand side

        Returns
        -------
        bool
            Return True if the dimensions on two sides are equal, otherwise return False.
        '''
        return self._dimchecker.is_dc(lhs,rhs)

    def print_quant(self,s: str, is_include_alias=False) -> None:
        '''Print the possible quantity by deriving the combination of quantities.
        '''
        self._dimchecker.possible_quant(s, is_include_alias)
        
    def display_all_expr(self) -> None:
        '''Display all expressions based on the configuration file.
        '''
        self._dimchecker.display_all_expr()

    def display_all_quant(self) -> None:
        '''Display all quantities based on the configuration file.
        '''
        self._dimchecker.display_all_quant()

    def clean(self) -> None:
        '''Clean the binary file of current object.
        '''
        if os.path.exists(self._pickle_file):
            os.system("rm "+ self._pickle_file)


    def _compare_serialized_hash_match(self,config_file: str):
        serialized_dimcheck_obj = None
        pickle_file = self._get_serialized_path(config_file=config_file)
        self._pickle_file = pickle_file
        if os.path.exists(pickle_file):
            with open(pickle_file,"rb") as f:
                data = pickle.load(f)
                serialized_hash = data["hash"]
                config_file_hash = self._cal_config_file_hash(config_file)
                if serialized_hash == config_file_hash:
                    serialized_dimcheck_obj = data["dimcheck_obj"]
        
        return serialized_dimcheck_obj
    
    def _cal_config_file_hash(self, config_file, hash_method=hashlib.md5):
        digestobj=None
        with open(config_file,"r") as f:
            digestobj=hash_method()
            while True:
                data = f.read(4096)
                if not data:
                    break
                digestobj.update(data.encode("utf-8"))

        return digestobj.hexdigest()

    def _get_serialized_path(self, config_file: str):
        file_info=_divide_file_path(config_file)
        filename=file_info["filename"]
        dir_path=file_info["dir_path"]
        pickle_file=dir_path+os.sep+filename+".pickle"
        return pickle_file
    



class Dimchecker():
    ''' Initialize a dimensional checker with a certain configuration file.
    '''
    def __init__(self,config_file="./si.json",env="si"):
        
        self._env=env
        self._lhs_list=[]
        self._rhs_list=[]
        self._discription_list=[]
        self._quant_set=set()
        self._alias_map={}
        self._quant2alias_map={}
        self._config_file=config_file

        self.base_unit_dict={"kg":kg, "m":m, "s":s, "A":A, "cd":cd, "K":K, "mol":mol}
        

        ## Load the configuration file and check repeated definition of symbols
        with open(config_file,"r") as f:
            json_str=json.load(f)
            formula=json_str["quantities"]
            for x in formula:
                quant = x["quantity"]
                self._append_formulas(quant,x["rhs"],x["discription"])
                

                alias=x.get("alias")
                if alias:
                    if isinstance(alias,str):
                        self._append_alias(alias,quant)
                        self._quant2alias_map[quant] = [alias]
                    else:
                        for q in alias:
                            self._append_alias(q,quant)
                        self._quant2alias_map[quant] = alias
                else:
                    self._quant2alias_map[quant] = None
                
        

        print("Starting to parse the configuration file...")
        dependence_queue=Queue()
        unsolved_num=0
        self._drv_unit_dict=copy.deepcopy(self.base_unit_dict)
        for i in range(len(self._lhs_list)):
            self._check_sym_format(self._lhs_list[i])
            striped_str=self._strip_sqr_bkt(self._rhs_list[i],is_defining_symbol=True)
            expr=sp.simplify(sympify(striped_str).subs(self._drv_unit_dict))
            quants=re.findall("_dimcheck_quant_\w+",str(expr))
            is_unsolved=False
            tmp=[]
            for quant in quants:
                if not self._drv_unit_dict.get(quant):
                    tmp.append(quant)
            dependence_queue.put((expr,i,tmp))
            if is_unsolved==True:
                continue
            
            self._drv_unit_dict.update({self._strip_sqr_bkt(self._lhs_list[i],is_defining_symbol=True):expr})

        unsolved_num=dependence_queue.qsize()
        for i in range(unsolved_num):
            num=dependence_queue.qsize()
            for j in range(num):
                expr,idx,unsolved_quants=dependence_queue.get()
                is_solvable=True
                for quant in unsolved_quants:
                    if self._drv_unit_dict.get(quant) is None:
                        is_solvable=False
                if is_solvable:
                    expr=sp.simplify(expr.subs(self._drv_unit_dict))
                    self._drv_unit_dict.update({self._strip_sqr_bkt(self._lhs_list[idx],is_defining_symbol=True):expr})
                else:
                    dependence_queue.put((expr,idx,unsolved_quants))

        unsolved_num=dependence_queue.qsize()
        if unsolved_num!=0:
            for i in range(unsolved_num):
                expr,idx,unsolved_quants=dependence_queue.get()
                # Debug("{},{},{}".format(expr,idx,self._quant2str(unsolved_quants)))
            raise DerivingQuantError("The quantities {} can not be totally derived!".format(self._quant2str(unsolved_quants)))

        self._expr_dict={}
        for quant,expr in self._drv_unit_dict.items():
            if self._is_base_unit(quant):
                continue
            quants=self._expr_dict.get(expr)
            if not quants:
                self._expr_dict[expr]=[quant]
            else:
                self._expr_dict[expr].append(quant)

        print("Successfully parsed!")

        # self.display_all_quant(is_include_alias=False)
        # self.display_all_expr(is_include_alias=False)

    def dimension(self,s):
        striped_str=self._strip_sqr_bkt(s)
        expr=sp.simplify(sympify(striped_str).subs(self._drv_unit_dict))
        return expr


    def is_dc(self,lhs,rhs):
        l=self.dimension(lhs)
        r=self.dimension(rhs)
        if l==r:
            print("DC!")
            print("{} == {} ({})".format(lhs, rhs, l))
            return True
        else:
            print("There are some differences on two sides. {} != {} ({} != {})".format(lhs, rhs, l, r))
            missing_dim   = sp.simplify(r/l)
            missing_quant = self._possible_quant(missing_dim,is_print=False)
            
                
            if missing_quant:
                missing_term = self._quant2str(missing_quant[0])
                print("According to the configuration file, there is a direct connection on two sides: {} * {} = {},".format(lhs, missing_term, rhs))
                print("where {} = {}.".format(missing_term, missing_dim))
            else:
                missing_dim = 1/missing_dim
                missing_quant = self._possible_quant(missing_dim,is_print=False)
                if missing_quant:
                    missing_term = self._quant2str(missing_quant[0])
                    print("According to the configuration file, there is a direct connection on two sides: {}  = {} * {},".format(lhs, missing_term, rhs))
                    print("where {} = {}.".format(missing_term, missing_dim))

            return False

    def possible_quant(self, s, is_include_alias: bool=False):
        return self._possible_quant(dim=self.dimension(s), s=s, is_include_alias=is_include_alias)
        
    def display_all_expr(self, is_include_alias: bool=False):
        print("All Expression:")
        for expr,quants in self._expr_dict.items():
            print("{} {}".format(str(expr).center(60," "), self._quant2str(quants, is_include_alias=is_include_alias)))

    def display_all_quant(self, is_include_alias: bool=False):
        print("All Quantities:")
        for i in range(len(self._lhs_list)):
            quant=self._lhs_list[i]
            discription=self._discription_list[i]
            var=self._strip_sqr_bkt(self._lhs_list[i])
            expr=self._drv_unit_dict[var]
            if is_include_alias and self._quant2alias_map[quant]:
                print("{} {} {}".format(quant.center(30," "), discription.center(60," "), str(expr)))
                print("{}".format(("A" +str(self._quant2alias_map[quant])).center(30, " ")))
            else:
                print("{} {} {}".format(quant.center(30," "), discription.center(60," "), expr))


    
    def _wrap_quant(self,quant):
        if isinstance(quant,str):
            return "_dimcheck_quant_"+quant
        else:
            wraped_list=[]
            for x in quant:
                wraped_list.append("_dimcheck_quant_"+x)
            return wraped_list
        

    def _unwrap_quant(self,quant):
        if isinstance(quant,str):
            if self._is_base_unit(quant):
                return quant
            else:
                return "["+quant.replace("_dimcheck_quant_","")+"]"
        else:
            unwraped_list=[]
            for x in quant:
                if self._is_base_unit(x):
                    s=x
                else:
                    s="["+x.replace("_dimcheck_quant_","")+"]"
                unwraped_list.append(s)
            return unwraped_list
    

    def _quant2str(self, quant, is_include_alias=False):
        ## translate a quantity or some quantities to the string format with its/their alias
        if isinstance(quant,str):
            quant_s = self._unwrap_quant(quant)
            if is_include_alias and self._quant2alias_map[quant_s]:
                return quant_s + " A" +str(self._quant2alias_map[quant_s])
            else:
                return quant_s
        else:
            s = ""
            for x in self._unwrap_quant(quant):
                if is_include_alias and self._quant2alias_map[x]:
                    s += x + " A"+ str(self._quant2alias_map[x]) + " "
                else:
                    s += x + " "
            return s

    def _is_base_unit(self, s: str):
        if s in self.base_unit_dict:
            return True
        else:
            return False
        
    def _check_and_replace_alias(self, quant):
        real_quant=self._alias_map.get(quant)
        if real_quant:
            return real_quant
        else:
            return quant


    def _strip_sqr_bkt(self, s: str, is_defining_symbol: bool = False):
        def check_quant_in_symbols(quant,origin_s):
            if not self._is_in_existing_symbols(quant):
                raise UnknowQuantFoundError("Unknown quantity {} is found in {}".format(quant,origin_s))

        origin_s=copy.copy(s)

        ## Find [m],[eps_0]... form variables.
        quants=re.findall(r"\[\s*[a-zA-Z][\w_]*\s*\]",s)
        for quant in quants:
            if not is_defining_symbol:
                check_quant_in_symbols(quant,origin_s)
            tmp=self._wrap_quant(self._check_and_replace_alias(quant).strip("[]").strip())
            s=s.replace(quant,tmp)

        ## Find [sgm*c/E]... form quantiables.
        combined_quants=re.findall(r"\[.*?\]",s)
        for combined_quant in combined_quants:
            new_quants=combined_quant.strip("[]").strip()
            quants=re.findall(r"\b[a-zA-Z][\w_]*\b",new_quants)
            for quant in quants:
                if not is_defining_symbol:
                    if not is_defining_symbol:
                        check_quant_in_symbols("["+quant+"]",origin_s)
                tmp=self._wrap_quant(self._check_and_replace_alias("["+quant+"]").strip("[]").strip())
                new_quants=re.sub(r"\b"+quant+r"\b",tmp,new_quants)
            new_quants="("+new_quants+")"
            s=s.replace(combined_quant,new_quants)


        ## Check whether brakets match or not.
        bkts=re.findall("[\[\]]",s)
        if bkts:
            raise ExprBktError("Brakets does not match for {}".format(origin_s))
        
        ## Check brakets for quantities
        quants=re.findall(r"\b[a-zA-Z][\w_]*\b",s)
        for quant in quants:
            if not self._is_base_unit(quant):
                raise MissingBktError("Missing a braket for the quantity {}".format(quant))
            
        return s

        

        


    def _possible_quant(self, dim: Quantity, s: str="", is_print=True, is_include_alias: bool=False):
        quant=self._expr_dict.get(dim)
        if quant:
            if is_print:
                print("Possible quantity for {} ({}): {}".format(s, dim, self._quant2str(quant, is_include_alias=is_include_alias)))
            return quant
        else:
            if is_print:
                print("No proper quantity found for {} ({}).".format(s, dim))
            return None 


    def _append_formulas(self,quant,rhs,discription):
        if not self._is_in_existing_symbols(quant):
            self._quant_set.add(quant)
        else:
            raise MultiSymDefError("The symbol {} has already been defined before, please choose another symbol!".format(quant))
        self._lhs_list.append(quant)
        self._rhs_list.append(rhs)
        self._discription_list.append(discription)    

    def _append_alias(self,alias,quant):
        self._check_sym_format(alias)
        if not self._is_in_existing_symbols(alias):
            self._alias_map[alias]=quant
        else:
            raise MultiSymDefError("The symbol {} has already been defined before, please choose another symbol!".format(alias))

    def _is_in_existing_symbols(self,quant):
        if quant in self._quant_set or quant in self._alias_map:
            return True
        else:
            return False

    def _check_sym_format(self,s):
        ## Check the format of a quantity or alias
        match_obj=re.fullmatch(r"\[[a-zA-Z][\w_]*\]",s)
        if not match_obj:
            raise InvalidSymFormatError("Please use the correct format for the quantity or alias {} as the regular expression '\[[a-zA-Z][\w_]\]', e.g. '[sigma_xx]', '[m_11]'".format(s))

    
        





# module_path=_split_dir_file(dimcheck.)

## A singleton instance of si unit.
si=Dimcheck(config_file=DIMCHECK_MODULE_PATH+"/si.json",is_save=False)


