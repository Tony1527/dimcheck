'''A module for dimension check.'''


import copy
import json
import re
import hashlib
import pickle
import os
import sympy as sp
import numpy as np
from sympy.physics.units.quantities import Quantity
from sympy import sympify,SympifyError
from queue import Queue
from typing import Any


def _ext_sig_exp(num, accuracy=10):    
    ## Extract the significant digit and the exponent of the number
    if num == 0:
        return (0, 0)
    sign = 1
    if num < 0:
        sign = -1
    num = np.abs(num)
    base_digit = np.floor(np.log10(num))
    significant_digit = num / 10**base_digit

    significant_digit = np.round(significant_digit * 10**accuracy) / 10**accuracy
    return (sign * significant_digit, base_digit)

def _keep_sig(num, accuracy=10):
    ## Keep the significant digit of the number
    significant_digit, base_digit = _ext_sig_exp(num, accuracy)
    return significant_digit * 10**base_digit

def _equal_to_abs(x, y, accuracy=1e-10, abs_zero=1e-14):
    ## Check whether the difference between x and y is smaller than accuracy
    if abs(x) < abs_zero:
        x = 0

    if abs(y) < abs_zero:
        y = 0

    if y != 0:
        if abs((x - y) / y) > accuracy:
            return False
    elif abs((x - y)) > accuracy:
        return False
    return True

def _split_dir_file(s):
    ## Split the directory and filename
    dir_end = s.rfind(os.sep)
    if dir_end != -1:
        return s[:dir_end], s[dir_end + 1 :]
    else:
        return None, s

def _divide_file_path(path):
    ## Divide the file path into directory path, filename and postfix
    dir_path,filename = _split_dir_file(path)
    filename_parts=filename.split(".")
    if len(filename_parts)!=2:
        print("filename ({}) is not in XXX.xxx form.".format(filename))
    filename=filename_parts[0]
    postfix=filename_parts[-1]
    return {"filename":filename,"dir_path":dir_path,"postfix":postfix}

## The global settings
_DIMCHECK_MODULE_PATH = __file__
_dimcheck_path = _split_dir_file(_DIMCHECK_MODULE_PATH)[0]
with open(_dimcheck_path+os.sep+"setting.json") as f:
    json_str=json.load(f)
    if json_str["output"][-1]=="/":
        _DIMCHECK_DEFAULT_OUTPUT_DIR = json_str["output"]
    else:
        _DIMCHECK_DEFAULT_OUTPUT_DIR = json_str["output"] + os.sep
    _DIMCHECK_SI_INSTANCE = json_str["si_instance"]
    _DIMCHECK_IS_SAVE = json_str["is_save"]



class ConflictSymbolDefinitionError(Exception):
    '''Exception raised when users choose the same symbol with different definitions.
    '''
    pass

class SquareBracketMismatchError(Exception):
    '''Exception raised when the left and right braket in the expression does not match.
    '''
    pass

class MissingSquareBracketError(Exception):
    '''Exception raised when users input a quantity without a braket.
    '''
    pass

class QuantityModeError(Exception):
    '''Exception raised when a square bracket appears around quantities in the quantity mode.
    '''
    pass



class MissingParenthesesError(Exception):
    '''Exception raised when users input a e.g. "a**1/2" without a parentheses. Since in python, the priority of "**" is higher than "/", so the expression will be parsed as "(a**1)/2".
    Please use "a**(1/2)" instead.
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

class DerivingQuantityError(Exception):
    '''Exception raised when the quantities cannot be totally derived.
    '''
    pass

class UnknownQuantityError(Exception):
    '''Exception raised when a quantity not defined is found.
    '''
    pass


class InvalidSyntaxError(Exception):
    '''Exception raised when users input the expression with the incorrect syntax.
    '''
    pass











class Dimcheck:
    '''Create a Dimcheck object based on the quant_def_file. This is the interface for users to use.
    '''
    def __init__(self, quant_def_file="./si.json", unit_system="si", is_save=False, setting_file=""):
        '''Initialize the object with quantity definition file

        Parameters
        ----------
        quant_def_file : str, optional
            Filepath of definitions of quantities, by default "/path/to/your/pip/lib/dimcheck/si.json"
        unit_system : str, optional
            Environment of the unit system. So far only "si" is allowed, by default "si"
        is_save : bool, optional
            Whether to save the object in a binary form or not. A hash method is used to check the quantity definition file , by default False
        '''
        self._quant_def_file = quant_def_file
        self._unit_system = unit_system
        self._is_save = is_save
        self._serialized_file = None
        self._is_pretty = False
        self._is_quant = True
        self._setting_file=setting_file

        serialized_dimchecker, is_pretty, is_quant, self._serialized_file, self._hash = self._compare_serialized_hash_match(quant_def_file=quant_def_file)
        if serialized_dimchecker:
            self._dimchecker = serialized_dimchecker
            self._is_save = False
            self._is_pretty = is_pretty
            self._is_quant = is_quant
        else:
            self._dimchecker = Dimchecker(quant_def_file=quant_def_file, unit_system=unit_system)

        if self._is_save:
            self.save()

    @property
    def is_pretty(self) -> bool:
        '''
        The output is setting to be more easy to read (based on some symbols in UTF-8) or not.
        '''
        return self._is_pretty      
    
    @is_pretty.setter
    def is_pretty(self, is_pretty: bool) -> None:
        '''
        Set the output to be pretty or not.
        '''
        self._is_pretty = is_pretty


    @property
    def is_quant(self) -> bool:
        '''
        If True, then the symbols that users input will be regarded as quantities, then no square bracket is needed anymore.
        '''
        return self._is_quant      
    
    @is_quant.setter
    def is_quant(self, is_quant: bool) -> None:
        '''
        Set the output to be pretty or not.
        '''
        self._is_quant = is_quant


    @property
    def quant_def_file(self) -> str:
        """
        Returns the path of the file of quantity definition.
        """
        return self._quant_def_file
    
    @property
    def setting_file(self) -> str:
        """
        Returns the path of the global setting file.
        """
        return self._setting_file
    
    
    @property
    def serialized_file(self) -> str:
        '''
        Returns the path of the serialized file.
        '''
        return self._serialized_file
    
    @property
    def unit_system(self) -> str:
        '''
        Returns the unit system.
        '''
        return self._unit_system

    @property
    def hash(self) -> str:
        '''
        Returns the hash of the quantity definition file.
        '''
        return self._hash

    def unit(self,s: str) -> str:
        '''Get the unit of the input

        Parameters
        ----------
        s : str
            Any quantity or the combination of the quantities

        Returns
        -------
        expr : str
            unit of the input
        '''
        return str(self._dimchecker.unit(s,is_pretty=self._is_pretty,is_quant=self._is_quant))
    
    def dimension(self,s: str) -> str:
        '''Get the dimension of the input

        Parameters
        ----------
        s : str
            Any quantity or the combination of the quantities

        Returns
        -------
        expr : str
            unit of the input
        '''
        return str(self._dimchecker.dimension(s,is_pretty=self._is_pretty,is_quant=self._is_quant))
    
    def dim(self,s: str) -> str:
        '''Alias of Dimcheck.dimension.
        
        Get the dimension of the input.

        Parameters
        ----------
        s : str
            Any quantity or the combination of the quantities

        Returns
        -------
        expr : str
            dimension of the input
        '''
        return self.dimension(s)

    def is_dc(self,lhs: str, rhs: str, is_print: bool=False) -> bool:
        '''To see whether the quantities on the two sides are equal or not.

        If return False, a straightforward relation based on the quantity definition file will be checked. 

        Parameters
        ----------
        lhs : str
            A quantity or a quantity combinationon on the left hand side
        rhs : str
            A quantity or a quantity combinationon on the right hand side
        is_print : bool
            Whether to print the process or not

        Returns
        -------
        bool
            Return True if the dimensions on two sides are equal, otherwise return False.
        '''
        return self._dimchecker.is_dc(lhs,rhs, is_print=is_print, is_pretty=self._is_pretty,is_quant=self._is_quant)

    def quant(self,s: str, is_print=False) -> Any:
        '''Print the possible quantity by deriving the combination of quantities.
        '''
        quant = self._dimchecker.quant(s, is_print=is_print, is_pretty=self._is_pretty,is_quant=self._is_quant)
        if quant:
            if len(quant)==1:
                return quant[0]
            else:
                return quant
        else:
            return quant
        
    def all_expr(self, is_inc_alias: bool=False, is_sorted: bool=True):
        '''Display all expressions based on the quantity definition file.
        '''
        self._dimchecker.display_all_expr(is_inc_alias=is_inc_alias, is_pretty=self._is_pretty, is_sorted=is_sorted, is_quant=self._is_quant)

    def all_quant(self, is_inc_alias: bool=False, is_sorted: bool=True):
        '''Display all quantities based on the quantity definition file.
        '''
        self._dimchecker.display_all_quant(is_inc_alias=is_inc_alias, is_pretty=self._is_pretty, is_sorted=is_sorted, is_quant=self._is_quant)

    def save_all_quant(self, file_path=_DIMCHECK_DEFAULT_OUTPUT_DIR+"Quantities.csv", is_sorted: bool=True) -> bool:
        '''Save all quantities in ".csv" format based on the quantity definition file.

        Parameters
        ----------
        file_path : str, optional
            so far, only ".csv" format is allowed, by default "./Quantities.csv"
        '''
        return self._dimchecker.save_all_quant(file_path=file_path, is_pretty=self._is_pretty, is_sorted=is_sorted, is_quant=self._is_quant)

    def save_all_expr(self, file_path=_DIMCHECK_DEFAULT_OUTPUT_DIR+"Expressions.csv", is_sorted: bool=True) -> bool:
        '''Save all expression in ".csv" format based on the quantity definition file.

        Parameters
        ----------
        file_path : str, optional
            so far, only ".csv" format is allowed, by default "./Quantities.csv"
        '''
        return self._dimchecker.save_all_expr(file_path=file_path, is_pretty=self._is_pretty, is_sorted=is_sorted, is_quant=self._is_quant)

    def save(self) -> bool:
        '''Serialize the current object
        '''
        return self._serialize()

    def clean(self) -> bool:
        '''Clean the current object.
        '''
        if os.path.exists(self._serialized_file):
            self._hash = ""
            os.remove(self._serialized_file)
            return True
        else:
            print("{} does not exist!".format(self._serialized_file))
            return False

    def omit_quant(self, lhs: str, rhs: str, omit_quant: list) -> str:
        '''Restore the units of lhs based on the omitted quantities and rhs.

        Parameters
        ----------
        lhs : str
            A quantity or a quantity combinationon on the left hand side
        rhs : str
            A quantity or a quantity combinationon on the right hand side
        omit_quant : list
            The omitted quantities, by default [], you can also input a single quantity in the string format

        Returns
        -------
        expr : str
            The formula based on the rhs and omit_quant
        '''
        return self._dimchecker.omit_quant(lhs=lhs, rhs=rhs, omit_quant=omit_quant, is_pretty=self._is_pretty,is_quant=self._is_quant)

    def formula(self, lhs: str, parameters: list) -> str:
        '''Print the formula of the target based on the parameters.

        Parameters
        ----------
        lhs : str
            The target quantity or quantity combination
        parameters : list,
            parameters of the system
        
        Returns
        -------
        expr : str
            The formula based on the parameters
        '''
        return self.omit_quant(lhs=lhs, rhs="1", omit_quant=parameters)

    def _serialize(self) -> None:
        ## Serialize the current object
        self._hash = self._cal_quant_def_file_hash(self._quant_def_file)
        data = {}
        with open(self._serialized_file,"wb") as f:
            data["unit_system"] = self._unit_system
            data["quant_def_file"] = self._quant_def_file
            data["is_pretty"] = self._is_pretty
            data["is_quant"] = self._is_quant

            data["hash"] = self._hash
            data["dimcheck_obj"] = self._dimchecker
            
            pickle.dump(data,f)
            return True
        

    def _compare_serialized_hash_match(self,quant_def_file: str):
        ## Compare the hash of the quantity definition file and the serialized object
        serialized_dimcheck_obj = None
        is_pretty = False
        is_quant = False
        quant_def_file_hash = ""
        pickle_file = self._get_serialized_path(quant_def_file=quant_def_file)
        
        if os.path.exists(pickle_file):
            try:
                with open(pickle_file,"rb") as f:
                    data = pickle.load(f)
                    serialized_hash = data["hash"]
                    is_pretty = data["is_pretty"]
                    is_quant = data["is_quant"]
                    quant_def_file_hash = self._cal_quant_def_file_hash(quant_def_file)
                    if serialized_hash == quant_def_file_hash:
                        serialized_dimcheck_obj = data["dimcheck_obj"]
            except Exception as e:
                print("The serialized file {} is outdated! It will remove the outdated data automatically.".format(pickle_file))
                os.remove(pickle_file)
        else:
            quant_def_file_hash = ""
        
        return serialized_dimcheck_obj, is_pretty, is_quant, pickle_file, quant_def_file_hash
    
    def _cal_quant_def_file_hash(self, quant_def_file, hash_method=hashlib.md5):
        ## Calculate the hash of the quantity definition file
        digestobj=None
        with open(quant_def_file,"r") as f:
            digestobj=hash_method()
            while True:
                data = f.read(4096)
                if not data:
                    break
                digestobj.update(data.encode("utf-8"))

        return digestobj.hexdigest()

    def _get_serialized_path(self, quant_def_file: str):
        ## Get the path of the serialized object
        file_info=_divide_file_path(quant_def_file)
        filename=file_info["filename"]
        dir_path=file_info["dir_path"]
        pickle_file=dir_path+os.sep+filename+".pickle"
        return pickle_file




class Dimchecker():
    class Unit():
        def __init__(self, id: str, alias: str, expr: str, quant_form: str, definition: str, discription: str, dim_array: np.ndarray):
            self._id = id
            self._dim_array  = dim_array
            self._quant_form = quant_form
            self._definition = definition
            self._discription = discription
            self._alias = alias
            self._expr = expr
            


        def __str__(self):
            return self._quant_form
        
        @property
        def quant_form(self):
            return self._quant_form
        
        @property
        def id(self):
            return self._id
        
        @property
        def alias(self):
            return self._alias
        
        @property
        def expr(self):
            return self._expr
        
        @property
        def dim_array(self):
            return self._dim_array
        
        @property
        def definition(self):
            return self._definition
        
        @property
        def discription(self):
            return self._discription
        
        def print(self):
            print(self._id, 
                self._quant_form,
                self._alias,
                self._expr,
                self._dim_array,
                self._definition,
                self._discription)
    

    ''' Initialize a dimensional checker with a certain quantity definition file.
    '''
    def __init__(self,quant_def_file="./si.json",unit_system="si"):
        '''Initialize a dimensional checker with a certain quantity definition file.

        Parameters
        ----------
        quant_def_file : str, optional
            quantity definition file.
        unit_system : str, optional
            environment, by default "si". So far only "si" is allowed.

        Raises
        ------
        DerivingQuantityError
            Raise an error if the quantities cannot be totally derived.
        '''
        
        self._unit_system=unit_system
        self._quant_def_file=quant_def_file

        ## Definition of the symbols
        self._lhs_list=[]
        self._rhs_list=[]
        self._discription_list=[]

        self._dimcheck_quant_list=[]

        ## Filter the repeated definition
        self._quant_set=set()

        ## Two mutual maps between a quantity and its alias
        self._alias2quant_map={}
        self._quant2alias_map={}
        
        kg, m, s, A, cd, K, mol  = sp.symbols("kg, m, s, A, cd, K, mol",positive=True)
        M, L, T, I, J, Theta, N  = sp.symbols("M, L, T, I, J, Theta, N",positive=True)

        self._base_unit_dict={"kg":kg, "m":m, "s":s, "A":A, "cd":cd, "K":K, "mol":mol}
        self._unit2dimension_str_dict={"kg":"M", "m":"L", "s":"T", "A":"I", "cd":"J", "K":"Theta", "mol":"N"}
        self._unit2dimension_sym_dict={kg:M, m:L, s:T, A:I, cd:J, K:Theta, mol:N}
        self._dimension_dict={"M":M, "L":L, "T":T, "I":I, "J":J, "Theta":Theta, "N":N}

        # self._base_unit_dict={"kg":"kg", "m":"m", "s":"s", "A":"A", "cd":"cd", "K":"K", "mol":"mol"}
        self._base_unit_order={"kg":0,"m":1,"s":2,"A":3,"cd":4,"K":5,"mol":6}
        self._base_unit_num=len(self._base_unit_order)

        self._keywords = ["sqrt","cbrt"]




        self._superscripts = {
            "1": "¹",
            "2": "²",
            "3": "³",
            "4": "⁴",
            "5": "⁵",
            "6": "⁶",
            "7": "⁷",
            "8": "⁸",
            "9": "⁹",
            "0": "⁰",
            "-": "⁻",
            ".": "·",
            "/": "ʴ"    #"ʹ"
        }
        

        ## Load the quantity definition file and check the repeated definitions of symbols
        with open(quant_def_file,"r") as f:
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
                
        

        print("Starting to parse the quantity definition file: \n{}".format(self._quant_def_file))
        dependence_queue=Queue()
        unsolved_num=0

        ## {A quantity:An expression}
        self._drv_unit_dict=copy.deepcopy(self._base_unit_dict)


        ## {An expression:A quantity or a quantity combination(which may have multiple symbols)}. E.g. m**(-1):["[k]","[nabla]"]
        self._expr_dict={}
        self._dim_arr_dict={}

        ## Defining symbols according to the quantity definition file
        for i in range(len(self._lhs_list)):
            self._check_sym_format(self._lhs_list[i])
            striped_str=self._strip_sqr_bkt(self._rhs_list[i],is_def_sym=True)

            ## Substute the symbols according to self._drv_unit_dict and simplify the whole expression
            try:
                expr=sp.powsimp(sympify(striped_str).subs(self._drv_unit_dict))
            except SympifyError or SympifyError:
                raise InvalidSyntaxError("Invalid syntax: {}".format(self._rhs_list[i]))
            quants=re.findall("_dimcheck_quant_\w+",str(expr))
            tmp=[]
            for quant in quants:
                if not self._drv_unit_dict.get(quant):
                    tmp.append(quant)
            dependence_queue.put((expr,i,tmp))
            
            ## Assign uncompleted expression for a quantity
            self._drv_unit_dict.update({self._strip_sqr_bkt(self._lhs_list[i],is_def_sym=True):expr})




        ## Use BFS algorithm to find all undefined quantities
        max_iter=99
        round=0
        for r in range(max_iter+1):
            round=r
            unsolved_num=dependence_queue.qsize()
            for i in range(unsolved_num):
                num=dependence_queue.qsize()
                for j in range(num):
                    expr,idx,unsolved_quants=dependence_queue.get()

                    ## Set is_solvable to False if the de
                    is_solvable=True
                    for quant in unsolved_quants:
                        testing_expr = self._drv_unit_dict.get(quant)
                        
                        if testing_expr is None:
                            is_solvable=False
                        else: 
                            uncomplete_quants = re.findall("_dimcheck_quant_\w+",str(testing_expr))
                            if len(uncomplete_quants)!=0:
                                is_solvable=False

                    if is_solvable:
                        try:
                            expr=sp.powsimp(expr.subs(self._drv_unit_dict))
                        except SympifyError:
                            raise InvalidSyntaxError("Invalid syntax: {}".format(self._drv_unit_dict))
                        self._drv_unit_dict.update({self._strip_sqr_bkt(self._lhs_list[idx],is_def_sym=True):expr})
                    else:
                        dependence_queue.put((expr,idx,unsolved_quants))


        ## Check the dependence_queue, if there is any unsolve quantity, then raise an error
        unsolved_num=dependence_queue.qsize()
        if unsolved_num!=0:
            for i in range(unsolved_num):
                expr,idx,unsolved_quants=dependence_queue.get()

            error_msg=""
            if round==max_iter:
                error_msg+="Exceed maximum recursive depth! There is probably a CIRCULAR DEFINITION!\n"

            error_msg += "The quantities {} can not be totally derived!".format(self._quant2str(unsolved_quants))
            
            raise DerivingQuantityError(error_msg)

        
        for quant,expr in self._drv_unit_dict.items():
            if self._is_base_unit(quant):
                continue
            expr_s = str(expr)
            quants=self._expr_dict.get(expr_s)
            dim_arr=self._generate_dim_array(expr)
            if not quants:
                self._expr_dict[expr_s]=[quant]
                self._dim_arr_dict[dim_arr.tobytes()]=[quant]
            else:
                self._expr_dict[expr_s].append(quant)
                self._dim_arr_dict[dim_arr.tobytes()].append(quant)


        # print(self._expr_dict)
        # print(self._drv_unit_dict)


        ## self._all_units: {id:unit}; the map to all defined units
        self._all_units=self._generate_all_units()
        # for key,values in self._all_units.items():
        #     values.print()
        
        print("Successfully parsed!")


    def _save_default_str(self, x, default=""):
        ## Save the default string for the csv file
        if (not isinstance(x, str)) and x!=None and hasattr(x, "__iter__"):
            s=""
            for i in range(len(x)):
                v=x[i]

                if i==0:
                    s += "{}".format(v)
                else:
                    s += "  {}".format(v)
            return s
        elif x:
            return str(x).replace(",",".")
        else:
            return default

    def save_all_quant(self, file_path: str, is_pretty: bool=False, is_sorted: bool=True, is_quant: bool=False) -> bool:
        ## Save all quantities in the csv file
        with open(file_path, mode="w", encoding="utf-8") as f:
            f.write("{},{},{},{},{},{}".format("Quantity", "Discription", "Definition", "Unit", "Alias", os.linesep))
            quant_list = []
            if is_sorted:
                quant_list = sorted(self._dimcheck_quant_list)
            else:
                quant_list = self._dimcheck_quant_list
            
            for quant in quant_list:
                one_line = ""
                unit = self._all_units[quant]
                if is_pretty:
                    expr_s = self._pretty_sp_expr(unit.expr)
                    quant_s = self._pretty_sp_expr(unit.quant_form)
                    def_s = self._pretty_sp_expr(unit.definition)
                    pretty_aliases = []
                    if isinstance(unit.alias,list):
                        for alias in unit.alias:
                            pretty_aliases.append(self._pretty_sp_expr(alias))
                    else:
                        pretty_aliases = self._pretty_sp_expr(unit.alias)
                    alias_s = pretty_aliases
                else:
                    expr_s = unit.expr
                    quant_s = unit.quant_form
                    alias_s = unit.alias
                    def_s = unit.definition

                
                # if is_quant:
                #     quant_s = quant_s.strip("[]")
                #     if alias_s:
                #         alias_s = [a.strip("[]") for a in alias_s]
                one_line += "{},{},{},{},{},{}".format(self._save_default_str(quant_s), self._save_default_str(unit.discription), self._save_default_str(def_s), self._save_default_str(expr_s), self._save_default_str(alias_s),os.linesep)
                f.write(one_line)

        return True


    def save_all_expr(self, file_path: str, is_pretty: bool=False, is_sorted: bool=True, is_quant: bool=False) -> bool:
        ## Save all expressions in the csv file
        with open(file_path, mode="w", encoding="utf-8") as f:
            f.write("{},{},{}".format("Unit", "Quantity", os.linesep))
            keys = []
            if is_sorted:
                keys = sorted(self._expr_dict)
            else:
                keys = self._expr_dict.keys()
            for expr in keys:
                quants = self._expr_dict[expr]
                one_line = ""
                
                if is_pretty:
                    expr_s = self._pretty_sp_expr(expr)
                    pretty_quants = []
                    for quant in quants:
                        pretty_quants.append(self._pretty_sp_expr(quant))
                    quants_s = self._quant2str(pretty_quants, is_quant=False)
                else:
                    expr_s = expr
                    quants_s = self._quant2str(quants, is_quant=False)
                f.write("{},{},{}".format(self._save_default_str(expr_s), self._save_default_str(quants_s), os.linesep))
                f.write(one_line)

        return True

    def omit_quant(self, lhs: str, rhs: str, omit_quant: list=[], is_pretty = False, is_quant: bool=False):
        def is_solvable(A, b):
            x=np.dot(np.linalg.pinv(A),b)
            if not np.allclose(np.dot(A,x),b):
                return False
            else:
                return True
            
        def search_min_array(chosen_units, remained_units, b):
            r=copy.copy(remained_units)
            rtn=None
            idx_list, chosen_units_list= chosen_units
            for i in range(len(r)):
                idx,v=r.pop()
                new_idx_list=idx_list+[idx]
                tmp_array = chosen_units_list+[v]
                A=np.vstack(tmp_array).T
                if len(tmp_array)<len(omit_quant) and is_solvable(A,b):
                    rtn = (new_idx_list,tmp_array)
                    break

                rtn = search_min_array(chosen_units = (new_idx_list,tmp_array),remained_units=r, b=b)
                if rtn!=None:
                    break
            return rtn

        

        ## Restore the units of lhs based on the omitted quantities.
        ldim=str(self.unit(lhs,is_pretty=False, is_quant=is_quant))
        rdim=str(self.unit(rhs,is_pretty=False, is_quant=is_quant))

        ldim_array=self._generate_dim_array(ldim)
        rdim_array=self._generate_dim_array(rdim)
        units_array=[]
        
        if isinstance(omit_quant,str):
            omit_quant = [omit_quant]
        
        for omitted_unit in omit_quant:
            units_array.append(self._generate_dim_array(self.unit(omitted_unit, is_pretty=False, is_quant=is_quant)))

        ## in the Transpose form Ax=b, x=pinv(AT).b
        A=np.vstack(units_array)
        b=ldim_array - rdim_array
        
        mask_A=A!=0
        mask_b=b!=0
        l=[np.any((mask_A==mask_b)[:,i]) for i in range(len(b))]

        ## Check if it is solvable
        for i,value in enumerate(mask_b):
            if value and l[i]==False:
                print("Based on the omitted quantities, no formula between lhs and rhs is found.")
                return None
            
        min_units_array=[]
        min_idx_list=[]
        min_A=None
        
        min_rtn=search_min_array(([],[]),remained_units=[(i,v) for i,v in enumerate(units_array)],b=b)
        if min_rtn!=None:
            min_idx_list, min_units_array = min_rtn
            min_omit_quant=[]
            for idx in min_idx_list:
                min_omit_quant.append(omit_quant[idx])
            min_A=np.vstack(min_units_array).T
            # print(min_idx_list, min_units_array)
        else:
            min_units_array = units_array
            min_idx_list = []
            min_omit_quant = omit_quant
            min_A=A.T
        

        if len(min_units_array)<len(units_array):
            print("Redudant units are found.")
    
        x=np.dot(np.linalg.pinv(min_A),b)
        
        is_rhs_unit_s=False
        if rhs == "1":
            is_rhs_unit_s=True
            s="{} = ".format(lhs)
        else:
            s="{} = {}".format(lhs,rhs)

        is_1st_nonzero=True
        for i in range(len(min_omit_quant)):
            omitted_unit = min_omit_quant[i]
            exp = _keep_sig(x[i])
            single_word_match=re.match(r"^[a-zA-Z][\w_]*\r?$",omitted_unit)

            if exp == int(exp):
                exp = int(exp)
                exp_str = str(exp)
            else:
                exp_str = "({:.9})".format(exp)

            if single_word_match or not is_quant:
                omitted_unit_str = omitted_unit
            else:
                omitted_unit_str = "("+omitted_unit+")"

            

            
            if _equal_to_abs(exp,0):
                continue
            elif _equal_to_abs(exp,1):
                if is_1st_nonzero and is_rhs_unit_s:
                    is_1st_nonzero=False
                else:
                    s+="*"
                s+="{}".format(omitted_unit)    
            else:
                if is_1st_nonzero and is_rhs_unit_s:
                    is_1st_nonzero=False
                else:
                    s+="*"
                s+="{}**{}".format(omitted_unit_str,exp_str) 
                
        # print("Based on the omitted quantities, a formula between lhs and rhs is found:")
        # print(s)
                    
        if is_pretty:
            s=self._pretty_sp_expr(s)

        return s
        # else:
        #     print("Based on the omitted quantities, no formula between lhs and rhs is found.")
        #     return None
        # else:
        #     common_factor=0
        #     is_linear_independent=False
        #     omitted_unit_array=units_array[0]
        #     for i in range(len(ldim_array)):
        #         lx = omitted_unit_array[i]
        #         rx = rdim_array[i]-ldim_array[i]
        #         if lx==0 and rx==0:
        #             continue
        #         elif (lx==0 and rx!=0) or (lx!=0 and rx==0):
        #             is_linear_independent=True
        #             break
        #         else:
        #             factor = rx/lx
        #             if common_factor==0:
        #                 common_factor = factor
        #             else:
        #                 if not np.allclose(np.array([common_factor]), np.array([factor])):
        #                     is_linear_independent=True
        #                     break

        #     if is_linear_independent:
        #         print("Linear indepent")
        #     else:
        #         print("{}*{}**({})={}".format(lhs,omit_quant[0],int(common_factor),rhs))


    def _generate_all_units(self):
        ## Generate all units based on the variables that are extracted from the quantity definition file
        units={}
        for i in range(len(self._lhs_list)):
            quant = self._lhs_list[i]
            definition = self._rhs_list[i]
            discription = self._discription_list[i]
            var = self._strip_sqr_bkt(self._lhs_list[i])
            self._dimcheck_quant_list.append(var)
            expr = str(self._drv_unit_dict[var])
            dim_array = self._generate_dim_array(expr)
            units.update({var:Dimchecker.Unit(id=var,quant_form=quant,alias=self._quant2alias_map[quant],expr=expr,definition=definition,discription=discription,dim_array=dim_array)})

        return units

    def _generate_dim_array(self,expr):
        ## Generate the dimension array based on the expression
        def convert_dim_array(s):
            s=s.strip("-")
            variables = s.replace("**","^").split("*")
            total_dim_arr=np.zeros(self._base_unit_num, dtype=np.double)
            for variable in variables:
                terms=variable.split("^")
                dim_arr=np.zeros(self._base_unit_num, dtype=np.double)
                if len(terms)==2:
                    unit,exp=terms
                    exp = exp.strip("()")
                    float_exp = float(exp)
                    int_exp = int(float_exp)
                    if float_exp == int_exp:
                        exp = int_exp
                    else:
                        raise DerivingQuantityError("When comparing the dimension of two quantities, the exponent of the dimension must be an integer rather than {}.".format(s))
                    dim_arr[self._base_unit_order[unit]]=1*exp
                elif len(terms)==1:
                    unit=terms[0]
                    if unit=="1":
                        continue
                    if not unit in self._base_unit_order:
                        raise DerivingQuantityError("the dimension array of {} can't be obtained.".format(s))    
                    dim_arr[self._base_unit_order[unit]]=1
                else:
                    raise DerivingQuantityError("the dimension array of {} can't be obtained.".format(s))
                total_dim_arr+=dim_arr
            return total_dim_arr


        dim_array=None
        expr = str(sp.powsimp(sp.sympify(expr).subs(self._drv_unit_dict)))

        ## self.keywords is required to be checked first
        keywords_obj = re.search(r"(cbrt|sqrt)",expr)
        pow_obj = re.search(r"(\*\*\(-*\d+\/\d+\)|\^\(-*\d+\/\d+\))",expr)
        if pow_obj or keywords_obj:
            raise DerivingQuantityError("When comparing the unit of two quantities, the exponent of the unit must be an integer rather than {}.".format(expr))

        parts=expr.split("/")
        if len(parts)==2:
            dim_array = convert_dim_array(parts[0].strip("()"))-convert_dim_array(parts[1].strip("()"))
        elif len(parts)==1:
            dim_array = convert_dim_array(parts[0].strip("()"))
        elif len(parts)==0:
            dim_array = convert_dim_array(parts[0].strip("()"))
        else:
            raise InvalidSymFormatError("the unit array of {} can't be obtained.".format(expr))
        
        return dim_array

    def unit(self,s: str, is_pretty: bool=False, is_quant: bool=False):
        ## Get the unit of the input
        striped_str=self._strip_sqr_bkt(str(s),is_quant=is_quant)

        try:
            expr=str(sp.powsimp(sympify(striped_str).subs(self._drv_unit_dict)))
        except SympifyError:
            raise InvalidSyntaxError("Invalid syntax: {}".format(s))
        if is_pretty:
            expr=str(self._pretty_sp_expr(expr))
            
        return expr
    
    def dimension(self,s: str, is_pretty: bool=False, is_quant: bool=False):
        ## Get the dimension of the input
        unit = self.unit(s,is_pretty=is_pretty,is_quant=is_quant)
        for key,value in self._unit2dimension_str_dict.items():
            unit = unit.replace(key,value)
        return unit
    
    def dimension_array(self, expr: str):
        ## Get the dimension array of the input
        return self._generate_dim_array(expr)


    def is_dc(self, lhs, rhs, is_print: bool=True, is_pretty: bool=False, is_quant: bool=False):
        ## To see whether the quantities on the two sides are equal or not.
        l=self.unit(lhs,is_pretty=False, is_quant=is_quant)
        r=self.unit(rhs,is_pretty=False, is_quant=is_quant)

        l_s=None
        r_s=None

        if is_pretty:
            lhs = self._pretty_sp_expr(lhs)
            rhs = self._pretty_sp_expr(rhs)
            l_s = self._pretty_sp_expr(l)
            r_s = self._pretty_sp_expr(r)
        else:
            l_s = l
            r_s = r

        if np.all(self.dimension_array(str(l))==self.dimension_array(str(r))):
            if is_print:
                print("DC!")
                print("{} == {} ({})".format(lhs, rhs, r_s))
            return True
        else:
            if is_print:
                print("There are some differences on two sides. {} != {} ({} != {})".format(lhs, rhs, l_s, r_s))
            missing_dim   = self.unit("("+str(r)+")/("+str(l)+")",is_pretty=False)
            missing_quant = self._quant(missing_dim,is_print=False,is_pretty=is_pretty)
            
                
            if missing_quant:
                missing_term = self._quant2str(missing_quant[0], is_quant=is_quant)
                if is_print:
                    print("According to the quantity definition file, there is a direct connection on two sides: {} * {} = {},".format(lhs, missing_term, rhs))
                    print("where {} = {}.".format(missing_term, self._pretty_sp_expr(missing_dim)))
            else:
                missing_dim = "1/"+str(missing_dim)
                missing_quant = self._quant(missing_dim,is_print=False,is_pretty=is_pretty)
                if missing_quant:
                    missing_term = self._quant2str(missing_quant[0], is_quant=is_quant)
                    if is_print:
                        print("According to the quantity definition file, there is a direct connection on two sides: {}  = {} * {},".format(lhs, missing_term, rhs))
                        print("where {} = {}.".format(missing_term, self._pretty_sp_expr(missing_dim)))

            return False

    def quant(self, s, is_print: bool=False, is_pretty: bool=False, is_quant: bool=False):
        ## Print the possible quantity by deriving the combination of quantities.
        quant_list = self._unwrap_quant(self._quant(dim=self.unit(s,is_pretty=False,is_quant=is_quant), s=s, is_print=is_print, is_pretty=is_pretty), is_quant=is_quant)
        if is_pretty and quant_list:
            for i in range(len(quant_list)):
                quant_list[i] = self._pretty_sp_expr(quant_list[i])
        return quant_list
        
    def display_all_expr(self, is_inc_alias: bool=False, is_pretty: bool=False, is_sorted: bool=True, is_quant: bool=False):
        ## Display all expressions based on the quantity definition file.
        print("All Expression:")
        keys = []
        if is_sorted:
            keys = sorted(self._expr_dict)
        else:
            keys = self._expr_dict.keys()
        for expr in keys:
            quants = self._expr_dict[expr]
            s = self._pretty_sp_expr(str(expr)).ljust(60," ")
            quant = self._quant2str(quants, is_inc_alias=is_inc_alias, is_quant=is_quant)
            if is_pretty:
                s = self._pretty_sp_expr(s)
                quant = self._pretty_sp_expr(quant)
            print("{} {}".format(s, quant))

    def display_all_quant(self, is_inc_alias: bool=False, is_pretty: bool=False, is_sorted: bool=True, is_quant: bool=False):
        ## Display all quantities based on the quantity definition file.
        print("All Quantities:")

        quant_list = []
        if is_sorted:
            quant_list = sorted(self._dimcheck_quant_list)
        else:
            quant_list = self._dimcheck_quant_list
        
        for var in quant_list:
            unit = self._all_units[var]
            quant = unit.quant_form
            discription = unit.discription
            expr = unit.expr
            alias = unit.alias

            if is_quant:
                quant = quant.strip("[]")
                if alias:
                    alias = [a.strip("[]") for a in alias]

            
            if is_pretty:
                expr_s = self._pretty_sp_expr(expr)
            else:
                expr_s = expr
            if is_inc_alias and alias:
                print("{} {} {}".format(quant.ljust(30," "), discription.ljust(60," "), expr_s))
                print("{}".format(("A" +str(alias)).ljust(30, " ")))
            else:
                print("{} {} {}".format(quant.ljust(30," "), discription.ljust(60," "), expr_s))


    def  _pretty_sp_expr(self,sp_expr:Any):
        ## Pretty print the sympy expression
        sp_expr = sp_expr.__str__()
        
        bkt_exp_list=re.findall(r"(\*\*\(-*[0-9]+[\./]{0,1}[0-9]*\)|\^\(-*[0-9]+[\./]{0,1}[0-9]*\))",sp_expr)
        no_bkt_exp_list=re.findall(r"(\*\*-*[0-9]+[\.]{0,1}[0-9]*|\^-*[0-9]+[\.]{0,1}[0-9]*)",sp_expr)

        for exp in bkt_exp_list:
            new_exp=""
            for s in exp:
                if s not in ["(",")","*","^"]:
                    new_exp+=self._superscripts[s]

            sp_expr = re.sub(re.escape(exp),new_exp,sp_expr)
            # sp_expr = re.sub(re.escape(exp)+r"\b",new_exp,sp_expr)
            
        for exp in no_bkt_exp_list:
            new_exp=""
            for s in exp:
                if s not in ["*","^"]:
                    new_exp+=self._superscripts[s]
            sp_expr = re.sub(re.escape(exp),new_exp,sp_expr)
            # sp_expr = re.sub(re.escape(exp)+r"\b",new_exp,sp_expr)
        
        return sp_expr

    
    def _wrap_quant(self,quant):
        ## Wrap a quantity or some quantities with a prefix "_dimcheck_quant_"
        if isinstance(quant,str):
            return "_dimcheck_quant_"+quant
        else:
            wraped_list=[]
            for x in quant:
                wraped_list.append("_dimcheck_quant_"+x)
            return wraped_list
        

    def _unwrap_quant(self,quant,is_quant: bool=False):
        ## Unwrap a quantity or some quantities with a prefix "_dimcheck_quant_"
        if not quant:
            return None
        
        if isinstance(quant,str):
            if self._is_base_unit(quant) or self._is_reserved_keywords(quant):
                return quant
            else:
                if is_quant:
                    return quant.replace("_dimcheck_quant_","")
                else:
                    return "["+quant.replace("_dimcheck_quant_","")+"]"
        else:
            unwraped_list=[]
            for x in quant:
                if self._is_base_unit(x) or self._is_reserved_keywords(quant):
                    s=x
                else:
                    if is_quant:
                        s=x.replace("_dimcheck_quant_","")
                    else:
                        s="["+x.replace("_dimcheck_quant_","")+"]"
                unwraped_list.append(s)
            return unwraped_list
    

    def _quant2str(self, quant, is_inc_alias: bool=False, is_quant: bool=False):
        ## translate a quantity or some quantities to the string format with its/their alias
        if isinstance(quant,str):
            quant_s = self._unwrap_quant(quant)
            if is_inc_alias and self._quant2alias_map[quant_s]:
                alias = self._quant2alias_map[quant_s]
                alias_s = str([a.strip("[]") for a in alias])
                if is_quant:
                    return quant_s.strip() + " A" + alias_s
                else:
                    return quant_s + " A" + alias_s
            else:
                return quant_s
        else:
            s = ""
            for x in self._unwrap_quant(quant):
                
                if is_quant:
                    x_s = x.strip("[]")
                else:
                    x_s = x

                if is_inc_alias and self._quant2alias_map[x]:
                    alias = self._quant2alias_map[x]
                    alias_s = str([a.strip("[]") for a in alias])
                    s += x_s + " A"+ alias_s + " "
                else:
                    s += x_s + " "
            return s

    def _is_base_unit(self, s: str):
        ## Check whether the input is a base unit or not
        if s in self._base_unit_dict:
            return True
        else:
            return False
        
    def _is_reserved_keywords(self, s: str):
        if s in self._keywords:
            return True
        else:
            return False
        
    def _check_and_replace_alias(self, alias):
        ## Check whether the alias is defined or not, if yes, replace it with the quantity
        quant=self._alias2quant_map.get(alias)
        if quant:
            return quant
        else:
            return alias


    def _strip_sqr_bkt(self, s: str, is_def_sym: bool = False, is_quant: bool = False):
        ## Strip the square brakets of the input
        def check_quant_in_symbols(quant,origin_s):
            if not self._is_in_existing_symbols(quant):
                raise UnknownQuantityError("Unknown quantity {} is found in {}".format(quant,origin_s))
            
        def rm_mul_div_space(match):
            return match.group(1)

        origin_s=copy.copy(s)

        prefix_num_s=re.findall("\s*[^_][\d]+[a-zA-Z][\w_]*\s*",s)
        if prefix_num_s:
            raise InvalidSymFormatError("Invalid symbol format {}".format(s))

        invalid_pow_obj = re.search(r"(\*\*\s*-*\d+\s*/\s*\d+\s*|\^\s*-*\d+\s*/\s*\d+\s*)",s)
        if invalid_pow_obj:
            raise MissingParenthesesError("Missing parentheses for {}. Since in Python, the priority of '**' is higher than '/', then e.g. 'a**1/2' will be parsed as '(a**1)/2' in Python.".format(s))

        

        if not is_quant:
            bkts=re.findall("\[\]",s)
            if bkts:
                raise MissingQuantError("Missing a value for the bracket {}".format(s))
            
            ## Find [m][eps_0] form and replace them with [m][eps_0].
            quants=re.findall(r"\]\s*\[",s)
            for quant in quants:
                s=s.replace(quant,"]*[")

            ## Find [m],[eps_0]... form variables.
            quants=re.findall(r"\[\s*[a-zA-Z][\w_]*\s*\]",s)
            for quant in quants:
                if not is_def_sym:
                    check_quant_in_symbols(quant,origin_s)
                tmp=self._wrap_quant(self._check_and_replace_alias(quant).strip("[]").strip())
                s=s.replace(quant,tmp)

            ## Find [sgm*c/E]... form quantiables.
            combined_quants=re.findall(r"\[.*?\]",s)
            for combined_quant in combined_quants:
                new_quants=combined_quant.strip("[]").strip()
                ## Remove the space at the beginning and end of the operator"/","*"
                new_quants=re.sub(r"\s*([/\*])\s*",rm_mul_div_space,new_quants)
                ## Substitute the space with "*"
                new_quants=re.sub(r"\s+","*",new_quants)

                quants=re.findall(r"\b[a-zA-Z][\w_]*\b",new_quants)
                for quant in quants:
                    if self._is_reserved_keywords(quant):
                        continue
                    if not is_def_sym:
                        check_quant_in_symbols("["+quant+"]",origin_s)
                    tmp=self._wrap_quant(self._check_and_replace_alias("["+quant+"]").strip("[]").strip())
                    ## Replace the quantity with the wrapped quantity, be careful about the re.sub method. 
                    ## It is only recommended to be used in the case where the word that is going to be replace is not at two boundaries.
                    new_quants=re.sub(r"\b"+quant+r"\b",tmp,new_quants) 
                new_quants="("+new_quants+")"
                s=s.replace(combined_quant,new_quants)


            ## Check whether brakets match or not.
            bkts=re.findall("[\[\]]",s)
            if bkts:
                raise SquareBracketMismatchError("Brakets does not match for {}".format(origin_s))
            
            ## Check brakets for quantities
            quants=re.findall(r"\b[a-zA-Z][\w_]*\b",s)
            for quant in quants:
                ## There is also a case that the quantity is a reserved keyword, e.g. sqrt, so we need to check it.
                if not self._is_base_unit(quant) and not self._is_reserved_keywords(quant): 
                    raise MissingSquareBracketError("Missing a braket for the quantity {}".format(quant))
        else:
            bkts=re.findall(r"[\[\]]",s)
            if bkts:
                raise QuantityModeError("Square bracktes around {} should be removed in the quantity mode".format(s))
            
            quants=re.findall(r"[a-zA-Z][\w_]*",s)
            for quant in quants:
                if self._is_reserved_keywords(quant):
                    continue
                if not is_def_sym:
                    check_quant_in_symbols("["+quant+"]",origin_s)
                tmp=self._wrap_quant(self._check_and_replace_alias("["+quant+"]").strip("[]").strip())
                # s=s.replace(quant,tmp)
                s = re.sub(r"\b"+quant+r"\b",tmp,s)

            ## Check brakets for quantities
            quants=re.findall(r"\b[a-zA-Z][\w_]*\b",s)
            for quant in quants:
                ## There is also a case that the quantity is a reserved keyword, e.g. sqrt, so we need to check it.
                if not self._is_reserved_keywords(quant): 
                    raise UnknownQuantityError("Unknown quantity {} is found in {}".format(quant,origin_s))

        ## Find all spaces between two quantities and replace them with "*".
        while True:
            new_s=re.sub(r"(_*[a-zA-Z][\w_]*|\d+|\))\s+(_*[a-zA-Z][\w_]*|\()",r"\1*\2",s)
            if new_s==s:
                break
            else:
                s = new_s
        
        # ## Remove all useless space.
        s.replace(" ","")



        
        
        
            
        return s

        

        


    def _quant(self, dim: Quantity, s: str="", is_print=False, is_pretty: bool=False, is_quant: bool=False):
        ## Print the possible quantity by deriving the combination of quantities.
        # quant=self._expr_dict.get(dim)
        quant=self._dim_arr_dict.get(self._generate_dim_array(dim).tobytes())
        if quant:
            if is_quant:
                quant = [q.strip("[]") for q in quant]
            if is_print:
                print("Possible quantity for {} ({}): {}".format(s, dim, self._quant2str(quant, is_inc_alias=True, is_quant=is_quant)))
            return quant
        else:
            if is_pretty:
                dim=self._pretty_sp_expr(dim)
                s=self._pretty_sp_expr(s)
            if is_print:
                print("No proper quantity found for {} ({}).".format(s, dim))
            return None 


    def _append_formulas(self,quant,rhs,discription):
        ## Append the definition of a quantity
        if not self._is_in_existing_symbols(quant):
            self._quant_set.add(quant)
        else:
            raise ConflictSymbolDefinitionError("The symbol {} has already been defined before, please choose another symbol!".format(quant))
        self._lhs_list.append(quant)
        self._rhs_list.append(rhs)
        self._discription_list.append(discription)    

    def _append_alias(self,alias,quant):
        ## Append the alias of a quantity
        self._check_sym_format(alias)
        if not self._is_in_existing_symbols(alias):
            self._alias2quant_map[alias]=quant
        else:
            raise ConflictSymbolDefinitionError("The symbol {} has already been defined before, please choose another symbol!".format(alias))

    def _is_in_existing_symbols(self,quant):
        ## Check whether the quantity is defined or not
        if quant in self._quant_set or quant in self._alias2quant_map:
            return True
        else:
            return False

    def _check_sym_format(self,s):
        ## Check the format of a quantity or alias
        match_obj=re.fullmatch(r"\[[a-zA-Z][\w_]*\]",s)
        if not match_obj:
            raise InvalidSymFormatError("Please use the correct format for the quantity or alias '{}' as the regular expression '\[[a-zA-Z][\w_]\]', e.g. '[sigma_xx]', '[m_11]'".format(s))

 
        










# ## A singleton instance of si unit.
if _DIMCHECK_SI_INSTANCE:
    si=Dimcheck(quant_def_file=_dimcheck_path+os.sep+"si.json",is_save=_DIMCHECK_IS_SAVE, setting_file=_dimcheck_path+os.sep+"setting.json")
else:
    si=None
