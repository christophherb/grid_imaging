#reads in an instrument file and returs a params_dic containing all parameters of the simulation
import re
from sys import excepthook
def instr2params(filename):
    """takes a .instr mcStas file and returns a dictionary of default
    parameters {param: value}

    Args:
        filename (str): name of the .instr file 

    Returns:
        dict: {param: default_value, ...}
    """
    note = 1
    params_dic = dict()
    with open(filename) as file:
        for line in file:
            if 'DEFINE' in line:
                note = 1
            if 'DECLARE' in line:
                break
            if note:
                #split the line and fill the values into the dict for 
                #later use
                if '=' in line:
                    param = re.findall('\S* = \S*(?=,)',line)
                    
                    #print(param)
                    try:
                        param = param[0].split(' = ')
                        par,val = param[0],param[1]
                        #print(par,val)
                        try:
                            params_dic[par] = int(val)
                        except ValueError:
                            try:
                                params_dic[par] = float(val)
                            except ValueError:
                                params_dic[par] = str(val)#if nothing else makes sense we just put the string there
                    except IndexError:
                        #print(param)
                        pass
                #print(values.group())
                #scan_params = re.findall('(?<=\s)\w+(?=\s)',line)
    return params_dic

def param_dict2file(param_dict,filename):
    """writes a param_dict to params.txt

    Args:
        param_dict (dict): dictinary containing all param,value pair
    """
    with open(filename, "w") as outF:
        for key in param_dict.keys():
            outF.write(str(key)+' = '+str(param_dict[key])+',')
            outF.write('\n')

