import subprocess
import os
import numpy as np
from instr2params import instr2params, param_dict2file
import matplotlib.pyplot as plt
from matplotlib import patches

class McstasSimulation:
    """allows to run a mcstas simulation
    """
    def __init__(self, folder, filename, num_neutrons=10000):
        """initiates a mcstas simulation

        Args:
            num_neutrons (int): number of neutrons
            folder (str): which folder is the simulation in
            filename (str): name of the .instr file
        """
        self.num_neutrons = num_neutrons
        self.gravity = False #gravity is false by definition and can be set true
        self.folder = folder
        self.filename = filename
        self.params_dict = None
        self.new_compile = True


    def init_params_dict(self, file=None):
        """update the params dict of the simulation
        """
        if file==None:
            file = self.filename
        self.params_dict = instr2params(self.folder+self.filename)


    def update_param(self, param, value):
        """updates the params dict given a param and a value

        Args:
            param (str): name of the parameter
            value (float): value of the parameter
        """
        if param in self.params_dict:
            self.params_dict[param] = value
        else:
            print('{} is not a valid parameter, the change will be ignored'.format(param))


    def return_params_dict(self):
        """updates and returns the params_dict

        Returns:
            dict: {param: value} for all mcstas parameters
        """
        if self.params_dict is None:
            self.init_params_dict()
        return self.params_dict


    def return_string_from_dict(self, params_dict=None):
        """returns the command line string for a given params dict

        Args:
            params_dict (dict): dict of parameters and values

        Returns:
            str: command line string
        """
        if params_dict is None:
            if self.params_dict:
                params_dict = self.params_dict
            else:
                params_dict = self.return_params_dict()
        g_string = ' --gravity' if self.gravity else ''
        mcstring = 'mcrun '+ self.folder + self.filename + g_string
        for key in params_dict:
            mcstring += ' '+(key+'='+str(params_dict[key]))
        return mcstring


    def run_simulation(self, mcstring, output=False):
        """runs the mcstas simulation and returns the commandline output

        Args:
            mcstring (str): command string
            compile (bool, optional):
            whether the c-file shall be compiled before, defaults to False
            output (bool, optional):
            whether command line output shall be captured. Defaults to False.

        Returns:
            [type]: [description]
        """
        if self.new_compile:
            mcstring += ' -c -n {}'.format(self.num_neutrons)
        else:
            mcstring += ' -n {}'.format(self.num_neutrons)
        print(mcstring)
        command = mcstring.split()
        sim = subprocess.run(command, cwd=self.folder, capture_output=output)
        return sim


    def return_last_folder(self):
        """returns the last created folder

        Returns:
            str: name of the last created folder
        """
        folders = []
        for root, _, _ in os.walk(self.folder, topdown=False):
            if self.filename.split('.')[0] in root:
                folders += [root]
        return sorted(folders)[-1]


    def delete_last_folder(self):
        """deletes the last created folder
        """
        folder = self.return_last_folder()
        del_string = 'rm -r '+folder
        subprocess.run(del_string.split())


    def return_image_files(self, folder=None):
        """returns all paths to image files in the last folder

        Args:
            folder (str, optional): a path to folder with data. Defaults to None.

        Returns:
            list: list of file strings
        """
        if folder is None:
            folder = self.return_last_folder()
        for _, _, filenames in os.walk(folder):
            return [(folder+'/'+name) for name in filenames if '.dat' in name]

    def compare_image_intensities(self, compare_array, norm_array, folder=None):
        """returns the ratio of intensity for two image keys

        Args:
            compare_file (str): key of the file to compare, e.g. 'psd_sample.dat'
            norm_file (str): key of the norm file, e.g. 'psd_source.dat'
            folder (str, optional): folder of the images. Defaults to None.

        Returns:
            float: ratio of the intensity of the sample file divided by the intensity of the norm
        """
        ratio = np.sum(compare_array)/np.sum(norm_array)
        return ratio


    def return_images_data(self, folder=None) -> dict:
        """returns a dictionary with {path: data}

        Args:
            folder (str, optional): a path to folder with data. Defaults to None.

        Returns:
            dict: dictionary with all image files in the given folder
        """
        data_dict = dict()
        for path in self.return_image_files(folder):
            data_dict[path.split('/')[-1]] = np.genfromtxt(path)
        return data_dict

    def return_images_metadata(self, folder=None) -> dict:
        """returns a dictionary containing meta data for all images in one folder

        Args:
            folder (str, optional): folder in which the images sit. Defaults to None.

        Returns:
            dict: dictionary {'image_name': {parameter: value}}
        """
        meta_dicts = dict()
        for path in self.return_image_files(folder):
            meta_dict = dict()
            with open(path, 'r') as file:
                for line in file:
                    line = line.replace('#', '').replace('\n', '')#.replace(' ', '')
                    if 'Param' in line:
                        key, value = (line.split(':')[-1]).split('=')
                        meta_dict[key] = value
                    else:
                        try:
                            key, value = line.split(':')[-2:]
                            meta_dict[key] = value
                        except ValueError:
                            break
            meta_dicts[path.split('/')[-1]]=meta_dict
        return meta_dicts


    def plot_last_images(self, meta_dicts= None, image_data_dicts=None, folder=None, norm=1.0) -> dict:
        """Plots all images of the last folder

        Args:
            meta_dict (dict, optional): dict comprising the meta data of the images to plot. Defaults to None.
            image_data_dicts (dict, optional): dict comprising the image data of the images to plot. Defaults to None.
            folder (str, optional): folder with the images. Defaults to None.
            norm (float, optional): Weight by which the image is divided. Defaults to 1.

        Returns:
            dict: {file_name: (fig, ax)}
        """
        if meta_dicts == None:
            meta_dicts = self.return_images_metadata(folder)
        if image_data_dicts == None:
            image_data_dicts = self.return_images_data(folder)
        figures_dict = dict()
        for key in meta_dicts.keys():
            data = image_data_dicts[key]
            shape = data.shape
            p, pp, N = data[:shape[1]][::-1]/norm, data[shape[1]:shape[1]*2][::-1]/norm, data[shape[1]*2:][::-1]
            extent = [float(k) for k in meta_dicts[key][' xylimits'].split(' ')[1:]]
            print(extent)
            fig, ax = plt.subplots()
            im = ax.imshow(p, extent=extent, interpolation='none')
            ax.set_xlabel(meta_dicts[key][' xlabel'])
            ax.set_ylabel(meta_dicts[key][' ylabel'])
            ax.set_title(meta_dicts[key][' component'])
            if 'psd' in key:
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            cb = fig.colorbar(im)
            figures_dict[key] = (fig, ax)
        return figures_dict



    def run_and_return(self, params_dict=None):
        """runs a simulation with the current or chosen dict and
        returns a dictionary with the resulting images

        Args:
            params_dict (dict, optional): dict with all parameters. Defaults to None.

        Returns:
            dict: {filename: image_array}
        """
        mcstring = self.return_string_from_dict(params_dict)
        self.run_simulation(mcstring, output=False)
        images_dict = self.return_images_data()
        return images_dict


    def save_params_dict(self, filename='params.txt'):
        """writes the paramsdict parameter to a file in the same folder

        Args:
            filename (str, optional): filename to which to write.
            Defaults to 'params.txt'.
        """
        param_dict2file(self.params_dict, filename=filename)

def read_image_data(array):
    """takes an array and splices it in three parts

    Args:
        array (array): DATA ARRAY OF SHAPE 3*n, n

    Returns:
        tuple: p, pp, N
    """
    height, width = array.shape
    p, pp, N = array[:width], array[width:2*width], array[2*width:]
    return p, pp, N