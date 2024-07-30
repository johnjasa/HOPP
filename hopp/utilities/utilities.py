import os
import yaml


class Loader(yaml.SafeLoader):

    def __init__(self, stream):

        self._root = os.path.split(stream.name)[0]

        super().__init__(stream)

    def include(self, node):

        filename = os.path.join(self._root, self.construct_scalar(node))

        with open(filename, 'r') as f:
            return yaml.load(f, self.__class__)


Loader.add_constructor('!include', Loader.include)


class SmartLoader(yaml.Loader):

    def __init__(self, stream):

        self._root = os.path.split(stream.name)[0]

        super().__init__(stream)

    def include(self, node):
        if os.path.split(node.value)[0] == "":
            filename = os.path.join(self._root, self.construct_scalar(node))
        else:
            if os.path.isfile(node.value):
                filename = self.construct_scalar(node)

        with open(filename, 'r') as f:
            return yaml.load(f, self.__class__)


SmartLoader.add_constructor('!include', SmartLoader.include)

class SuperSmartLoader(yaml.Loader):

    def __init__(self, stream):

        self._root = os.path.split(stream.name)[0]

        super().__init__(stream)

    def include(self, node):
        if os.path.split(node.value)[0] == "":
            filename = os.path.join(self._root, self.construct_scalar(node))
        else:
            if os.path.isfile(node.value):
                filename = self.construct_scalar(node)
            else:
                # filename = os.path.join(INPUT_DIR,self.construct_scalar(node))
                if os.path.isfile(os.path.join(INPUT_DIR,self.construct_scalar(node))):
                    filename = os.path.join(INPUT_DIR,self.construct_scalar(node))
        # print(filename)
        with open(filename, 'r') as f:
            return yaml.load(f, self.__class__)


SuperSmartLoader.add_constructor('!include', SuperSmartLoader.include)

def load_yaml(filename, loader=Loader) -> dict:
    if isinstance(filename, dict):
        return filename  # filename already yaml dict
    with open(filename) as fid:
        return yaml.load(fid, loader)

def write_yaml(filename,data):
    if not '.yaml' in filename:
        filename = filename +'.yaml'

    with open(filename, 'w+') as file:
        yaml.dump(data, file,sort_keys=False,encoding = None,default_flow_style=False)
    return filename

def check_create_folder(filepath):
    if not os.path.isdir(filepath):
        os.makedirs(filepath,exist_ok=True)
        already_exists = False
    else:
        already_exists = True

    return already_exists