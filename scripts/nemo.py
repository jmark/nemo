# Copyright (c) 2019 Johannes Markert <johannes.markert@gmail.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function
import sys,os,re,shutil,filecmp,tempfile,glob
from collections import MutableSet
import fileinput
import fnmatch

# Python 3 compatibility hack
try:
    unicode('')
except NameError:
    unicode = str

def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def remove_directory(path):
    if os.path.exists(path):
        shutil.rmtree(path)

def purge_directory(path):
    remove_directory(path)
    create_directory(path)

def merge_dicts(*dict_args):
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def lazycopy(src,dest):
    """Copy file from 'src' to 'dest' only if they differ."""

    if os.path.exists(dest) and filecmp.cmp(src, dest):
        return False
    shutil.copy(src, dest)
    return True

class LazyFile(object):
    """Append lines to a buffer and write them to 'path' only if
    the newly created file differs from the old file."""

    def __init__(self,path):
        self.path = path
        self.buffer = [] 
        self.ending = "\n"

    def append(self,line=''):
        self.buffer.append(line)
        return line

    def slurp(self,srcpath,transform=lambda x: x):
        with open(srcpath,'r') as fh:
            for line in fh:
                self.append(transform(line.rstrip()))
        return self

    def commit(self):
        with tempfile.NamedTemporaryFile(mode='w+t') as fh:
            fh.writelines(str(line) + self.ending for line in self.buffer)
            fh.flush()
            if os.path.exists(self.path) and filecmp.cmp(self.path, fh.name):
                return None
            shutil.copy(fh.name, self.path)
        return self.path

class OrderedSet(MutableSet):
    """OrderedSet (Python recipe) by Raymond Hettinger
    Activestate Code (http://code.activestate.com/recipes/576694/)"""

    def __init__(self, iterable=None):
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

class Unit(object):
    def __new__(cls,path):
        if cls is Unit:
            if re.search(r'_mod\.f90$',path,flags=re.I):
                return FortranModule(path)
            if re.search(r'_mod\.f03$',path,flags=re.I):
                return FortranModule(path)
            if re.search(r'_mod\.f08$',path,flags=re.I):
                return FortranModule(path)
            if re.search(r'_prog\.f90$',path,flags=re.I):
                return FortranProgram(path)

            raise RuntimeError('Unknown compilation unit: {0}'.format(path))

        return super(Unit,cls).__new__(cls)

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self,other):
        return cmp(self.name,other.name)

    def __eq__(self,other):
        if isinstance(other,Unit):
            return self.name == other.name
        if isinstance(other,str):
            return self.name == other
        raise TypeError('Unknown type in comparison: {}'.format(type(other)))

    def __lt__(self,other):
        return str(self.name) < str(other.name)

    def __gt__(self,other):
        return str(self.name) > str(other.name)

    def __le__(self,other):
        return str(self.name) <= str(other.name)

    def __ge__(self,other):
        return str(self.name) >= str(other.name)

    def __repr__(self):
        return '{0}({1})'.format(self.__class__.__name__,self.name)

class Units(OrderedSet):
    def add(self,unit):
        super(Units,self).discard(unit)
        super(Units,self).add(unit)
        return unit

class FortranModule(Unit):
    regex_name = re.compile(r'^\s*module\s+(\w+_mod)')
    regex_deps = re.compile(r'^\s*use\s+(\w+)\s*,?')

    def __init__(self, path):
        self.path           = path
        self.basename       = os.path.basename(path)
        self.name           = self._scan_module_name(path)
        self.classname      = type(self).__name__.lower()
        self.dependencies   = self._scan_dependencies(path)
        _,self.extension    = os.path.splitext(path)

        if not self.name == os.path.splitext(self.basename)[0].lower():
            raise RuntimeError("Module name does not match basename of the file: '{0}' in '{1}'".format(self.name,self.path))

    @classmethod
    def _scan_module_name(cls,path):
        with open(path,'r') as fh:
            for line in fh:
                match = cls.regex_name.match(line.lower())
                if match:
                    modname = match.group(1)
                    return modname
        raise RuntimeError("Scan for module name failed: {0}".format(path))

    @classmethod
    def _scan_dependencies(cls,path):
        deps = set()
        with open(path,'r') as fh:
            for line in fh:
                match = cls.regex_deps.match(line.lower())
                if match: deps.add(match.group(1))
        return deps

class FortranProgram(FortranModule):
    regex_name = re.compile(r'^\s*program\s+(\w+)')

    @property
    def executable(self):
        "Return the stripped name of the executable."
        return re.sub('_prog$','',self.name)

class Hook(object):
    def __init__(self,build,name):
        self.name = name.lower()
        self.regex_hook = re.compile(r'^\s*subroutine\s+hook_{NAME}_([0-9]+)[\s!]'.format(NAME=self.name))
        self.subscribers = []

        self.filename = 'hook_{NAME}_mod.f90'.format(NAME=self.name)
        self.filepath = os.path.join(build,self.filename)

    def match(self,unit,line):
        match = self.regex_hook.match(line.lower())
        if match:
            self.subscribers.append(dict(unit=unit,priority=match.group(1)))

    def commit(self):
        lf = LazyFile(self.filepath)
        lf.append('module hook_{NAME}_mod'.format(NAME=self.name))
        lf.append()
        lf.append('contains')
        lf.append()
        
        lf.append('subroutine hook_{NAME}'.format(NAME=self.name))
        lf.append()
        for subscriber in sorted(self.subscribers,key=lambda d: int(d['priority'])):
            lf.append('   use {MODULE}, only: hook_{NAME}_{PRIORITY}'.format(
                MODULE=subscriber['unit'].name,NAME=self.name,PRIORITY=subscriber['priority']))
        lf.append()

        for subscriber in sorted(self.subscribers,key=lambda d: int(d['priority'])):
            lf.append('   call hook_{NAME}_{PRIORITY} !! {MODULE}'.format(
                MODULE=subscriber['unit'].name,NAME=self.name,PRIORITY=subscriber['priority']))
        lf.append()

        lf.append('end subroutine')
        lf.append()

        lf.append('end module')

        lf.commit()

def generate_hook_modules(mk):
    regex_hook = re.compile(r'^\s*use\s+hook_([_\w]+)_mod\s*,?')

    hooks = set()
    for unit in mk.units:
        with open(unit.path,'r') as fh:
            for line in fh:
                match = regex_hook.match(line.lower())
                if match: hooks.add(match.group(1))

    hooks = tuple(Hook(mk.build,name) for name in hooks)

    for unit in mk.units:
        with open(unit.path,'r') as fh:
            for line in fh:
                for hook in hooks:
                    hook.match(unit,line)

    for hook in hooks:
        hook.commit()
        mk.add(hook.filepath)

def generate_type_module(mk):
    regex_pld = re.compile(r'^\s*type\s*::\s*payload_(\w+)_t')
    regex_aux = re.compile(r'^\s*type\s*::\s*auxiliary_(\w+)_t')

    filename = 'types_mod.f90'
    filepath = os.path.join(mk.build,filename)

    payload = []
    auxiliary = []

    for unit in mk.units:
        with open(unit.path,'r') as fh:
            for line in fh:
                match = regex_pld.match(line.lower())
                if match: payload.append({'unit': unit, 'name': match.group(1)})
                match = regex_aux.match(line.lower())
                if match: auxiliary.append({'unit': unit, 'name': match.group(1)})

    lf = LazyFile(filepath)
    lf.append('module types_mod')
    lf.append()

    for member in sorted(payload,key=lambda x: x['name']):
        lf.append('use {MODULE}, only: payload_{NAME}_t'.format(
            MODULE=member['unit'].name,NAME=member['name']))
    lf.append()

    for member in sorted(auxiliary,key=lambda x: x['name']):
        lf.append('use {MODULE}, only: auxiliary_{NAME}_t'.format(
            MODULE=member['unit'].name,NAME=member['name']))
    lf.append()

    lf.append('type :: payload_t')
    for member in sorted(payload,key=lambda x: x['name']):
        lf.append('   type(payload_{NAME}_t) :: {NAME}'.format(NAME=member['name']))
    lf.append('end type')
    lf.append()

    lf.append('type :: auxiliary_t')
    for member in sorted(auxiliary,key=lambda x: x['name']):
        lf.append('   type(auxiliary_{NAME}_t) :: {NAME}'.format(NAME=member['name']))
    lf.append('end type')
    lf.append()

    lf.append('type(payload_t) :: payload_dummy')
    lf.append('type(auxiliary_t) :: auxiliary_dummy')
    lf.append()

    lf.append('integer(kind=8), parameter :: payload_size = STORAGE_SIZE(payload_dummy) / 8')
    lf.append('integer(kind=8), parameter :: auxiliary_size = STORAGE_SIZE(auxiliary_dummy) / 8')
    lf.append()

    lf.append('end module')
    lf.commit()

    mk.add(filepath)

def scan_definitions(mk):
    regex_definition = re.compile(r'^\s*!#\s*define\s+(\w+)\s+(.+)')

    defs = dict()
    for unit in mk.units:
        with open(unit.path,'r') as fh:
            for line in fh:
                match = regex_definition.match(line)
                if match: defs[match.group(1).strip()] = match.group(2).strip()

    mk.flags = merge_dicts(defs,mk.flags)

class Makefile(object):
    def __init__(self, root, flags={}, build='build'):
        self.units      = Units()
        self.root       = root
        self.build      = create_directory(build)
        self.flags      = flags
        self.program    = self.add(os.path.join(root,'source','nemo_prog.f90'))

    @staticmethod
    def _glob(pathexpr):
        # recursive search
        if pathexpr.find('**') > 0:
            matches = []
            for root,dirnames,filenames in os.walk(os.path.dirname(pathexpr)):
                for filename in fnmatch.filter(filenames,os.path.basename(pathexpr)):
                    matches.append(os.path.join(root,filename))
        # current level search
        else:
            matches = glob.glob(pathexpr)

        return list(map(os.path.abspath,matches))

    def add(self,path,alias=None):
        units = []
        paths = Makefile._glob(path)

        if len(paths) < 1:
            raise RuntimeError('No modules found! Given path: {0}'.format(path))
        if len(paths) > 1 and not alias is None:
            raise RuntimeError('Multiple modules found but \'alias\' given! {0}'.format(path))

        if alias is None:
            for p in paths:
                units.append(self.units.add(Unit(p)))
            return units[0] if len(units) == 1 else units
        else:
            unit = Unit(paths[0])
            path = os.path.abspath(os.path.join(create_directory(os.path.join(self.build,'unfiltered')),alias+unit.extension))
            LazyFile(path).slurp(unit.path,transform=lambda line: line.replace(unit.name,alias)).commit()
            return self.units.add(Unit(path))

    def discard(self,unit):
        if isinstance(unit,str):
            self.units.discard(unit)
        else:
            raise RuntimeError('Unknown unit type: {}'.formrat(unit))

    def _is_internal_dependency(self,name):
        return any(unit.name == name for unit in self.units)
   
    def _generate_makefile(self):
        lf = LazyFile(os.path.join(self.build, 'Makefile'))

        ## Header ##

        lf.append('all: filter compile')
        lf.append()

        ## Compiler & FLAGS ##

        for line in fileinput.input(['-']):
            lf.append(str(line).rstrip())
        lf.append()

        lf.append('FLAGS=')
        for key in sorted(self.flags):
            lf.append("FLAGS += -D{KEY}='{VALUE}'".format(KEY=key,VALUE=self.flags[key]))
        lf.append()

        ## Source Code Filtering ##
        def _quote(value):
            value = unicode(str(value))
            if value.isnumeric():
                return value
            return '"{0}"'.format(value)

        lf.append('PP = {0}/scripts/fypp -F -m time -m math {1}'.format(self.root,' '.join(
            "-D{0}='{1}'".format(key,_quote(self.flags[key])) for key in sorted(self.flags)
        )))
        lf.append()

        for unit in sorted(self.units):
            lf.append('{TARGET}: {DEPENDENCY}'.format(TARGET=unit.basename,DEPENDENCY='Makefile'))
            lf.append('{TARGET}: {DEPENDENCY}'.format(TARGET=unit.basename,DEPENDENCY=unit.path))
            lf.append("\t$(PP) {0} $< $@".format("-D{0}='\"{1}\"'".format('PP_MODULE',unit.name)))
            lf.append()

        for unit in sorted(self.units):
            lf.append('filter: {DEPENDENCY}'.format(DEPENDENCY=unit.basename))
        lf.append()

        ## Compilation ##

        for unit in sorted(self.units):
            for dep in sorted(filter(self._is_internal_dependency,unit.dependencies)):
                lf.append('{TARGET}.o: {DEPENDENCY}.o'.format(TARGET=unit.name,DEPENDENCY=dep))

            lf.append('{TARGET}.o: {DEPENDENCY}'.format(TARGET=unit.name,DEPENDENCY=unit.basename))
            lf.append("\t$(FC) -c -o $@ $< $(FCFLAGS) $(FLAGS)")
            lf.append()

        ## Linking ##

        for unit in sorted(self.units):
            lf.append('{TARGET}: {DEPENDENCY}.o'.format(TARGET=self.program.executable,DEPENDENCY=unit.name))
        lf.append("\t$(FC) -o $@ $^ $(FCFLAGS) $(FLAGS)")
        lf.append("\t@echo")
        lf.append("\t@echo '!! --- [SUCCESS] --- !!'")
        lf.append("\t@echo")
        lf.append()

        ## Footer ##

        lf.append('compile: {PROGRAM}'.format(PROGRAM=self.program.executable))
        lf.append()
        lf.append('.PHONY: all filter compile')

        lf.commit()

    def generate(self):
        scan_definitions(self)
        generate_type_module(self)
        generate_hook_modules(self)
        self._generate_makefile()
