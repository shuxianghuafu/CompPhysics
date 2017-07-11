# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.37
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        try:
            fp, pathname, description = imp.find_module('_myModule', [dirname(__file__)])
            _mod = imp.load_module('_myModule', fp, pathname, description)
        finally:
            if fp is not None: fp.close()
        return _mod
    _myModule = swig_import_helper()
    del swig_import_helper
else:
    import _myModule
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class Convert_MyArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Convert_MyArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Convert_MyArray, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _myModule.new_Convert_MyArray()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _myModule.delete_Convert_MyArray
    __del__ = lambda self : None;
    def my2py(self, *args): return _myModule.Convert_MyArray_my2py(self, *args)
    def py2my(self, *args): return _myModule.Convert_MyArray_py2my(self, *args)
    def my2py_copy(self, *args): return _myModule.Convert_MyArray_my2py_copy(self, *args)
    def py2my_copy(self, *args): return _myModule.Convert_MyArray_py2my_copy(self, *args)
    __swig_setmethods__["npy_size"] = _myModule.Convert_MyArray_npy_size_set
    __swig_getmethods__["npy_size"] = _myModule.Convert_MyArray_npy_size_get
    if _newclass:npy_size = _swig_property(_myModule.Convert_MyArray_npy_size_get, _myModule.Convert_MyArray_npy_size_set)
    __swig_setmethods__["int_size"] = _myModule.Convert_MyArray_int_size_set
    __swig_getmethods__["int_size"] = _myModule.Convert_MyArray_int_size_get
    if _newclass:int_size = _swig_property(_myModule.Convert_MyArray_int_size_get, _myModule.Convert_MyArray_int_size_set)
    def set_npy_size(self, *args): return _myModule.Convert_MyArray_set_npy_size(self, *args)
    def set_int_size(self, *args): return _myModule.Convert_MyArray_set_int_size(self, *args)
    def dump(self, *args): return _myModule.Convert_MyArray_dump(self, *args)
    def set_pyfunc(self, *args): return _myModule.Convert_MyArray_set_pyfunc(self, *args)
Convert_MyArray_swigregister = _myModule.Convert_MyArray_swigregister
Convert_MyArray_swigregister(Convert_MyArray)

class TestCpp(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TestCpp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TestCpp, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _myModule.new_TestCpp(*args)
        try: self.this.append(this)
        except: self.this = this
    def print_(self): return _myModule.TestCpp_print_(self)
    __swig_destroy__ = _myModule.delete_TestCpp
    __del__ = lambda self : None;
TestCpp_swigregister = _myModule.TestCpp_swigregister
TestCpp_swigregister(TestCpp)


