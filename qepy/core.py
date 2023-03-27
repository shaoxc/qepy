import os
import types
from functools import wraps


class Logger(type):
    def __new__(cls, name, bases, attrs):

        for attr_name, attr_value in attrs.items():
            if isinstance(attr_value, types.FunctionType):
                attrs[attr_name] = cls.stdout2file(attr_value)

        return super().__new__(cls, name, bases, attrs)

    @staticmethod
    def stdout2file(function, fileobj = None):
        return stdout2file(function, fileobj=fileobj)

    @staticmethod
    def print2file(fileobj = None):
        return print2file(fileobj=fileobj)

def print2file(fileobj = None):
    def decorator(function):
        return stdout2file(function, fileobj=fileobj)
    return decorator

def stdout2file(function, fileobj = None):
    @wraps(function)
    def wrapper(*args, **kwargs):
        fobj = fileobj
        if not fobj :
            if len(args) > 0 and hasattr(args[0], 'fileobj'): fobj = args[0].fileobj
        stdout = None
        if fobj is not None and not fobj.closed :
            if os.fstat(1).st_ino != os.fstat(fobj.fileno()).st_ino :
                stdout = os.dup(1)
                os.dup2(fobj.fileno(), 1)
        results = function(*args, **kwargs)
        if stdout is not None :
            os.dup2(stdout, 1)
            os.close(stdout)
        return results
    return wrapper


env = {
    'STDOUT' : None,  # file descriptor of output
    'DRIVER' : None, # save the instance of driver class
}
