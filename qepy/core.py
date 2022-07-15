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
    def stdout2file(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            if hasattr(args[0], 'fileobj'):
                fobj = args[0].fileobj
            else :
                fobj = None
            stdout = None
            if fobj is not None :
                if os.fstat(1).st_ino != os.fstat(fobj.fileno()).st_ino :
                    stdout = os.dup(1)
                    os.dup2(fobj.fileno(), 1)
            results = function(*args, **kwargs)
            if stdout is not None :
                os.dup2(stdout, 1)
            return results
        return wrapper

def stdout2file(fileobj = None):
    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            if hasattr(args[0], 'fileobj'):
                fobj = args[0].fileobj
            else :
                fobj = None
            stdout = None
            if fobj is not None :
                if os.fstat(1).st_ino != os.fstat(fobj.fileno()).st_ino :
                    stdout = os.dup(1)
                    os.dup2(fobj.fileno(), 1)
            results = function(*args, **kwargs)
            if stdout is not None :
                os.dup2(stdout, 1)
            return results
        return wrapper
    return decorator
